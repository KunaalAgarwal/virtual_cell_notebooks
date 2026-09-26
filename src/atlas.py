"""Anonymous, cached access to the Arc Virtual Cell Atlas.

Covers all three datasets in the bucket: the 2025 VCC set, scBaseCount, and Tahoe-100M.
Reads are free -- the bucket is public, so there is no GCP project, credential, or billing
account involved. What limits you is local disk and bandwidth, not cost.

Downloads are cached under data/atlas_cache/ mirroring the remote layout. A file already
present is not re-fetched: by default its size is checked against the remote (one cheap
metadata call), and --no-verify skips even that so a warm cache touches the network zero times.

CLI (run from the repo root):
    python src/atlas.py ls vcc2025
    python src/atlas.py ls scbasecount h5ad/GeneFull_Ex50pAS/Homo_sapiens
    python src/atlas.py sizes tahoe100M metadata
    python src/atlas.py get tahoe100M metadata/sample_metadata.parquet
    python src/atlas.py get vcc2025 train/adata_Training.h5ad --dry-run

Library:
    import sys; sys.path.insert(0, 'src')      # or '../src' from notebooks/
    from atlas import fetch, open_remote, ls
    p = fetch("vcc2025", "train/pert_counts_Training.csv")   # local Path, cached
    with open_remote("vcc2025", "train/adata_Training.h5ad") as fo:  # no download
        ...
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import gcsfs

BUCKET = "arc-institute-virtual-cell-atlas"

# Roots for each dataset. Pinned to specific snapshot dates so runs are reproducible;
# bump deliberately rather than tracking "latest".
DATASETS = {
    "vcc2025": f"{BUCKET}/virtual-cell-challenge/2025",
    "scbasecount": f"{BUCKET}/scbasecount/2026-01-12",
    "tahoe100M": f"{BUCKET}/tahoe100M/2025-02-25",
}

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = REPO_ROOT / "data" / "atlas_cache"

_fs = None

# Above this many categories, read only the ones a slice references.
_CATEGORY_FULL_READ_MAX = 100_000


def fs() -> gcsfs.GCSFileSystem:
    """Anonymous filesystem handle. token='anon' is what avoids needing credentials."""
    global _fs
    if _fs is None:
        _fs = gcsfs.GCSFileSystem(token="anon")
    return _fs


def remote(dataset: str, subpath: str = "") -> str:
    if dataset not in DATASETS:
        raise KeyError(f"unknown dataset {dataset!r}; choose from {sorted(DATASETS)}")
    root = DATASETS[dataset]
    return f"{root}/{subpath.strip('/')}" if subpath else root


def ls(dataset: str, subpath: str = "", detail: bool = False):
    return fs().ls(remote(dataset, subpath), detail=detail)


def sizes(dataset: str, subpath: str = "") -> list[dict]:
    """Entries under a prefix with sizes, largest first. Metadata only, no transfer."""
    out = []
    for it in fs().ls(remote(dataset, subpath), detail=True):
        out.append(
            {
                "name": it["name"].rsplit("/", 1)[-1] or it["name"].rstrip("/").rsplit("/", 1)[-1],
                "is_dir": it["type"] == "directory",
                "bytes": int(it.get("size") or 0),
            }
        )
    return sorted(out, key=lambda r: -r["bytes"])


def local_path(dataset: str, subpath: str) -> Path:
    return CACHE_DIR / dataset / subpath.strip("/")


def open_remote(dataset: str, subpath: str, mode: str = "rb"):
    """Open a remote file without downloading it.

    Combined with h5py this gives range-request access, so you can inspect a 25 GB
    file in seconds and read only the slices you need.
    """
    return fs().open(remote(dataset, subpath), mode)


def peek(dataset: str, subpath: str) -> dict:
    """Shape/format/columns of a remote .h5ad without downloading it."""
    import h5py

    with open_remote(dataset, subpath) as fo:
        f = h5py.File(fo, "r")
        n_cells, n_genes = (int(v) for v in f["X"].attrs["shape"])
        nnz = int(f["X/data"].shape[0])
        return {
            "cells": n_cells,
            "genes": n_genes,
            "nnz": nnz,
            "density": round(nnz / (n_cells * n_genes), 3),
            "format": f["X"].attrs["encoding-type"],
            "dtype": str(f["X/data"].dtype),
            "obs_cols": list(f["obs"].attrs.get("column-order", [])),
            "var_index": f["var"].attrs.get("_index", "_index"),
            "var_cols": [k for k in f["var"].keys()],
            "layers": list(f["layers"].keys()) if "layers" in f else [],
        }


def _decode(arr):
    if getattr(arr, "dtype", None) is not None and arr.dtype.kind in ("S", "O"):
        return [x.decode() if isinstance(x, bytes) else x for x in arr]
    return arr


def _read_col_slice(elem, start, stop):
    """Read rows [start:stop) of one encoded obs/var column."""
    import h5py
    import pandas as pd

    if isinstance(elem, h5py.Group):
        enc = elem.attrs.get("encoding-type")
        if enc == "categorical":
            import numpy as np

            codes = elem["codes"][start:stop]
            cat_ds = elem["categories"]
            ordered = bool(elem.attrs.get("ordered", False))
            # Reading every category is fine when there are few, but some columns are
            # near-unique per cell (Tahoe's BARCODE has 3M categories), so for those
            # fetch only the ones this slice actually references.
            if cat_ds.shape[0] > _CATEGORY_FULL_READ_MAX:
                used = np.unique(codes[codes >= 0])
                cats = _decode(cat_ds[used]) if used.size else []
                remap = np.full(int(cat_ds.shape[0]), -1, dtype=np.int64)
                remap[used] = np.arange(used.size)
                codes = np.where(codes >= 0, remap[codes], -1)
            else:
                cats = _decode(cat_ds[:])
            return pd.Categorical.from_codes(codes, categories=cats, ordered=ordered)
        values = _decode(elem["values"][start:stop])
        series = pd.Series(values)
        if "mask" in elem:
            series[elem["mask"][start:stop]] = pd.NA
        return series.to_numpy()
    return _decode(elem[start:stop])


def read_obs_slice(group, start: int, stop: int):
    """Rows [start:stop) of an encoded h5ad dataframe group.

    Slices each column at the dataset level. `read_elem` would materialize the whole
    table first, which is fatal on Tahoe (obs is 8M rows x 16 columns).
    """
    import pandas as pd

    cols = [c for c in group.attrs.get("column-order", [])]
    data = {c: _read_col_slice(group[c], start, stop) for c in cols}
    index_key = group.attrs.get("_index", "_index")
    index = _read_col_slice(group[index_key], start, stop)
    return pd.DataFrame(data, index=pd.Index(index))


def slice_cells(dataset: str, subpath: str, start: int, stop: int):
    """AnnData for rows [start:stop) of a remote CSR .h5ad, without downloading it.

    CSR only. A csc_matrix file is column-major, so pulling a few cells would touch
    every column block -- for those, fetch() the whole file instead (scBaseCount
    per-accession files are small enough that this is fine).
    """
    import anndata as ad
    import h5py
    import scipy.sparse as sp
    from anndata.io import read_elem

    with open_remote(dataset, subpath) as fo:
        f = h5py.File(fo, "r")
        fmt = f["X"].attrs["encoding-type"]
        if fmt != "csr_matrix":
            raise ValueError(f"{subpath} is {fmt}; slice_cells needs csr_matrix -- use fetch() instead")
        n_genes = int(f["X"].attrs["shape"][1])
        indptr = f["X/indptr"][start : stop + 1]
        lo, hi = int(indptr[0]), int(indptr[-1])
        X = sp.csr_matrix(
            (f["X/data"][lo:hi], f["X/indices"][lo:hi], indptr - indptr[0]),
            shape=(stop - start, n_genes),
        )
        obs = read_obs_slice(f["obs"], start, stop)
        var = read_elem(f["var"])
    return ad.AnnData(X=X, obs=obs, var=var)


def fetch(
    dataset: str,
    subpath: str,
    force: bool = False,
    verify: bool = True,
    min_free_ratio: float = 1.15,
) -> Path:
    """Download to the cache and return the local path; a no-op if already cached.

    verify=True compares the cached size against the remote (one metadata call).
    verify=False trusts any existing file and makes no network call at all.
    """
    dest = local_path(dataset, subpath)

    if dest.exists() and not force:
        if not verify:
            return dest
        remote_size = int(fs().info(remote(dataset, subpath))["size"])
        if dest.stat().st_size == remote_size:
            return dest
        print(f"cached {dest.name} is {dest.stat().st_size} B, remote is {remote_size} B -- refetching")

    remote_size = int(fs().info(remote(dataset, subpath))["size"])
    dest.parent.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(dest.parent).free
    if free < remote_size * min_free_ratio:
        raise OSError(
            f"need ~{remote_size / 1e9:.1f} GB (plus headroom) for {subpath}, "
            f"only {free / 1e9:.1f} GB free"
        )

    # Download to .part first so an interrupted transfer never looks like a complete file.
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"fetching {subpath} ({remote_size / 1e9:.3f} GB) -> {dest}")
    fs().get(remote(dataset, subpath), str(tmp))
    got = tmp.stat().st_size
    if got != remote_size:
        tmp.unlink(missing_ok=True)
        raise IOError(f"size mismatch for {subpath}: got {got} B, expected {remote_size} B")
    tmp.replace(dest)
    return dest


def _fmt(n: int) -> str:
    return f"{n / 1e9:.3f} GB" if n >= 1e9 else f"{n / 1e6:.2f} MB"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    for name in ("ls", "sizes"):
        s = sub.add_parser(name)
        s.add_argument("dataset", choices=sorted(DATASETS))
        s.add_argument("subpath", nargs="?", default="")

    g = sub.add_parser("get")
    g.add_argument("dataset", choices=sorted(DATASETS))
    g.add_argument("subpaths", nargs="+")
    g.add_argument("--force", action="store_true", help="refetch even if cached")
    g.add_argument("--no-verify", action="store_true", help="trust the cache; make no network call")
    g.add_argument("--dry-run", action="store_true", help="report sizes and disk headroom, download nothing")

    a = p.parse_args()

    if a.cmd in ("ls", "sizes"):
        rows = sizes(a.dataset, a.subpath)
        total = sum(r["bytes"] for r in rows)
        n_dir = sum(r["is_dir"] for r in rows)
        for r in rows[:40]:
            tag = "d" if r["is_dir"] else "f"
            print(f"  [{tag}] {r['name'][:60]:60s} {_fmt(r['bytes']) if r['bytes'] else ''}")
        if len(rows) > 40:
            print(f"  ... +{len(rows) - 40} more")
        print(f"\n{len(rows) - n_dir} files, {n_dir} dirs, total {_fmt(total)}")
        return

    if a.dry_run:
        total = 0
        for sp in a.subpaths:
            n = int(fs().info(remote(a.dataset, sp))["size"])
            cached = local_path(a.dataset, sp).exists()
            print(f"  {sp}  {_fmt(n)}{'  [cached]' if cached else ''}")
            if not cached:
                total += n
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        free = shutil.disk_usage(CACHE_DIR).free
        print(f"\nwould download {_fmt(total)}; {_fmt(free)} free")
        return

    for sp in a.subpaths:
        out = fetch(a.dataset, sp, force=a.force, verify=not a.no_verify)
        print(f"  ready: {out}  ({_fmt(out.stat().st_size)})")


if __name__ == "__main__":
    main()
