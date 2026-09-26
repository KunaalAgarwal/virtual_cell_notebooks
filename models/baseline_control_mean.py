"""Control-mean baseline: predict "no change" for every perturbation.

For each context, the prediction for all 300 targets is the per-gene mean of that
context's non-targeting control cells, rounded to integer counts.

The output matrix is ~2.6e9 nonzeros (~21 GB as in-memory CSR), so both the mean
computation and the write stream in chunks and never materialize it.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem

CONTEXTS = ("A", "B", "C")
STREAM_ELEMS = 1 << 24  # elements per read chunk when accumulating the mean
WRITE_ROWS = 2000  # rows per write chunk
H5_CHUNK = 1 << 20  # HDF5 chunk size; large so gzip sees many repeats of the row pattern


def _read_single_column_csv(path: Path) -> list[str]:
    return pd.read_csv(path).iloc[:, 0].astype(str).tolist()


def _control_mean_counts(h5ad_path: Path, n_genes: int) -> np.ndarray:
    """Per-gene mean over the control cells, rounded to integer counts.

    Streams X/data and X/indices instead of loading the matrix (each context is
    ~110M nonzeros).
    """
    with h5py.File(h5ad_path, "r") as f:
        shape = tuple(f["X"].attrs["shape"])
        if shape[1] != n_genes:
            raise ValueError(f"{h5ad_path.name}: has {shape[1]} genes, gene list has {n_genes}")
        data, indices = f["X/data"], f["X/indices"]
        nnz = data.shape[0]
        sums = np.zeros(n_genes, dtype=np.float64)
        for start in range(0, nnz, STREAM_ELEMS):
            stop = min(start + STREAM_ELEMS, nnz)
            sums += np.bincount(
                indices[start:stop].astype(np.intp),
                weights=data[start:stop],
                minlength=n_genes,
            )
        return np.rint(sums / shape[0]).astype(np.int32)


def _check_gene_order(h5ad_path: Path, gene_names: list[str]) -> None:
    with h5py.File(h5ad_path, "r") as f:
        var = f["var/_index"]
        raw = var["values"][:] if isinstance(var, h5py.Group) else var[:]
    file_genes = [g.decode() if isinstance(g, bytes) else str(g) for g in raw]
    if file_genes != gene_names:
        raise ValueError(
            f"{h5ad_path.name}: gene order differs from the gene list "
            f"(first mismatch at {next(i for i, (a, b) in enumerate(zip(file_genes, gene_names)) if a != b)})"
        )


def predict(context_paths, gene_names_path, perts_path, out_path, cells_per_pert=400, perts_limit=None):
    """Write a control-mean prediction .h5ad covering every target in every context.

    context_paths: mapping of context label -> control .h5ad path
    """
    t0 = time.perf_counter()
    context_paths = {k: Path(v) for k, v in dict(context_paths).items()}
    gene_names = _read_single_column_csv(Path(gene_names_path))
    perts = _read_single_column_csv(Path(perts_path))
    if perts_limit:
        perts = perts[:perts_limit]
    n_genes = len(gene_names)
    out_path = Path(out_path)

    rows = {}
    for ctx, path in context_paths.items():
        _check_gene_order(path, gene_names)
        mean_counts = _control_mean_counts(path, n_genes)
        idx = np.flatnonzero(mean_counts).astype(np.int32)
        rows[ctx] = (idx, mean_counts[idx].astype(np.float32))
        print(
            f"context {ctx}: {idx.size} nonzero genes, library size {int(mean_counts.sum())}",
            flush=True,
        )

    rows_per_ctx = len(perts) * cells_per_pert
    n_cells = rows_per_ctx * len(context_paths)
    nnz_per_ctx = {c: int(rows[c][0].size) for c in context_paths}
    total_nnz = sum(n * rows_per_ctx for n in nnz_per_ctx.values())
    print(f"writing {n_cells} cells x {n_genes} genes, {total_nnz} nonzeros -> {out_path}", flush=True)

    obs = pd.DataFrame(
        {
            "target_gene": pd.Categorical(
                np.tile(np.repeat(np.array(perts, dtype=object), cells_per_pert), len(context_paths))
            ),
            "context": pd.Categorical(
                np.repeat(np.array(list(context_paths), dtype=object), rows_per_ctx)
            ),
        },
        index=[
            f"{ctx}_{pert}_{i:03d}"
            for ctx in context_paths
            for pert in perts
            for i in range(cells_per_pert)
        ],
    )
    var = pd.DataFrame(index=pd.Index(gene_names))

    # indptr must be int64: total_nnz exceeds the int32 range.
    indptr = np.zeros(n_cells + 1, dtype=np.int64)
    per_row = np.concatenate([np.full(rows_per_ctx, nnz_per_ctx[c], dtype=np.int64) for c in context_paths])
    np.cumsum(per_row, out=indptr[1:])

    with h5py.File(out_path, "w") as f:
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        write_elem(f, "obs", obs)
        write_elem(f, "var", var)
        for key in ("layers", "obsm", "obsp", "varm", "varp", "uns"):
            write_elem(f, key, {})

        x = f.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["encoding-version"] = "0.1.0"
        x.attrs["shape"] = np.array([n_cells, n_genes], dtype=np.int64)
        x.create_dataset("indptr", data=indptr, compression="gzip", compression_opts=1)
        data_ds = x.create_dataset(
            "data", shape=(total_nnz,), dtype=np.float32,
            chunks=(H5_CHUNK,), compression="gzip", compression_opts=1,
        )
        idx_ds = x.create_dataset(
            "indices", shape=(total_nnz,), dtype=np.int32,
            chunks=(H5_CHUNK,), compression="gzip", compression_opts=1,
        )

        cursor = 0
        for ctx in context_paths:
            idx, vals = rows[ctx]
            width = idx.size
            for start in range(0, rows_per_ctx, WRITE_ROWS):
                n_rows = min(WRITE_ROWS, rows_per_ctx - start)
                span = n_rows * width
                data_ds[cursor:cursor + span] = np.tile(vals, n_rows)
                idx_ds[cursor:cursor + span] = np.tile(idx, n_rows)
                cursor += span
            print(f"  context {ctx} written ({cursor}/{total_nnz} nonzeros)", flush=True)

    elapsed = time.perf_counter() - t0
    size_gb = out_path.stat().st_size / 1e9
    print(f"done in {elapsed:.1f}s | {out_path} is {size_gb:.3f} GB", flush=True)
    return out_path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--controls-dir", type=Path, default=Path("data/validation_2026"),
                   help="directory holding context_{A,B,C}.h5ad, gene_names.csv, pert_counts.csv")
    p.add_argument("-g", "--genes", type=Path, default=None, help="gene_names.csv (default: <controls-dir>/gene_names.csv)")
    p.add_argument("--perts", type=Path, default=None, help="pert_counts.csv (default: <controls-dir>/pert_counts.csv)")
    p.add_argument("-o", "--output", type=Path, default=Path("submissions/pred.h5ad"))
    p.add_argument("--cells-per-pert", type=int, default=400)
    p.add_argument("--limit-perts", type=int, default=None, help="use only the first N targets (smoke tests)")
    args = p.parse_args()

    context_paths = {c: args.controls_dir / f"context_{c}.h5ad" for c in CONTEXTS}
    missing = [str(v) for v in context_paths.values() if not v.exists()]
    if missing:
        p.error("missing control files: " + ", ".join(missing) + "\nRun: vcc datasets download controls")

    predict(
        context_paths,
        args.genes or args.controls_dir / "gene_names.csv",
        args.perts or args.controls_dir / "pert_counts.csv",
        args.output,
        cells_per_pert=args.cells_per_pert,
        perts_limit=args.limit_perts,
    )


if __name__ == "__main__":
    main()
