"""Validate a combined prediction with the official vcc CLI and package it as a .vcc.

`vcc prep` holds the whole matrix in memory (about 1.3x of ~16 GB for a full 2026
prediction), which does not fit on a 16 GB machine. This script gets the same guarantees
in pieces:

  1. each context is cut out to its own .h5ad and checked by the official
     `vcc prep --dry-run --contexts <X>` (genes, targets, 400 cells each, integer counts,
     no controls, per-cell cap);
  2. the combined obs is checked by the CLI's own validate_contexts / validate_targets;
  3. the .vcc is written exactly as vcc's _write_vcc does: tar(meta.json, pred.h5ad.zst),
     with the minimal obs (target_gene, context; index "0".."n-1") that prep would write.

On a machine with >= 32 GB RAM, `vcc prep pred.h5ad -g gene_names.csv --perts
pert_counts.csv -o pred.vcc` on the same .h5ad is the one-step equivalent.

Run with the vcc tool's own Python so the checks are the CLI's code:
    ~/.local/share/uv/tools/vcc-cli/bin/python package_vcc.py --pred pred.h5ad \
        --data vcc2026-data/ --out pred.vcc
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import zstandard as zstd
from vcc import __version__ as vcc_version
from vcc import prep as vprep
from vcc.vccfile import META_MEMBER, META_SCHEMA, PRED_MEMBER, validate_vcc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--vcc-bin", default=shutil.which("vcc") or os.path.expanduser("~/.local/bin/vcc"))
    ap.add_argument("--tmp", default=None)
    a = ap.parse_args()
    data = Path(a.data)
    genes_csv, perts_csv = str(data / "gene_names.csv"), str(data / "pert_counts.csv")

    # 0. minimal obs, as prep writes it
    with h5py.File(a.pred, "r") as f:
        obs = ad.io.read_elem(f["obs"])
        nnz = int(f["X/indptr"][-1])
        n_obs, n_vars = (int(v) for v in f["X"].attrs["shape"])
    minimal_obs = pd.DataFrame({"target_gene": obs["target_gene"].astype(str).to_numpy(),
                                "context": obs["context"].astype(str).to_numpy()},
                               index=np.arange(n_obs).astype(str))

    # 1. official dry run per context
    contexts = sorted(minimal_obs["context"].unique())
    backed = ad.read_h5ad(a.pred, backed="r")
    for ctx in contexts:
        rows = np.flatnonzero(minimal_obs["context"].to_numpy() == ctx)
        assert np.all(np.diff(rows) == 1), "contexts must be contiguous"
        with tempfile.TemporaryDirectory(dir=a.tmp) as tmp:
            part = backed[rows[0]:rows[-1] + 1].to_memory()
            assert part.X.has_sorted_indices
            p = os.path.join(tmp, f"context_{ctx}.h5ad")
            part.write_h5ad(p)
            del part
            r = subprocess.run([a.vcc_bin, "prep", p, "-g", genes_csv, "--perts", perts_csv,
                                "--contexts", ctx, "--dry-run", "--json"], capture_output=True, text=True)
            print(f"context {ctx}: vcc prep --dry-run exit {r.returncode}\n{r.stdout.strip()}\n{r.stderr.strip()[-2000:]}",
                  flush=True)
            if r.returncode != 0:
                raise SystemExit(f"context {ctx} failed vcc prep")

    # 2. official obs-level checks on the combined file
    obs_only = ad.AnnData(obs=minimal_obs, var=pd.DataFrame(index=backed.var_names))
    per_ctx = vprep.validate_contexts(obs_only, "context", "target_gene", "non-targeting", ("A", "B", "C"))
    vprep.validate_targets(obs_only, vprep.read_pert_counts(perts_csv, "target_gene", "context", "non-targeting",
                                                            default_n_cells=400),
                           "context", "target_gene", "non-targeting", ("A", "B", "C"), check_cell_counts=True)
    assert backed.var_names.tolist() == vprep.read_gene_list(genes_csv)
    assert n_obs <= 400_000 and nnz <= 4_750_000_000
    print("combined: contexts", per_ctx, "targets and cell counts verified; nnz", f"{nnz:,}", flush=True)
    backed.file.close()

    # 3. package
    # the compressed payload goes next to the output (a local disk may not hold both copies)
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.out))) as tmp:
        zst = os.path.join(tmp, "pred.h5ad.zst")
        # rewrite obs to the minimal form in a copy-free way: obs is tiny, X is untouched
        with h5py.File(a.pred, "r+") as f:
            del f["obs"]
            ad.io.write_elem(f, "obs", minimal_obs)
        cctx = zstd.ZstdCompressor(level=3, threads=max(1, len(os.sched_getaffinity(0))))
        with open(a.pred, "rb") as src, open(zst, "wb") as dst:
            cctx.copy_stream(src, dst)
        meta = os.path.join(tmp, META_MEMBER)
        with open(meta, "w", encoding="utf-8") as fh:
            json.dump({"schema": META_SCHEMA, "nnz": nnz, "n_obs": n_obs, "n_vars": n_vars,
                       "cli_version": vcc_version}, fh, sort_keys=True)

        def _normalize(ti):
            ti.uid = ti.gid = 0; ti.uname = ti.gname = ""; ti.mtime = 0
            return ti
        with tarfile.open(a.out, "w") as tar:
            tar.add(meta, arcname=META_MEMBER, filter=_normalize)
            tar.add(zst, arcname=PRED_MEMBER, filter=_normalize)
    validate_vcc(a.out)
    print(f"wrote {a.out} ({os.path.getsize(a.out) / 2**30:.2f} GiB); validate_vcc passed")


if __name__ == "__main__":
    main()
