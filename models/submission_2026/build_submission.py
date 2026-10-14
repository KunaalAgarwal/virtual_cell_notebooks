"""Build a Virtual Cell Challenge 2026 submission from the public-screen effect-borrowing method.

For each context (A, B, C) the method of public-screens/RESULTS.md (Round 4) predicts a log
change for every one of the 300 targets:

  * shared response: H1's mean knockdown effect (Nadig et al. 2025 / VCC 2025 H1 data),
  * specific part: H1 effects of the m genes most similar to the target in the Replogle 2022
    K562 genome-wide screen, times an amplitude C (c8_strength: C_k = 8 * s_k^2 with s_k the
    K562 effect strength, m = 10; c1_calibrated, the default: C = 1, m = 20),
  * the target's own K562 profile converted to an approximate log fold change (weight 1),
  * no change on genes detected in < 1% of this context's control cells,
  * the target gene itself drops to its K562 remaining fraction.

Cells are made in counts space, as the 2026 scorer requires: each predicted cell is a real
control cell from the same context whose counts are scaled by exp(predicted log change) gene by
gene, with integer counts kept by binomial thinning (fold < 1) or adding Poisson counts
(fold > 1). Every target gets exactly 400 cells per context; no control cells are written.

The combined .h5ad is written in a streaming way (CSR appended row block by row block), so the
whole matrix (~2e9 stored counts, ~17 GB) never sits in memory.

    python build_submission.py --data vcc2026-data/ --h1 h1_pseudobulk.npz \
        --k562 K562_gwps_normalized_bulk_01.h5ad --method-dir public-screens/ --out pred.h5ad [--preset c1_calibrated]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

CTRL = "non-targeting"
CELLS_PER_PERT = 400
PRESETS = {
    # Round 4 recommended default (tuned for cell-eval 0.8 discrimination on held-out lines)
    "c8_strength": dict(sim="corr", m=10, alpha=0.5, C=8.0, shared="h1", amp="strength", gamma=2.0,
                        gate=0.01, w_direct=1.0, self_drop="k562"),
    # Round 4 "tuned for Pearson" (best MAE on both held-out lines). Its effect sizes match real
    # H1 knockdowns (median norm 2.0 vs 2.0), where c8_strength overshoots 3x (15x at p90);
    # the 2026 metrics are scale-free for discrimination and penalize overshoot in mse / lfc_nmae.
    "c1_calibrated": dict(sim="corr", m=20, alpha=0.5, C=1.0, shared="h1", amp="fixed", gamma=0.5,
                          gate=0.01, w_direct=1.0, self_drop="k562"),
}
PARAMS = PRESETS["c8_strength"]
FOLD_CLIP = (np.exp(-3.0), np.exp(3.0))
SELF_FLOOR = 0.05


def context_stats(X):
    """Control mean of log1p(CP10k) per gene and detection rate, from a raw-count CSR."""
    lib = np.asarray(X.sum(1)).ravel()
    L = X.multiply(1e4 / lib[:, None]).tocsr()
    L.data = np.log1p(L.data)
    return np.asarray(L.mean(0)).ravel(), np.asarray((X > 0).mean(0)).ravel()


def predict_folds(ctx, X, genes, targets, h1, k562, xb):
    """Per-target fold factors (targets x genes) on the count scale, plus the log changes."""
    ctrl_mean, det = context_stats(X)
    line = dict(name=ctx, genes=np.asarray(genes), perts=list(targets), ctrl_mean=ctrl_mean, det=det)
    F = xb.build_features(line, h1, k562)
    delta = xb.predict(line, F, PARAMS)
    fold = np.clip(np.exp(delta), *FOLD_CLIP)
    gi = {g: i for i, g in enumerate(genes)}
    fe = F["meta"]["fold_expr"].to_numpy()
    for i, t in enumerate(targets):  # the target's own drop, on the count scale
        fold[i, gi[t]] = np.clip(fe[i], SELF_FLOOR, 1.0) if np.isfinite(fe[i]) else np.exp(xb.H1_SELF_DROP)
    info = pd.DataFrame({"target": targets, "in_k562": F["in_k"], "k562_fold_expr": fe,
                         "n_genes_changed_10pct": (np.abs(np.log(fold)) > np.log(1.1)).sum(1)})
    return fold.astype(np.float32), delta, info


def make_cells(X, fold_row, rng):
    """400 control cells, counts scaled by fold_row with integer-preserving sampling."""
    idx = rng.choice(X.shape[0], CELLS_PER_PERT, replace=False)
    S = X[idx].tocsr()
    c = S.data.astype(np.int64)
    f = fold_row[S.indices]
    down = f < 1
    out = c.copy()
    out[down] = rng.binomial(c[down], f[down])
    up = ~down
    out[up] = c[up] + rng.poisson(c[up] * (f[up] - 1.0))
    S.data = out.astype(np.float32)
    S.eliminate_zeros()
    return S


class CSRWriter:
    """Append CSR row blocks to an h5ad X group without holding the matrix in memory."""

    def __init__(self, f, n_genes):
        g = f.create_group("X")
        g.attrs["encoding-type"] = "csr_matrix"
        g.attrs["encoding-version"] = "0.1.0"
        self.g, self.n_genes, self.nnz, self.rows = g, n_genes, 0, 0
        g.create_dataset("data", shape=(0,), maxshape=(None,), dtype="float32", chunks=(1 << 20,))
        g.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32", chunks=(1 << 20,))
        self.indptr = [np.zeros(1, dtype=np.int64)]

    def append(self, S):
        n = S.nnz
        for name, arr in (("data", S.data.astype(np.float32)), ("indices", S.indices.astype(np.int32))):
            d = self.g[name]
            d.resize((self.nnz + n,))
            d[self.nnz:self.nnz + n] = arr
        self.indptr.append(S.indptr[1:].astype(np.int64) + self.nnz)
        self.nnz += n
        self.rows += S.shape[0]

    def close(self):
        self.g.create_dataset("indptr", data=np.concatenate(self.indptr))
        self.g.attrs["shape"] = (self.rows, self.n_genes)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="folder with context_*.h5ad, gene_names.csv, pert_counts.csv")
    ap.add_argument("--h1", required=True, help="h1_pseudobulk.npz (log1p CP10k effects of 150 H1 knockdowns)")
    ap.add_argument("--k562", required=True, help="Replogle K562_gwps_normalized_bulk_01.h5ad")
    ap.add_argument("--method-dir", required=True, help="folder holding xline_benchmark.py")
    ap.add_argument("--out", required=True, help="combined prediction .h5ad")
    ap.add_argument("--diag", default=None, help="folder for diagnostics (log changes, per-target info)")
    ap.add_argument("--contexts", default="A,B,C")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--preset", default="c1_calibrated", choices=list(PRESETS))
    a = ap.parse_args()
    global PARAMS
    PARAMS = PRESETS[a.preset]
    sys.path.insert(0, str(Path(a.method_dir).resolve()))
    import xline_benchmark as xb  # noqa: E402

    data = Path(a.data)
    genes = pd.read_csv(data / "gene_names.csv")["gene_name"].astype(str).tolist()
    targets = pd.read_csv(data / "pert_counts.csv")["target_gene"].astype(str).tolist()
    contexts = a.contexts.split(",")

    # H1 effects moved onto the 2026 gene axis (3 H1 genes are not in it; 456 2026 genes get 0)
    z = np.load(a.h1, allow_pickle=True)
    col = pd.Series(np.arange(len(z["genes"])), index=z["genes"].astype(str))
    pos = col.reindex(genes).to_numpy()
    have = np.isfinite(pos)
    D = np.zeros((len(z["perts"]), len(genes)), dtype=np.float32)
    D[:, have] = z["delta"][:, pos[have].astype(int)]
    h1 = dict(delta=D, perts=z["perts"].astype(str))
    k562 = xb.load_k562(a.k562)

    diag = Path(a.diag) if a.diag else None
    if diag:
        diag.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    obs_t, obs_c = [], []
    with h5py.File(a.out, "w") as f:
        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"
        w = CSRWriter(f, len(genes))
        for ctx in contexts:
            X = ad.read_h5ad(data / f"context_{ctx}.h5ad").X.tocsr()
            fold, delta, info = predict_folds(ctx, X, genes, targets, h1, k562, xb)
            print(f"context {ctx}: {info.in_k562.sum()} targets in K562; median genes changed >10%: "
                  f"{int(info.n_genes_changed_10pct.median())}", flush=True)
            if diag:
                np.save(diag / f"log_change_{ctx}.npy", delta.astype(np.float32))
                np.save(diag / f"fold_{ctx}.npy", fold)
                info.to_csv(diag / f"targets_{ctx}.csv", index=False)
            for i, t in enumerate(targets):
                w.append(make_cells(X, fold[i], rng))
                obs_t += [t] * CELLS_PER_PERT
                obs_c += [ctx] * CELLS_PER_PERT
            print(f"  wrote {w.rows} cells, {w.nnz:,} stored counts so far", flush=True)
            del X
        w.close()
        obs = pd.DataFrame({"target_gene": pd.Categorical(obs_t, categories=targets),
                            "context": pd.Categorical(obs_c, categories=contexts)},
                           index=[f"pred_{i}" for i in range(len(obs_t))])
        ad.io.write_elem(f, "obs", obs)
        ad.io.write_elem(f, "var", pd.DataFrame(index=genes))
        for k in ("obsm", "varm", "obsp", "varp", "layers", "uns"):
            ad.io.write_elem(f, k, {})
    json.dump({"preset": a.preset, "params": PARAMS, "seed": a.seed, "cells_per_pert": CELLS_PER_PERT, "n_cells": w.rows,
               "nnz": int(w.nnz), "contexts": contexts},
              open(Path(a.out).with_suffix(".json"), "w"), indent=1)
    print(f"done: {w.rows} cells x {len(genes)} genes, {w.nnz:,} stored counts")


if __name__ == "__main__":
    main()
