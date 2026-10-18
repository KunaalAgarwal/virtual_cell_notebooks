"""Pseudobulk sampling noise of 400 cells in the scorer's space (sum counts, normalize 5e4, log1p),
random vs mean-matched subsets, context A."""
import numpy as np, anndata as ad, scipy.sparse as sp
X = ad.read_h5ad("/mnt/project-files/vcc2026-data/context_A.h5ad").X.tocsr()
lib = np.asarray(X.sum(1)).ravel()
def pb(rows):
    s = np.asarray(X[rows].sum(0)).ravel(); return np.log1p(s / s.sum() * 5e4)
full = pb(np.arange(X.shape[0]))
cpm = np.asarray(X.sum(0)).ravel() / lib.sum() * 1e6
for name, mask in (("all genes", np.ones_like(cpm, bool)), ("genes >= 5 CPM", cpm >= 5)):
    rng = np.random.default_rng(0)
    rnd = [np.linalg.norm((pb(rng.choice(X.shape[0], 400, replace=False)) - full)[mask]) for _ in range(20)]
    # matched: draw 20 random subsets, keep the one closest to the full mean (cheap stratification)
    best = []
    for _ in range(5):
        cands = [rng.choice(X.shape[0], 400, replace=False) for _ in range(20)]
        best.append(min(np.linalg.norm((pb(c) - full)[mask]) for c in cands))
    # per-cell log-normalized noise (what each cell carries) for scale: median cell norm vs full
    print(f"{name}: {mask.sum()} genes; random 400-cell pseudobulk deviation norm median {np.median(rnd):.2f}; "
          f"best of 20 random subsets {np.median(best):.2f}")
