"""Robustness check: learn the control-cell network on metacells instead of single cells.

Single-cell counts are sparse and noisy. A metacell averages each cell with its k nearest
neighbours (in PCA space), which removes most sampling noise. If the network only failed
because single cells are too noisy, the metacell version should predict knockdowns better.

Writes a grn_ridge-style .npz so eval_h1.py can score it.

    python metacell_grn.py --h5ad h1_controls.h5ad --group-col batch --k 30 --out grn_h1_meta30.npz
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
import grn_ridge  # noqa: E402


def metacell_smooth(X, k=30, n_pcs=30, seed=0, log=print):
    """Replace each cell by the mean of itself and its k nearest neighbours.

    Neighbours are found in PCA space of the standardized matrix; the averaging is done on
    the (normalized, log) expression itself. Neighbour sets overlap, so the result is a
    smoothed version of the same cells, not a smaller set.
    """
    from sklearn.neighbors import NearestNeighbors
    from sklearn.utils.extmath import randomized_svd

    t0 = time.perf_counter()
    Xs = (X - X.mean(0)) / (X.std(0) + 1e-8)
    u, s, _ = randomized_svd(Xs, n_components=n_pcs, random_state=seed)
    pcs = u * s
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=-1).fit(pcs)
    _, idx = nn.kneighbors(pcs)                      # includes the cell itself
    n = X.shape[0]
    rows = np.repeat(np.arange(n), k + 1)
    A = sp.csr_matrix((np.full(rows.size, 1.0 / (k + 1), dtype=np.float32), (rows, idx.ravel())), shape=(n, n))
    Xm = (A @ X).astype(np.float32)
    log(f"  metacells: k={k}, {n_pcs} PCs ({time.perf_counter() - t0:.0f}s)")
    return Xm


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--group-col", default="batch")
    p.add_argument("--k", type=int, default=30)
    p.add_argument("--n-pcs", type=int, default=0, help="PCs to regress out AFTER smoothing (0 = none)")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--min-frac", type=float, default=0.2)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    X, genes, groups, covars = grn_ridge.load(a.h5ad, a.group_col, a.min_frac)
    print(f"{X.shape[0]} cells x {X.shape[1]} genes, {len(np.unique(groups))} groups")
    X = metacell_smooth(X, a.k)
    if a.n_pcs:
        X = grn_ridge.remove_cell_state(X, covars, a.n_pcs)
    res = grn_ridge.fit(X, groups, a.folds)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, genes=genes.astype(str), **res)
    grn_ridge.top_edges(res, genes).to_csv(out.with_suffix(".top_edges.csv"), index=False)
    print(f"lambda {res['lam']:g}, mean held-out R2 {res['r2'].mean():.3f}; wrote {out}")


if __name__ == "__main__":
    main()
