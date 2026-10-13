"""Can we predict a held-out knockdown's effect from the *other* knockdowns?

For each knocked-down gene k, predict its effect as a similarity-weighted average of the
effects of the other 149 knockdowns:

    pred_k = sum_j w_kj * delta_j / sum_j w_kj,   w_kj = similarity(k, j)^+ over the m most similar j

Similarity between genes k and j comes from one of:
  uniform      every j weighted equally  (= the mean-effect baseline)
  coexpr       correlation of k's and j's expression across control cells
  grn_row      correlation of their regulator profiles in the control-cell network (W[:, k] vs W[:, j])
  grn_col      correlation of their target profiles (W[k, :] vs W[j, :])
  effect       correlation of their true effects (an oracle ceiling, uses the answer)

This is the "use control cells to decide which known effects to borrow" idea. Scored like
eval_h1.py: Pearson between predicted and real change over the network's genes, excluding k.

    python loo_transfer_h1.py --pseudobulk h1_pseudobulk.npz --grn grn_h1_pc0.npz \
        --controls h1_controls.h5ad --out loo_transfer.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import grn_ridge  # noqa: E402
from eval_h1 import score  # noqa: E402


def corr_rows(M):
    Z = M - M.mean(1, keepdims=True)
    Z /= np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12
    return Z @ Z.T


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--grn", required=True)
    ap.add_argument("--controls", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", nargs="+", type=int, default=[5, 15, 40])
    a = ap.parse_args()

    pb = np.load(a.pseudobulk, allow_pickle=True)
    delta_all, pb_genes, perts = pb["delta"].astype(np.float64), list(pb["genes"]), list(pb["perts"])
    g = np.load(a.grn, allow_pickle=True)
    genes = list(g["genes"]); W = g["coef_mean"]
    gi = {x: i for i, x in enumerate(genes)}
    cols = np.array([pb_genes.index(x) for x in genes])
    in_net = [p for p in perts if p in gi]
    P = np.array([perts.index(p) for p in in_net])
    K = np.array([gi[p] for p in in_net])
    T = delta_all[P][:, cols]                                  # perts x network genes (true effects)
    print(f"{len(in_net)} knockdowns in the {len(genes)}-gene network")

    # Only the knocked-down genes' expression is needed: read sparse, normalize, slice columns.
    import anndata as ad
    import scipy.sparse as sp
    ctrl = ad.read_h5ad(a.controls)
    Xs = sp.csr_matrix(ctrl.X, dtype=np.float32)
    lib = np.asarray(Xs.sum(1)).ravel()
    ci = {x: i for i, x in enumerate(ctrl.var_names)}
    Xk = Xs[:, [ci[p] for p in in_net]].toarray().astype(np.float64)
    Xk = np.log1p(Xk * (1e4 / np.maximum(lib, 1))[:, None])
    del ctrl, Xs
    sims = {
        "uniform": np.ones((len(in_net), len(in_net))),
        "coexpr": corr_rows(Xk.T),
        "grn_row": corr_rows(W[:, K].T),
        "grn_col": corr_rows(W[K, :]),
        "effect": corr_rows(T),
    }
    del Xk, W

    rows = []
    for name, S in sims.items():
        S = S.copy(); np.fill_diagonal(S, -np.inf)
        for m in ([len(in_net) - 1] if name == "uniform" else a.top):
            for i, p in enumerate(in_net):
                order = np.argsort(-S[i])[:m]
                w = np.clip(S[i, order], 0, None) if name != "uniform" else np.ones(m)
                if w.sum() == 0:
                    w = np.ones(m)
                pred = (w[:, None] * T[order]).sum(0) / w.sum()
                r = score(pred, T[i], K[i])
                rows.append({"method": name, "top": m, "pert": p, **r})
    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    print(df.groupby(["method", "top"])[["pearson_all", "pearson_top100", "sign_top50"]].mean().round(3).to_string())


if __name__ == "__main__":
    main()
