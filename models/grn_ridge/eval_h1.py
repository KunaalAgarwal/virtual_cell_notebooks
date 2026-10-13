"""Simulate H1 knockdowns through a control-learned GRN and score against the real ones.

For each knocked-down gene k that is in the network, the knockdown's change propagates
through the edges (stage 2 of grn_linear_model.pdf, in its well-posed "change" form):

    d = e_k * d_k + W' d      (iterated; gene k stays clamped at d_k)

W[j, i] is the ridge slope of regulator j on target i (log1p CP10k units). Because the
largest |eigenvalue| of W is < 1 the iteration converges; `steps` caps it (1 = direct
targets only).

d_k (how far the knocked-down gene itself drops) is either the observed value (oracle)
or the median over all H1 knockdowns (what we'd have to assume zero-shot).

Scored per knockdown over the network's genes, excluding k itself:
  pearson_all     Pearson(predicted change, real change) over all genes
  pearson_top100  same, on the 100 genes the knockdown changed most
  sign_top50      fraction of the 50 most-changed genes whose direction is right

Baselines:
  mean_effect     the average real change over the other 149 knockdowns (no GRN)
  coexpression    marginal slope cov(i, k) / var(k) from controls (no conditioning)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import grn_ridge  # noqa: E402


def propagate(W, idx, d_self, steps):
    """Change in every gene for each knockdown (columns), via clamped iteration."""
    n, p = W.shape[0], len(idx)
    D0 = np.zeros((n, p), dtype=np.float64)
    D0[idx, np.arange(p)] = d_self
    D = D0.copy()
    Wt = W.T.astype(np.float64)
    for _ in range(steps):
        D = D0 + Wt @ D
        D[idx, np.arange(p)] = d_self
    return D


def score(pred, true, k):
    mask = np.ones(len(true), bool)
    mask[k] = False
    p, t = pred[mask], true[mask]
    out = {}
    out["pearson_all"] = np.corrcoef(p, t)[0, 1] if p.std() > 0 else np.nan
    top = np.argsort(-np.abs(t))
    t100 = top[:100]
    out["pearson_top100"] = np.corrcoef(p[t100], t[t100])[0, 1] if p[t100].std() > 0 else np.nan
    t50 = top[:50]
    out["sign_top50"] = float(np.mean(np.sign(p[t50]) == np.sign(t[t50])))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grn", nargs="+", required=True, help="grn_ridge .npz files to evaluate")
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--controls", required=True, help="H1 controls .h5ad (for the co-expression baseline)")
    ap.add_argument("--out", required=True, help="output CSV of per-knockdown scores")
    a = ap.parse_args()

    pb = np.load(a.pseudobulk, allow_pickle=True)
    delta_all, pb_genes, perts = pb["delta"].astype(np.float64), list(pb["genes"]), pb["perts"]
    gi_all = {g: i for i, g in enumerate(pb_genes)}
    self_all = np.array([delta_all[i, gi_all[p]] for i, p in enumerate(perts)])
    d_fixed = float(np.median(self_all))
    print(f"median self-knockdown change {d_fixed:.3f} (log1p CP10k)")

    rows = []
    coexp_done = set()
    for path in a.grn:
        g = np.load(path, allow_pickle=True)
        genes = list(g["genes"])
        W = g["coef_mean"]
        name = Path(path).stem
        cols = np.array([gi_all[x] for x in genes])
        true = delta_all[:, cols].T                       # network genes x perts
        in_net = np.array([p in set(genes) for p in perts])
        idx = np.array([genes.index(p) for p in perts[in_net]])
        P = np.flatnonzero(in_net)
        T = true[:, P]
        print(f"{name}: {len(genes)} genes, {len(P)}/{len(perts)} knockdowns in network")

        preds = {}
        for steps in (1, 3, 20):
            preds[f"grn_oracle_{steps}step"] = propagate(W, idx, T[idx, np.arange(len(P))], steps)
        preds["grn_fixed_20step"] = propagate(W, idx, np.full(len(P), d_fixed), 20)

        # Leave-one-out mean effect over the other knockdowns, alone and plus the GRN shift.
        tot = true.sum(1, keepdims=True)
        mean_eff = (tot - true)[:, P] / (true.shape[1] - 1)
        preds["mean_effect+grn"] = mean_eff + preds["grn_oracle_20step"]
        if "mean_effect" not in coexp_done:
            preds["mean_effect"] = mean_eff
            # Co-expression: marginal slope of every gene on k in the raw (no state removal) controls.
            X, cgenes, _, _ = grn_ridge.load(a.controls, "batch", 0.0)
            ci = {x: i for i, x in enumerate(cgenes)}
            Xg = X[:, [ci[x] for x in genes]].astype(np.float64)
            Xg -= Xg.mean(0)
            Xk = Xg[:, idx]
            slope = (Xg.T @ Xk) / (Xk**2).sum(0)              # genes x perts
            preds["coexpression"] = slope * T[idx, np.arange(len(P))]
            del X, Xg
            coexp_done.add("mean_effect")

        for method, M in preds.items():
            label = method if method in ("mean_effect", "coexpression") else f"{name}:{method}"
            for c, pi in enumerate(P):
                r = score(M[:, c], T[:, c], idx[c])
                rows.append({"method": label, "pert": perts[pi], "effect_size": np.abs(T[:, c]).mean(), **r})

    df = pd.DataFrame(rows)
    df.to_csv(a.out, index=False)
    strong = df.effect_size >= df.groupby("method").effect_size.transform(lambda s: s.quantile(0.75))
    summ = df.groupby("method")[["pearson_all", "pearson_top100", "sign_top50"]].mean()
    summ["pearson_all_strongest25%"] = df[strong].groupby("method").pearson_all.mean()
    print(summ.round(3).sort_values("pearson_all", ascending=False).to_string())
    summ.round(4).to_csv(Path(a.out).with_suffix(".summary.csv"))


if __name__ == "__main__":
    main()
