"""Borrow known H1 knockdown effects, weighted by gene similarity learned from a public screen.

For each H1 knockdown k (leave-one-out over the 150), predict its change as a weighted average
of the other knockdowns' real H1 changes, where the weights come from how alike k and j look in
a public Perturb-seq screen (Replogle 2022 K562 genome-wide and RPE1, Nadig 2025 HepG2 and Jurkat; z-scored pseudobulk):

    pred_k = sum_j w_kj * delta_j / sum_j w_kj,   w_kj = sim(k, j)^+ over the m most similar j

Also tested:
  shrink     pred = (1-a) * mean_effect + a * similarity-weighted average, a in {0.5, 0.75}
  direct     mean_effect + beta * (k's own screen profile mapped onto H1 genes), beta fit on the
             other knockdowns (leave-one-out), so no H1 answer for k leaks in
  uniform    every other knockdown weighted equally (= the 0.19 mean-effect baseline)
  oracle     similarity from k's true H1 effect (the 0.46 ceiling)

Knockdowns absent from a screen fall back to uniform. Scored with eval_h1.score over the same
10,993 genes as the earlier RESULTS.md runs (target gene excluded).

    python screen_transfer_h1.py --pseudobulk h1_pseudobulk.npz --genes-from grn_h1_pc0.top_edges.csv \
        --screen k562=K562_gwps_normalized_bulk_01.h5ad rpe1=rpe1_normalized_bulk_01.h5ad \
        hepg2=hepg2_bulk.h5ad jurkat=jurkat_bulk.h5ad --out out_dir/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from h1_eval_utils import score  # noqa: E402


def corr_rows(M):
    Z = M - M.mean(1, keepdims=True)
    Z /= np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12
    return Z @ Z.T


def load_screen(path):
    """Return (profiles: perturbed gene -> z-profile, readout gene names). One row per target gene."""
    import anndata as ad
    a = ad.read_h5ad(path)
    X = np.asarray(a.X, dtype=np.float64)
    X[~np.isfinite(X)] = 0                              # z-scores are inf where a readout gene has zero control variance
    np.clip(X, -10, 10, out=X)
    tgt = a.obs.index.str.split("_").str[1].to_numpy()
    df = pd.DataFrame(X).groupby(tgt).mean()          # a few genes have 2 transcripts: average them
    return df, a.var["gene_name"].to_numpy()


def screen_similarity(df, sgenes, perts, exclude_targets=True):
    """Correlation of screen profiles between H1 knockdowns. NaN where either is missing."""
    have = [p for p in perts if p in df.index]
    M = df.loc[have].to_numpy().copy()
    if exclude_targets:                                 # don't let the knocked-down genes' own drops drive similarity
        M = M[:, ~np.isin(sgenes, have)]
    S_have = corr_rows(M)
    S = np.full((len(perts), len(perts)), np.nan)
    idx = np.array([perts.index(p) for p in have])
    S[np.ix_(idx, idx)] = S_have
    return S


def predict_from_sim(S, T, i, m, signed=False):
    s = S[i].copy(); s[i] = np.nan
    ok = np.flatnonzero(np.isfinite(s))
    if len(ok) == 0 or m == 0:
        return None
    key = np.abs(s[ok]) if signed else s[ok]
    order = ok[np.argsort(-key)[:m]]
    w = s[order] if signed else np.clip(s[order], 0, None)
    if np.abs(w).sum() == 0:
        return None
    return (w[:, None] * T[order]).sum(0) / np.abs(w).sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--genes-from", required=True, help="CSV whose 'target' column lists the genes to score on")
    ap.add_argument("--screen", nargs="+", required=True, metavar="NAME=PATH",
                    help="one or more screens, e.g. k562=K562_gwps_normalized_bulk_01.h5ad rpe1=...")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", nargs="+", type=int, default=[5, 10, 20, 40])
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    pb = np.load(a.pseudobulk, allow_pickle=True)
    delta_all, pb_genes, perts = pb["delta"].astype(np.float64), list(pb["genes"]), list(pb["perts"])
    genes = list(dict.fromkeys(pd.read_csv(a.genes_from, usecols=["target"])["target"]))
    gi_pb = {g: i for i, g in enumerate(pb_genes)}
    cols = np.array([gi_pb[g] for g in genes])
    gi = {g: i for i, g in enumerate(genes)}
    keep = [p for p in perts if p in gi]
    P = np.array([perts.index(p) for p in keep])
    T = delta_all[P][:, cols]                     # knockdowns x genes, real H1 changes
    K = np.array([gi[p] for p in keep])
    n = len(keep)
    print(f"{n} H1 knockdowns scored over {len(genes)} genes")

    screens = {nm: load_screen(pth) for nm, pth in (x.split("=", 1) for x in a.screen)}
    sims = {}
    for name, (df, sgenes) in screens.items():
        sims[name] = screen_similarity(df, sgenes, keep)
        print(f"{name}: {np.isfinite(np.diag(sims[name])).sum()} of {n} knockdowns in screen")
    if len(sims) > 1:    # average over the screens that have both genes, else whichever has them
        sims["+".join(sims)] = np.nanmean(np.stack(list(sims.values())), 0)
        if {"k562", "rpe1"} <= set(sims):
            sims["k562+rpe1"] = np.nanmean(np.stack([sims["k562"], sims["rpe1"]]), 0)
    sims["oracle"] = corr_rows(T)

    mean_loo = (T.sum(0)[None, :] - T) / (n - 1)
    rows = []

    def add(method, m, i, pred, covered):
        rows.append({"method": method, "top": m, "pert": keep[i], "in_screen": covered, **score(pred, T[i], K[i])})

    for i in range(n):
        add("uniform", n - 1, i, mean_loo[i], True)

    for name, S in sims.items():
        for m in a.top:
            for signed in (False, True):
                if signed and name == "oracle":
                    continue
                tag = name + ("_signed" if signed else "")
                for i in range(n):
                    pred = predict_from_sim(S, T, i, m, signed)
                    covered = pred is not None
                    if pred is None:
                        pred = mean_loo[i]
                    add(tag, m, i, pred, covered)
                    if not signed and name != "oracle":
                        for alpha in (0.5, 0.75):
                            add(f"{name}_shrink{alpha}", m, i, (1 - alpha) * mean_loo[i] + alpha * pred, covered)

    # Direct: k's own screen profile, mapped onto H1 genes, scaled by beta fit on the other knockdowns.
    for name, (df, sgenes) in screens.items():
        sg = {g: c for c, g in enumerate(sgenes)}
        shared = np.array([gi[g] for g in sgenes if g in gi])
        scol = np.array([sg[g] for g in sgenes if g in gi])
        prof = np.zeros((n, len(genes))); has = np.zeros(n, bool)
        for i, p in enumerate(keep):
            if p in df.index:
                prof[i, shared] = df.loc[p].to_numpy()[scol]; prof[i, K[i]] = 0; has[i] = True
        R = T - mean_loo                                  # residual after the shared response
        num = (prof * R).sum(1); den = (prof * prof).sum(1)
        for i in range(n):
            o = has & (np.arange(n) != i)
            beta = num[o].sum() / den[o].sum()
            add(f"direct_{name}", 0, i, mean_loo[i] + beta * prof[i], bool(has[i]))
        print(f"direct_{name}: beta ~ {num[has].sum() / den[has].sum():.4f}")

    df = pd.DataFrame(rows)
    df.to_csv(out / "screen_transfer.csv", index=False)
    cols_ = ["pearson_all", "pearson_top100", "sign_top50"]
    summ = df.groupby(["method", "top"])[cols_].mean()
    cov = df[df.in_screen].groupby(["method", "top"])["pearson_all"].agg(["mean", "size"])
    summ["pearson_all_in_screen"] = cov["mean"]; summ["n_in_screen"] = cov["size"]
    summ = summ.sort_values("pearson_all", ascending=False)
    summ.to_csv(out / "screen_transfer_summary.csv")
    print(summ.round(3).to_string())


if __name__ == "__main__":
    main()
