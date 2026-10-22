"""Sanity checks on a built 2026 prediction (build_submission.py output).

1. The 13 targets that were also knocked down in H1: does the predicted change (which never
   uses the target's own H1 effect) agree with that real H1 effect? Reported for the whole
   profile and for the specific part (each profile minus the mean over targets), next to what
   the shared response alone achieves.
2. Effect sizes per context: target-gene drop, genes moved by more than 10%.
3. Read-back: per-target pseudobulk of the written cells against the intended fold changes,
   and the obs / cell-count / integer checks that vcc prep also runs.

    python sanity_checks.py --pred pred.h5ad --diag diag/ --h1 h1_pseudobulk.npz \
        --data vcc2026-data/ --out sanity/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import bundle


def corr(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--diag", required=True)
    ap.add_argument("--h1", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contexts", default=None, help="comma-separated; default: the bundle's manifest.json")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    data, diag = Path(a.data), Path(a.diag)
    genes = pd.read_csv(data / "gene_names.csv")["gene_name"].astype(str).tolist()
    gi = {g: i for i, g in enumerate(genes)}
    contexts = bundle.contexts(data, a.contexts)
    targets = bundle.targets(data, contexts)[contexts[0]]   # the 2026 lists are shared by all contexts

    z = np.load(a.h1, allow_pickle=True)
    h1g = z["genes"].astype(str); h1p = list(z["perts"].astype(str))
    common = [g for g in h1g if g in gi]
    hc = np.array([list(h1g).index(g) for g in common]); pc = np.array([gi[g] for g in common])
    tgt_set = set(targets)
    keep = np.array([g not in tgt_set for g in common])          # drop all panel targets, as pds does
    hc, pc = hc[keep], pc[keep]
    H = z["delta"][:, hc]; Hmean = H.mean(0)
    shared = [t for t in targets if t in h1p]

    rows = []
    for ctx in contexts:
        L = np.load(diag / f"log_change_{ctx}.npy")[:, pc]
        Lmean = L.mean(0)
        for t in shared:
            i = targets.index(t); h = H[h1p.index(t)]
            rows.append({"context": ctx, "target": t,
                         "r_pred_vs_h1": corr(L[i], h),
                         "r_shared_only_vs_h1": corr(Hmean, h),
                         "r_specific_vs_h1_specific": corr(L[i] - Lmean, h - Hmean),
                         "pred_specific_norm": float(np.linalg.norm(L[i] - Lmean)),
                         "h1_specific_norm": float(np.linalg.norm(h - Hmean))})
    h1df = pd.DataFrame(rows)
    h1df.to_csv(out / "h1_covered_targets.csv", index=False)

    eff = []
    for ctx in contexts:
        info = pd.read_csv(diag / f"targets_{ctx}.csv")
        fold = np.load(diag / f"fold_{ctx}.npy")
        self_fold = np.array([fold[i, gi[t]] for i, t in enumerate(targets)])
        eff.append({"context": ctx, "targets_in_k562": int(info.in_k562.sum()),
                    "median_target_remaining": float(np.median(self_fold)),
                    "median_genes_moved_10pct": float(info.n_genes_changed_10pct.median()),
                    "max_genes_moved_10pct": int(info.n_genes_changed_10pct.max())})
    effdf = pd.DataFrame(eff)
    effdf.to_csv(out / "effect_sizes.csv", index=False)

    # read-back of the written file
    p = ad.read_h5ad(a.pred, backed="r")
    obs = p.obs
    counts = obs.groupby(["context", "target_gene"], observed=True).size()
    rb = {"n_cells": int(p.n_obs), "n_genes": int(p.n_vars), "genes_in_order": p.var_names.tolist() == genes,
          "contexts": sorted(obs["context"].astype(str).unique().tolist()),
          "all_targets_every_context": bool(all(set(obs.loc[obs.context == c, "target_gene"].astype(str)) == tgt_set
                                                for c in contexts)),
          "cells_per_target_min": int(counts.min()), "cells_per_target_max": int(counts.max()),
          "control_cells": int((obs["target_gene"].astype(str) == "non-targeting").sum())}
    rng = np.random.default_rng(1)
    fid = []
    for ctx in contexts:
        X = ad.read_h5ad(bundle.control_file(data, ctx)).X
        cmean = np.asarray(X.mean(0)).ravel()
        fold = np.load(diag / f"fold_{ctx}.npy")
        for t in rng.choice(targets, 5, replace=False):
            rows_ = np.flatnonzero(((obs.context == ctx) & (obs.target_gene == t)).to_numpy())
            S = p.X[rows_[0]:rows_[-1] + 1]
            d = S.data
            ok = bool(np.all(d >= 0) and np.all(d == np.round(d)))
            pm = np.asarray(S.mean(0)).ravel()
            expressed = cmean > 0.5
            fid.append({"context": ctx, "target": t, "integer_nonneg": ok,
                        "max_cell_total": float(np.asarray(S.sum(1)).max()),
                        "r_log_written_vs_intended": corr(np.log((pm[expressed] + 1e-3) / cmean[expressed]),
                                                          np.log(fold[targets.index(t), expressed]))})
    fiddf = pd.DataFrame(fid)
    fiddf.to_csv(out / "readback.csv", index=False)
    pd.Series(rb).to_csv(out / "file_checks.csv", header=["value"])

    pd.set_option("display.width", 200)
    print(pd.Series(rb).to_string())
    print(effdf.round(3).to_string(index=False))
    print(fiddf.round(3).to_string(index=False))
    print(h1df.round(3).to_string(index=False))
    print(h1df.groupby("context")[["r_pred_vs_h1", "r_shared_only_vs_h1", "r_specific_vs_h1_specific"]]
          .median().round(3).to_string())


if __name__ == "__main__":
    main()
