"""Do GRNBoost2 / GENIE3 networks help predict real H1 knockdowns?

Three uses of a network, each scored leave-one-out against the 147 real H1 knockdowns
(eval_h1.score: Pearson over genes excluding the target, top-100 Pearson, top-50 sign):

  A. propagate   the knocked-down gene's out-edges are its predicted effect:
                 pred_g = -imp[k, g] * sign(rho[k, g]) (1 step), or pushed through the
                 column-normalized signed network for 2-3 steps. Scored on the network's genes.
  B. borrow      gene similarity from the network (correlation of two knockdown targets'
                 out-edge or in-edge profiles) picks which other knockdowns' real effects
                 to average, exactly as screen_transfer_h1.py does with screen similarity.
                 Scored on the same 10,993 genes as RESULTS.md.
  C. combine     network similarity averaged with the public-screen similarity
                 (K562+RPE1+HepG2+Jurkat), to see whether it adds anything.

    python eval_grn_h1.py --net name=net.npz ... --out out_dir
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
for p in (HERE.parent / "public_screens", Path("/mnt/project-files/public-screens"), Path("/mnt/project-files/logo-shap/v2")):
    sys.path.insert(0, str(p))
try:
    from h1_eval_utils import score  # repo layout
except ImportError:
    from eval_h1 import score  # project-folder layout
from screen_transfer_h1 import corr_rows, load_screen, predict_from_sim, screen_similarity  # noqa: E402

PB = "/mnt/project-files/logo-shap/v2/h1/h1_pseudobulk.npz"
GENES = "/mnt/project-files/logo-shap/v2/h1/grn_h1_pc0.top_edges.csv"
SD = "/mnt/project-files/public-screens/data/"
SCREEN_FILES = {"k562": "K562_gwps_normalized_bulk_01.h5ad", "rpe1": "rpe1_normalized_bulk_01.h5ad",
                "hepg2": "hepg2_bulk.h5ad", "jurkat": "jurkat_bulk.h5ad"}


def signed(imp, rho):
    return imp * np.sign(rho) if rho is not None and rho.ndim == 2 else imp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", nargs="+", required=True, metavar="NAME=PATH")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", nargs="+", type=int, default=[5, 10, 20, 40])
    ap.add_argument("--pseudobulk", default=PB, help="h1_pseudobulk.npz (real H1 effects)")
    ap.add_argument("--genes-from", default=GENES, help="CSV whose 'target' column lists the 10,993 scoring genes")
    ap.add_argument("--screen-dir", default=SD, help="folder with the four screen .h5ad files")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    pb = np.load(a.pseudobulk, allow_pickle=True)
    delta_all, pb_genes, perts = pb["delta"].astype(np.float64), list(pb["genes"]), list(pb["perts"])
    genes = list(dict.fromkeys(pd.read_csv(a.genes_from, usecols=["target"])["target"]))
    gi_pb = {g: i for i, g in enumerate(pb_genes)}
    gi = {g: i for i, g in enumerate(genes)}
    keep = [p for p in perts if p in gi]
    n = len(keep)
    T = delta_all[[perts.index(p) for p in keep]][:, [gi_pb[g] for g in genes]]
    K = np.array([gi[p] for p in keep])
    mean_loo = (T.sum(0)[None, :] - T) / (n - 1)

    nets = {}
    for spec in a.net:
        nm, pth = spec.split("=", 1)
        d = np.load(pth, allow_pickle=True)
        rho = d["rho"] if d["rho"].ndim == 2 else None
        nets[nm] = (d["imp"].astype(np.float64), list(d["genes"]), rho)

    rows = []

    def add(method, m, i, pred, true, k, part):
        rows.append({"part": part, "method": method, "top": m, "pert": keep[i], **score(pred, true, k)})

    # ---- A. propagate through the network (scored on the network's own genes) ----
    for nm, (imp, ng, rho) in nets.items():
        ngi = {g: j for j, g in enumerate(ng)}
        cols = np.array([gi_pb[g] for g in ng])
        Tn = delta_all[[perts.index(p) for p in keep]][:, cols]
        Kn = np.array([ngi.get(p, -1) for p in keep])
        mean_n = (Tn.sum(0)[None, :] - Tn) / (n - 1)
        S = -signed(imp, rho)                                  # knockdown = drop in the regulator
        A = signed(imp, rho) / (np.abs(imp).sum(0, keepdims=True) + 1e-12)   # each target's in-weights sum to 1
        coex = rho if rho is not None else None
        for i in range(n):
            k = Kn[i]
            if k < 0:
                continue
            add("mean_effect", n - 1, i, mean_n[i], Tn[i], k, f"A:{nm}")
            p1 = S[k].copy()
            add("grn_1step", 1, i, p1, Tn[i], k, f"A:{nm}")
            d = np.zeros(len(ng)); d[k] = -1.0; acc = np.zeros(len(ng)); cur = d
            for s in range(3):
                cur = A.T @ cur; cur[k] = 0; acc += cur
            add("grn_3step", 3, i, acc, Tn[i], k, f"A:{nm}")
            if coex is not None:
                add("coexpression", 0, i, -coex[k], Tn[i], k, f"A:{nm}")
            # mean effect plus the network's knockdown-specific push, scale fit on the other knockdowns later
        # mean + beta * grn (beta fit leave-one-out by least squares on residuals)
        ok = Kn >= 0
        P1 = np.stack([S[Kn[i]] if ok[i] else np.zeros(len(ng)) for i in range(n)])
        for i in np.flatnonzero(ok):
            P1[i, Kn[i]] = 0
        R = Tn - mean_n
        num, den = (P1 * R).sum(1), (P1 * P1).sum(1)
        for i in np.flatnonzero(ok):
            o = ok & (np.arange(n) != i)
            beta = num[o].sum() / den[o].sum()
            add("mean+grn_1step", 1, i, mean_n[i] + beta * P1[i], Tn[i], Kn[i], f"A:{nm}")
        print(f"{nm}: propagation done ({len(ng)} genes, {ok.sum()} knockdowns); LOO beta ~ {num.sum() / den.sum():.3g}")

    # ---- B/C. similarity for borrowing (scored on the 10,993 genes) ----
    sims = {}
    for nm, (imp, ng, rho) in nets.items():
        ngi = {g: j for j, g in enumerate(ng)}
        have = [p for p in keep if p in ngi]
        idx = np.array([keep.index(p) for p in have]); kk = np.array([ngi[p] for p in have])
        Ssg = signed(imp, rho)
        for kind, M in (("out", Ssg[kk, :]), ("in", Ssg[:, kk].T)):
            M = M.copy()
            M[:, kk] = 0                                         # drop the targets' own entries
            S = np.full((n, n), np.nan); S[np.ix_(idx, idx)] = corr_rows(M)
            sims[f"{nm}_{kind}"] = S
        if rho is not None:
            S = np.full((n, n), np.nan); S[np.ix_(idx, idx)] = rho[np.ix_(kk, kk)]
            sims[f"{nm}_coexpr"] = S
    screens = {k: load_screen(str(Path(a.screen_dir) / v)) for k, v in SCREEN_FILES.items()}
    ss = [screen_similarity(df, sg, keep) for df, sg in screens.values()]
    sims["screens"] = np.nanmean(np.stack(ss), 0)
    for nm in list(sims):
        if nm != "screens" and not nm.endswith("coexpr"):
            sims[f"screens+{nm}"] = np.nanmean(np.stack([sims["screens"], sims[nm]]), 0)
    sims["oracle"] = corr_rows(T)

    for i in range(n):
        add("uniform", n - 1, i, mean_loo[i], T[i], K[i], "B")
    for nm, S in sims.items():
        for m in a.top:
            for i in range(n):
                pred = predict_from_sim(S, T, i, m)
                if pred is None:
                    pred = mean_loo[i]
                add(nm, m, i, pred, T[i], K[i], "B")
                if nm != "oracle":
                    add(nm + "_shrink0.5", m, i, 0.5 * mean_loo[i] + 0.5 * pred, T[i], K[i], "B")

    df = pd.DataFrame(rows)
    df.to_csv(out / "grn_h1_scores.csv", index=False)
    summ = df.groupby(["part", "method", "top"])[["pearson_all", "pearson_top100", "sign_top50"]].mean()
    summ.to_csv(out / "grn_h1_summary.csv")
    with pd.option_context("display.max_rows", 500, "display.width", 200):
        print(summ.round(3).to_string())


if __name__ == "__main__":
    main()
