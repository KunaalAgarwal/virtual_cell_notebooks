"""cell-eval on real H1 knockdowns for screen-similarity effect borrowing, next to the baselines.

Builds each submission the same way as logo-shap/v2/celleval_h1.py (resample 400 control
cells per knockdown and shift them by the predicted change), with these methods:

  ctrl_sampled   no change (chance)
  mean_effect    leave-one-out mean of the other 149 H1 knockdowns (the 0.19 Pearson baseline)
  screen_sim     similarity-weighted mean of the other knockdowns, similarity = correlation of
                 their profiles in public screens (default: Replogle K562 GWPS + RPE1), top-m,
                 shrunk halfway to mean_effect; falls back to mean_effect if k is in no screen
  oracle         the real change
  ctrl_mean      every cell = rounded control mean (the repo baseline)
  screen_sim_aC  mean_effect + C x (screen_sim - mean_effect): the knockdown-specific part amplified
                 C times (--amp). cell-eval rewards a confident, distinct change per knockdown.

The knocked-down gene itself is set to the H1 median drop in every non-oracle method.

    python celleval_screen_h1.py --real h1_real_test.h5ad --pseudobulk h1_pseudobulk.npz \
        --screens K562_gwps_normalized_bulk_01.h5ad rpe1_normalized_bulk_01.h5ad --out out/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from h1_eval_utils import SELF_DROP, build_pred, control_mean_pred  # noqa: E402
from screen_transfer_h1 import load_screen, predict_from_sim, screen_similarity  # noqa: E402


def screen_sim_deltas(delta_all, perts, screens, m, alpha):
    sims = [screen_similarity(df, sg, perts) for df, sg in screens]
    S = np.nanmean(np.stack(sims), 0) if len(sims) > 1 else sims[0]
    n = len(perts)
    mean_loo = (delta_all.sum(0)[None, :] - delta_all) / (n - 1)
    out, covered = {}, 0
    for i in range(n):
        p = predict_from_sim(S, delta_all, i, m)
        if p is None:
            out[perts[i]] = mean_loo[i]
        else:
            out[perts[i]] = (1 - alpha) * mean_loo[i] + alpha * p
            covered += 1
    print(f"screen_sim: {covered}/{n} knockdowns covered by a screen")
    return out, mean_loo


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--real", required=True)
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--screens", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--alpha", type=float, default=0.5)
    ap.add_argument("--cells-per-pert", type=int, default=400)
    ap.add_argument("--methods", nargs="+", default=["ctrl_sampled", "ctrl_mean", "mean_effect", "screen_sim", "oracle"])
    ap.add_argument("--amp", nargs="*", type=float, default=[], help="amplification factors for screen_sim_aC")
    ap.add_argument("--scale-amp", nargs="*", default=[], metavar="S:C",
                    help="variants S x mean_effect + C x (screen_sim - mean_effect), named screen_sS_aC")
    ap.add_argument("--shuffled-amp", nargs="*", type=float, default=[],
                    help="null control: give each knockdown another knockdown's screen-specific part, amplified C times")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    real = ad.read_h5ad(a.real)
    if np.modf(real.X.data[:100000])[0].max() > 1e-3:   # the 2025 support-set file stores log1p(counts)
        real.X.data = np.rint(np.expm1(real.X.data)).astype(np.float32)
        print("real file was log1p(counts): converted back to counts")
    real.obs["target"] = real.obs["target"].astype(str)
    test = sorted(p for p in real.obs["target"].unique() if p != "non-targeting")
    genes = real.var_names.to_numpy()
    gi = {g: i for i, g in enumerate(genes)}
    print(f"real: {real.n_obs} cells, {len(test)} knockdowns, {len(genes)} genes")

    pb = np.load(a.pseudobulk, allow_pickle=True)
    assert list(pb["genes"]) == list(genes), "gene order differs between pseudobulk and real file"
    perts = [str(p) for p in pb["perts"]]
    delta_all = pb["delta"].astype(np.float64)

    screens = [load_screen(s) for s in a.screens]
    sim, mean_loo = screen_sim_deltas(delta_all, perts, screens, a.top, a.alpha)

    def with_self(d, p):
        d = d.copy(); d[gi[p]] = SELF_DROP; return d

    deltas = {
        "ctrl_sampled": {p: np.zeros(len(genes)) for p in test},
        "mean_effect": {p: with_self(mean_loo[perts.index(p)], p) for p in test},
        "screen_sim": {p: with_self(sim[p], p) for p in test},
        "oracle": {p: delta_all[perts.index(p)] for p in test},
    }

    for c in a.amp:
        name = f"screen_sim_a{c:g}"
        deltas[name] = {p: with_self(mean_loo[perts.index(p)] + c * (sim[p] - mean_loo[perts.index(p)]), p) for p in test}
        if name not in a.methods:
            a.methods.append(name)

    for v in a.scale_amp:
        sc, c = (float(x) for x in v.split(":"))
        name = f"screen_s{sc:g}_a{c:g}"
        deltas[name] = {p: with_self(sc * mean_loo[perts.index(p)] + c * (sim[p] - mean_loo[perts.index(p)]), p)
                        for p in test}
        if name not in a.methods:
            a.methods.append(name)

    perm = np.random.default_rng(a.seed + 1).permutation(len(perts))
    for c in a.shuffled_amp:
        name = f"shuffled_a{c:g}"
        deltas[name] = {}
        for p in test:
            i = perts.index(p); j = perm[i] if perm[i] != i else (i + 1) % len(perts)
            deltas[name][p] = with_self(mean_loo[i] + c * (sim[perts[j]] - mean_loo[j]), p)
        if name not in a.methods:
            a.methods.append(name)

    from cell_eval import MetricsEvaluator
    rows = []
    for m in a.methods:
        if m == "ctrl_mean":
            pred = control_mean_pred(real, test, a.cells_per_pert)
        else:
            pred = build_pred(real, deltas[m], a.cells_per_pert, np.random.default_rng(a.seed))
        print(f"== {m}: pred {pred.n_obs} cells", flush=True)
        ev = MetricsEvaluator(adata_pred=pred, adata_real=real.copy(), control_pert="non-targeting",
                              pert_col="target", outdir=str(out / m), num_threads=a.threads)
        res, _ = ev.compute(profile="vcc", write_csv=True)
        res = res.to_pandas(); res["method"] = m
        rows.append(res)
        print(res.drop(columns=["method"]).mean(numeric_only=True).round(4).to_string(), flush=True)

    df = pd.concat(rows, ignore_index=True)
    df.to_csv(out / "per_pert_metrics.csv", index=False)
    summ = df.groupby("method").mean(numeric_only=True)
    summ.to_csv(out / "summary.csv")
    print("\n" + summ.round(4).to_string())


if __name__ == "__main__":
    main()
