"""cell-eval on a held-out cell line (HepG2 or Jurkat) for settings chosen by xline_benchmark.py.

The support-set files are log-normalized, so predictions stay in log space: each knockdown
gets as many sampled control cells as the real data has, shifted by the predicted change and
clipped at zero. Methods: resampled controls, control mean in every cell (the repo baseline,
in log space), the round-3 x8 reference, the settings tuned on the *other* line for each
objective (chosen_by_objective.json), the same gene's effect in the other line, and the oracle.

    python celleval_xline.py --line jurkat.h5 --other hepg2.h5 --h1 h1_pseudobulk.npz \
        --k562 K562_gwps_normalized_bulk_01.h5ad --chosen chosen_by_objective.json --out out/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
import xline_benchmark as xb  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", required=True)
    ap.add_argument("--other", required=True)
    ap.add_argument("--h1", required=True)
    ap.add_argument("--k562", required=True)
    ap.add_argument("--chosen", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    real = ad.read_h5ad(a.line)
    real.obs = pd.DataFrame({"target": real.obs["target_gene"].astype(str).to_numpy()}, index=real.obs_names)
    X = np.asarray(real.X, dtype=np.float32)
    real.X = sp.csr_matrix(X)
    t = real.obs["target"].to_numpy()
    ctrl = X[t == xb.CTRL]

    L = xb.load_line(a.line); O = xb.load_line(a.other)
    F = xb.build_features(L, np.load(a.h1, allow_pickle=True), xb.load_k562(a.k562))
    chosen = json.load(open(a.chosen))[L["name"]]
    deltas = {"ctrl_sampled": np.zeros((len(L["perts"]), len(L["genes"])))}
    for name, p in chosen.items():
        deltas[name] = xb.predict(L, F, p)
    deltas["same gene, other line"] = np.stack([O["real"][O["perts"].index(q)] for q in L["perts"]])
    deltas["oracle"] = L["real"]

    def build(delta):
        blocks, obs = [ctrl], [xb.CTRL] * len(ctrl)
        for i, p in enumerate(L["perts"]):
            n = int((t == p).sum())
            cells = ctrl[rng.choice(len(ctrl), n, replace=n > len(ctrl))] + delta[i].astype(np.float32)
            blocks.append(np.clip(cells, 0, None)); obs += [p] * n
        o = pd.DataFrame({"target": obs}); o.index = [f"c{i}" for i in range(len(o))]
        return ad.AnnData(X=sp.csr_matrix(np.vstack(blocks)), obs=o, var=real.var.copy())

    def ctrl_mean_pred():
        mean_row = np.log1p(np.expm1(ctrl).mean(0))
        blocks, obs = [ctrl], [xb.CTRL] * len(ctrl)
        for p in L["perts"]:
            n = int((t == p).sum()); blocks.append(np.tile(mean_row, (n, 1))); obs += [p] * n
        o = pd.DataFrame({"target": obs}); o.index = [f"c{i}" for i in range(len(o))]
        return ad.AnnData(X=sp.csr_matrix(np.vstack(blocks)), obs=o, var=real.var.copy())

    from cell_eval import MetricsEvaluator
    rows = []
    for name in ["ctrl_mean"] + list(deltas):
        cache = out / (name.replace(" ", "_").replace(",", "") + ".csv")
        if cache.exists():   # resume an interrupted run
            rows.append(pd.read_csv(cache)); continue
        pred = ctrl_mean_pred() if name == "ctrl_mean" else build(deltas[name])
        print(f"== {name}", flush=True)
        ev = MetricsEvaluator(adata_pred=pred, adata_real=real.copy(), control_pert=xb.CTRL, pert_col="target",
                              outdir=str(out / name.replace(" ", "_").replace(",", "")), num_threads=a.threads)
        res, _ = ev.compute(profile="vcc", write_csv=False)
        res = res.to_pandas(); res["method"] = name; rows.append(res)
        res.to_csv(cache, index=False)
        print(res.drop(columns=["method"]).mean(numeric_only=True).round(4).to_string(), flush=True)
    df = pd.concat(rows, ignore_index=True)
    df.to_csv(out / "per_pert_metrics.csv", index=False)
    summ = df.groupby("method").mean(numeric_only=True)
    summ.to_csv(out / "summary.csv")
    print(summ.round(4).to_string())


if __name__ == "__main__":
    main()
