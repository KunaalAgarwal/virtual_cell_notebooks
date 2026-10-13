"""Score candidate submissions on real H1 knockdowns with Arc's official cell-eval metrics.

Every method builds a prediction the way a real submission would: sample 400 control
cells per knockdown, shift each cell by a predicted change in log1p(CP10k) space, and
convert back to integer counts. Methods:

  ctrl_sampled   resampled control cells, no change (the "nothing happens" baseline)
  ctrl_mean      every cell = rounded control mean (the repo's baseline_control_mean)
  mean_effect    control cells + the leave-one-out mean effect of the other knockdowns
  network        control cells + the control-cell GRN's propagated change (fixed self-drop)
  oracle         control cells + the knockdown's real pseudobulk change (an upper bound
                 for "shift the controls" predictions; not achievable zero-shot)

The real data is the test .h5ad from h1_extract_perts.py. cell-eval normalizes raw
counts itself and computes DE with pdex. Metrics in the "vcc" profile: mae (lower is
better), discrimination_score_l1 and overlap_at_N (higher is better).

    python celleval_h1.py --real h1_real_test.h5ad --pseudobulk h1_pseudobulk.npz \
        --grn grn_h1_pc0.npz --out celleval_results/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from eval_h1 import propagate  # noqa: E402

SELF_DROP = -0.633  # H1 median change of the knocked-down gene itself, log1p(CP10k)


def to_lognorm(X):
    X = sp.csr_matrix(X, dtype=np.float32)
    lib = np.asarray(X.sum(1)).ravel()
    Xn = sp.diags((1e4 / np.maximum(lib, 1)).astype(np.float32)) @ X
    Xn.data = np.log1p(Xn.data)
    return Xn, lib


def shifted_counts(ctrl_log, lib, delta):
    """Apply a per-gene log change to control cells and return integer counts."""
    out = np.asarray(ctrl_log.todense(), dtype=np.float32)
    out += delta[None, :].astype(np.float32)
    np.clip(out, 0, None, out=out)
    counts = np.expm1(out) * (lib / 1e4)[:, None]
    return sp.csr_matrix(np.rint(counts).astype(np.float32))


def build_pred(real, deltas, cells_per_pert, rng):
    """deltas: dict pert -> change vector over real.var_names (zeros for 'no change')."""
    ctrl_mask = (real.obs["target"] == "non-targeting").to_numpy()
    ctrl = real[ctrl_mask]
    ctrl_log, lib = to_lognorm(ctrl.X)
    blocks, obs = [ctrl.X.copy()], [pd.DataFrame({"target": ["non-targeting"] * ctrl.n_obs})]
    for pert, d in deltas.items():
        pick = rng.choice(ctrl.n_obs, cells_per_pert, replace=cells_per_pert > ctrl.n_obs)
        blocks.append(shifted_counts(ctrl_log[pick], lib[pick], d))
        obs.append(pd.DataFrame({"target": [pert] * cells_per_pert}))
    obs = pd.concat(obs, ignore_index=True)
    obs.index = [f"p{i}" for i in range(len(obs))]
    return ad.AnnData(X=sp.vstack(blocks).tocsr(), obs=obs, var=real.var.copy())


def control_mean_pred(real, perts, cells_per_pert):
    ctrl = real[(real.obs["target"] == "non-targeting").to_numpy()]
    mean = np.rint(np.asarray(ctrl.X.mean(0)).ravel()).astype(np.float32)
    row = sp.csr_matrix(mean)
    blocks = [ctrl.X.copy()] + [sp.vstack([row] * cells_per_pert) for _ in perts]
    obs = pd.DataFrame({"target": ["non-targeting"] * ctrl.n_obs + [p for p in perts for _ in range(cells_per_pert)]})
    obs.index = [f"p{i}" for i in range(len(obs))]
    return ad.AnnData(X=sp.vstack(blocks).tocsr(), obs=obs, var=real.var.copy())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--real", required=True)
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--grn", default=None, help="grn_ridge .npz for the network method (skip if omitted)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cells-per-pert", type=int, default=400)
    ap.add_argument("--methods", nargs="+", default=["ctrl_sampled", "ctrl_mean", "mean_effect", "network", "oracle"])
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    real = ad.read_h5ad(a.real)
    real.obs["target"] = real.obs["target"].astype(str)
    perts = sorted(p for p in real.obs["target"].unique() if p != "non-targeting")
    genes = real.var_names.to_numpy()
    n_genes = len(genes)
    print(f"real: {real.n_obs} cells, {len(perts)} knockdowns, {n_genes} genes")

    pb = np.load(a.pseudobulk, allow_pickle=True)
    assert list(pb["genes"]) == list(genes), "gene order differs between pseudobulk and real file"
    pb_perts = list(pb["perts"])
    delta_all = pb["delta"].astype(np.float64)
    gi = {g: i for i, g in enumerate(genes)}

    deltas = {}
    zero = np.zeros(n_genes)
    deltas["ctrl_sampled"] = {p: zero for p in perts}
    # Leave-one-out mean effect over the *other* knockdowns, self gene set to the median drop.
    tot = delta_all.sum(0)
    me = {}
    for p in perts:
        i = pb_perts.index(p)
        d = (tot - delta_all[i]) / (len(pb_perts) - 1)
        d[gi[p]] = SELF_DROP
        me[p] = d
    deltas["mean_effect"] = me
    deltas["oracle"] = {p: delta_all[pb_perts.index(p)] for p in perts}

    if "network" in a.methods and a.grn:
        g = np.load(a.grn, allow_pickle=True)
        ngenes = list(g["genes"])
        W = g["coef_mean"]
        in_net = [p for p in perts if p in set(ngenes)]
        idx = np.array([ngenes.index(p) for p in in_net])
        D = propagate(W, idx, np.full(len(idx), SELF_DROP), 20)       # network genes x perts
        cols = np.array([gi[x] for x in ngenes])
        net = {}
        for c, p in enumerate(in_net):
            d = np.zeros(n_genes)
            d[cols] = D[:, c]
            net[p] = d
        for p in perts:
            if p not in net:  # not in the network: fall back to self-drop only
                d = np.zeros(n_genes); d[gi[p]] = SELF_DROP; net[p] = d
        deltas["network"] = net
        del W, g

    from cell_eval import MetricsEvaluator

    rows = []
    for m in a.methods:
        if m == "ctrl_mean":
            pred = control_mean_pred(real, perts, a.cells_per_pert)
        elif m in deltas:
            pred = build_pred(real, deltas[m], a.cells_per_pert, np.random.default_rng(a.seed))
        else:
            print(f"skipping {m}"); continue
        print(f"== {m}: pred {pred.n_obs} cells", flush=True)
        ev = MetricsEvaluator(adata_pred=pred, adata_real=real.copy(), control_pert="non-targeting",
                              pert_col="target", outdir=str(out / m), num_threads=4)
        res, agg = ev.compute(profile="vcc", write_csv=True)
        res = res.to_pandas()
        res["method"] = m
        rows.append(res)
        print(res.drop(columns=["method"]).describe().loc[["mean", "50%"]].round(4).to_string(), flush=True)

    df = pd.concat(rows, ignore_index=True)
    df.to_csv(out / "per_pert_metrics.csv", index=False)
    summ = df.groupby("method").mean(numeric_only=True)
    summ.to_csv(out / "summary.csv")
    print("\n" + summ.round(4).to_string())


if __name__ == "__main__":
    main()
