"""Approximate the vcc2026 DE metrics (fidelity_yield_raw, sig_jaccard, reach_raw, lfc_nmae)
from cell-eval 0.8 style DE tables (target, feature, log2_fold_change, fdr, p_value)."""
import sys, numpy as np, pandas as pd

def load(p):
    d = pd.read_csv(p, usecols=["target", "feature", "log2_fold_change", "p_value", "fdr"])
    return d

def metrics(pred, real, alpha=0.05, floor=0.9):
    m = real.merge(pred, on=["target", "feature"], suffixes=("_r", "_p"))
    m = m[m.target != m.feature]
    out = []
    for t, g in m.groupby("target", sort=False):
        lr, lp = g.log2_fold_change_r.to_numpy(), g.log2_fold_change_p.to_numpy()
        sr, sp = (g.fdr_r < alpha).to_numpy(), (g.fdr_p < alpha).to_numpy()
        adj = np.isfinite(lr) & (lr != 0)
        com = np.isfinite(lp) & (lp != 0)
        match = adj & com & (np.sign(lp) == np.sign(lr))
        n_conf = int(sr.sum())
        scored = sp & adj
        n_pred, k = int(scored.sum()), int((scored & match).sum())
        fid = k / max(n_pred, n_conf) if max(n_pred, n_conf) else np.nan
        inter, uni = (sr & sp).sum(), (sr | sp).sum()
        jac = inter / uni if uni else 1.0
        # reach: rank real-sig adjudicable genes by pred (fdr, p, -|lfc|)
        pool = g[sr & adj].assign(_m=match[sr & adj], _a=-np.abs(np.nan_to_num(lp[sr & adj], nan=-np.inf)))
        pool = pool.sort_values(["fdr_p", "p_value_p", "_a"])
        if n_conf:
            pur = np.cumsum(pool._m.to_numpy()) / np.arange(1, len(pool) + 1)
            ok = np.nonzero(pur >= floor)[0]
            reach = (ok[-1] + 1) / n_conf if len(ok) else 0.0
        else:
            reach = np.nan
        gate = sr & np.isfinite(lr) & np.isfinite(lp)
        nmae = (np.abs(lp[gate] - lr[gate]).mean() / np.abs(lr[gate]).mean()) if gate.sum() >= 10 else np.nan
        acc = match[scored].mean() if n_pred else np.nan
        out.append(dict(target=t, n_conf=n_conf, n_pred=n_pred, sign_acc_called=acc, fid=fid, jac=jac, reach=reach, nmae=nmae))
    return pd.DataFrame(out)

if __name__ == "__main__":
    rows = []
    for d in sys.argv[1:]:
        r = metrics(load(f"{d}/pred_de.csv"), load(f"{d}/real_de.csv"))
        s = r.drop(columns="target").median().add_prefix("med_")
        s2 = r[["fid", "jac", "reach", "nmae", "sign_acc_called"]].mean()
        rows.append(pd.concat([pd.Series({"method": d}), s2, s[["med_n_conf", "med_n_pred"]]]))
        print(rows[-1].to_dict(), flush=True)
    pd.DataFrame(rows).to_csv("de_metrics_out.csv", index=False)
