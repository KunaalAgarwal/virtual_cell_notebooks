"""If a prediction called the top-N genes (ranked by |shared response|) significant with the
shared response's sign, what would direction_fidelity_yield_raw be on real H1 knockdowns?"""
import numpy as np, pandas as pd, anndata as ad
PF = "/mnt/project-files"
real = pd.read_csv(f"{PF}/logo-shap/v2/h1/celleval/mean_effect/real_de.csv",
                   usecols=["target", "feature", "log2_fold_change", "fdr"])
real = real[real.target != real.feature]
z = np.load(f"{PF}/logo-shap/v2/h1/h1_pseudobulk.npz", allow_pickle=True)
genes = z["genes"].astype(str); perts = z["perts"].astype(str); D = z["delta"]

def zmean(path, ctrl_pat="non-targeting"):
    a = ad.read_h5ad(path)
    X = np.nan_to_num(np.asarray(a.X, dtype=np.float32), posinf=0, neginf=0).clip(-10, 10)
    tg = pd.Index(a.obs_names.str.split("_").str[1])
    X = X[~tg.str.contains(ctrl_pat)]
    return pd.Series(X.mean(0), index=a.var["gene_name"].astype(str).values).groupby(level=0).mean()

src = {"K562": zmean(f"{PF}/public-screens/data/K562_gwps_normalized_bulk_01.h5ad"),
       "RPE1": zmean(f"{PF}/public-screens/data/rpe1_normalized_bulk_01.h5ad"),
       "HepG2": zmean(f"{PF}/public-screens/data/hepg2_bulk.h5ad"),
       "Jurkat": zmean(f"{PF}/public-screens/data/jurkat_bulk.h5ad")}
def rz(s): return s / s.abs().median()
allg = sorted(set().union(*[s.index for s in src.values()]))
X = pd.concat({k: rz(v) for k, v in src.items()}, axis=1).reindex(allg)
src["cross-line avg"] = X.mean(1, skipna=True)
src["cross-line sign vote"] = np.sign(X).mean(1, skipna=True) * X.abs().mean(1)
gi = pd.Series(np.arange(len(genes)), index=genes)
targets = real.target.unique()
Ns = [100, 300, 1000, 2000, 4000, 8000]
rows = []
for t in targets:
    g = real[real.target == t]
    lr = g.log2_fold_change.to_numpy(); feat = g.feature.to_numpy()
    sig = (g.fdr < 0.05).to_numpy(); adj = np.isfinite(lr) & (lr != 0)
    n_conf = sig.sum()
    cands = {}
    if t in set(perts):
        loo = D[perts != t].mean(0)
        cands["H1 own-line mean (≈ baseline)"] = pd.Series(loo, index=genes)
    for k, v in src.items():
        cands[k] = v
    for name, s in cands.items():
        sv = s.reindex(feat).to_numpy()
        ok = np.isfinite(sv) & (sv != 0) & adj
        order = np.argsort(-np.abs(np.where(ok, sv, 0)))
        for N in Ns:
            pick = order[:N]; pick = pick[ok[pick]]
            kk = (np.sign(sv[pick]) == np.sign(lr[pick])).sum()
            rows.append(dict(target=t, src=name, N=N, n_conf=n_conf,
                             fid=kk / max(len(pick), n_conf), acc=kk / max(len(pick), 1),
                             jac=(sig[pick].sum()) / (n_conf + len(pick) - sig[pick].sum())))
r = pd.DataFrame(rows)
print("median N_conf", r.groupby("target").n_conf.first().median())
for m in ("fid", "acc", "jac"):
    print(m); print(r.pivot_table(index="src", columns="N", values=m, aggfunc="mean").round(3).to_string())
r.to_csv("sim_topn_h1.csv", index=False)
