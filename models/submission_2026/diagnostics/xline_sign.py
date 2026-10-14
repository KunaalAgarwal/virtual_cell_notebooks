"""Cross-line transfer of the shared knockdown response: sign accuracy of source line S's mean
response on target line L's individual knockdowns, over the top-N genes by |shared_S|."""
import numpy as np, pandas as pd, anndata as ad
PF = "/mnt/project-files"
def zmat(path):
    a = ad.read_h5ad(path)
    X = np.nan_to_num(np.asarray(a.X, dtype=np.float32), posinf=0, neginf=0).clip(-10, 10)
    tg = pd.Index(a.obs_names.str.split("_").str[1])
    keep = ~tg.str.contains("non-targeting")
    return pd.DataFrame(X[keep], columns=a.var["gene_name"].astype(str).values).T.groupby(level=0).mean().T
L = {"K562": zmat(f"{PF}/public-screens/data/K562_gwps_normalized_bulk_01.h5ad"),
     "RPE1": zmat(f"{PF}/public-screens/data/rpe1_normalized_bulk_01.h5ad"),
     "HepG2": zmat(f"{PF}/public-screens/data/hepg2_bulk.h5ad"),
     "Jurkat": zmat(f"{PF}/public-screens/data/jurkat_bulk.h5ad")}
z = np.load(f"{PF}/logo-shap/v2/h1/h1_pseudobulk.npz", allow_pickle=True)
L["H1"] = pd.DataFrame(z["delta"], columns=z["genes"].astype(str))
# strong knockdowns only in the big screens, to resemble a 300-target panel chosen for effects
for k in ("K562", "RPE1"):
    n = (L[k].abs() > 1).sum(1); L[k] = L[k][n >= n.quantile(0.8)]
shared = {k: v.mean(0) for k, v in L.items()}
rows = []
for tl, T in L.items():
    others = [k for k in L if k != tl]
    cands = {s: shared[s] for s in others}
    cands["avg(others)"] = pd.concat({s: shared[s] / shared[s].abs().median() for s in others}, axis=1).mean(1)
    cands["own (baseline)"] = shared[tl]
    for N in (1000, 4000, 8000):
        for s, v in cands.items():
            g = v.index.intersection(T.columns)
            top = v[g].abs().sort_values(ascending=False).index[:N]
            acc = (np.sign(T[top].to_numpy()) == np.sign(v[top].to_numpy())).mean()
            rows.append(dict(target_line=tl, source=s, N=N, acc=acc))
r = pd.DataFrame(rows)
print(r[r.N == 4000].pivot_table(index="source", columns="target_line", values="acc").round(3).to_string())
print(r[r.N == 8000].pivot_table(index="source", columns="target_line", values="acc").round(3).to_string())
