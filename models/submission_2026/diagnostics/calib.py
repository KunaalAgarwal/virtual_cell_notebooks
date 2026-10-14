"""Calibrate how big a per-gene fold must be for the official vcc2026 DE test (wilcoxon,
CPM+log1p, per-target BH) to call it, when cells are made the way build_submission does."""
from pathlib import Path
MODELS = str(Path(__file__).resolve().parents[2])  # repo models/
import sys, time, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from cell_eval2.de_compute import compute_de
sys.path.insert(0, MODELS + "/submission_2026")
from build_submission import make_cells
PF = "/mnt/project-files/vcc2026-data"
ctx = sys.argv[1] if len(sys.argv) > 1 else "A"
c = ad.read_h5ad(f"{PF}/context_{ctx}.h5ad"); X = c.X.tocsr()
lib = np.asarray(X.sum(1)).ravel()
N = X.multiply(1e6 / lib[:, None]).tocsr()            # CPM
mu = np.asarray(N.mean(0)).ravel()
m2 = np.asarray(N.multiply(N).mean(0)).ravel()
cv = np.sqrt(np.maximum(m2 - mu**2, 0)) / np.maximum(mu, 1e-12)
det = np.asarray((X > 0).mean(0)).ravel()
rng = np.random.default_rng(0)
GR = dict(rate=np.asarray(X.mean(0)).ravel(), mean_lib=lib.mean())
G = X.shape[1]
kappas = [3, 4, 6]
blocks, labels, truth = [], [], {}
sgn_g = rng.choice([-1, 1], size=G)
kappas = [4, 6, 8]
for kap in kappas:
    for rep in (0,):
        e = np.clip(kap * cv / np.sqrt(400), 0, 0.9)
        f = np.exp(sgn_g * np.log1p(e)).astype(np.float32)
        name = f"k{kap}"
        blocks.append(make_cells(X, f, rng, GR)); labels += [name] * 400
blocks.append(X); labels += ["non-targeting"] * X.shape[0]
A = ad.AnnData(sp.vstack(blocks).tocsr(), obs=pd.DataFrame({"target": labels}, index=[f"c{i}" for i in range(len(labels))]),
               var=pd.DataFrame(index=c.var_names))
t = time.time()
de = compute_de(A, backend="pdex", groupby="target", reference="non-targeting", mean_calc="arithmetic",
                epsilon=1e-9, input_type="counts", target_sum=1e6, filter_gene_min_cpm_cell=5.0, fdr_scope="per_pert")
print("DE secs", round(time.time() - t), de.columns, de.shape)
d = de.to_pandas()
d.to_csv(f"calib_de_{ctx}.csv.gz", index=False)
g = pd.DataFrame({"feature": c.var_names, "mu": mu, "cv": cv, "det": det})
d = d.merge(g, on="feature")
d["mu_bin"] = pd.qcut(np.log10(d.mu + 1e-3), 6)
print("genes tested per target", d.groupby("target").size().iloc[0])
d["sgn"] = pd.Series(sgn_g, index=c.var_names).reindex(d.feature).to_numpy()
d["right"] = (d.p_adj < 0.05) & (np.sign(d.log2_fold_change) == d.sgn)
d["e"] = 0
for col in ("right",):
    print(d.pivot_table(index="mu_bin", columns=["target", "sgn"], values=col, aggfunc="mean").round(2).to_string())
print(d.assign(sig=d.p_adj<0.05).groupby("target")[["sig","right"]].mean())
