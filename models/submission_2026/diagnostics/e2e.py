"""End to end on one context: c1_calibrated vs c1_wide / c4_wide cells for 12 targets, official DE."""
from pathlib import Path
MODELS = str(Path(__file__).resolve().parents[2])  # repo models/
import sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
R = MODELS
sys.path.insert(0, f"{R}/submission_2026"); sys.path.insert(0, f"{R}/public_screens")
import build_submission as bs, xline_benchmark as xb
from cell_eval2.de_compute import compute_de
PF = "/mnt/project-files/vcc2026-data"; ctx = sys.argv[1]
genes = pd.read_csv(f"{PF}/gene_names.csv")["gene_name"].astype(str).tolist()
targets = pd.read_csv(f"{PF}/pert_counts.csv")["target_gene"].astype(str).tolist()
z = np.load(f"{R}/submission_2026/data/h1_pseudobulk.npz", allow_pickle=True)
col = pd.Series(np.arange(len(z["genes"])), index=z["genes"].astype(str)); pos = col.reindex(genes).to_numpy(); have = np.isfinite(pos)
Dm = np.zeros((len(z["perts"]), len(genes)), np.float32); Dm[:, have] = z["delta"][:, pos[have].astype(int)]
h1 = dict(delta=Dm, perts=z["perts"].astype(str))
k562 = xb.load_k562("/mnt/project-files/public-screens/data/K562_gwps_normalized_bulk_01.h5ad")
zs = np.load(f"{R}/submission_2026/data/shared_response.npz", allow_pickle=True)
shared = pd.Series(zs["score"], index=zs["genes"].astype(str)).reindex(genes).fillna(0).to_numpy()
X = ad.read_h5ad(f"{PF}/context_{ctx}.h5ad").X.tocsr()
lib = np.asarray(X.sum(1)).ravel(); GR = dict(rate=np.asarray(X.mean(0)).ravel(), mean_lib=lib.mean())
rng = np.random.default_rng(0); pick = list(rng.choice(len(targets), 12, replace=False))
blocks, labels, deltas = [], [], {}
for preset in ("c1_wide", "c4_wide"):
    bs.PARAMS = bs.PRESETS[preset]
    fold, delta, info = bs.predict_folds(ctx, X, genes, targets, h1, k562, xb, shared)
    gr = GR if bs.PARAMS.get("up_mode") == "depth" else None
    for i in pick:
        blocks.append(bs.make_cells(X, fold[i], rng, gr)); labels += [f"{preset}|{targets[i]}"] * 400
        deltas[f"{preset}|{targets[i]}"] = np.log(fold[i])
    print(preset, "median genes moved >10%:", int(info.n_genes_changed_10pct.median()), flush=True)
blocks.append(X); labels += ["non-targeting"] * X.shape[0]
A = ad.AnnData(sp.vstack(blocks).tocsr(), obs=pd.DataFrame({"target": labels}, index=[f"c{i}" for i in range(len(labels))]), var=pd.DataFrame(index=genes))
print("nnz per cell pred vs ctrl:", sp.vstack(blocks[:-1]).nnz / (len(labels) - X.shape[0]), X.nnz / X.shape[0])
de = compute_de(A, backend="pdex", groupby="target", reference="non-targeting", mean_calc="arithmetic", epsilon=1e-9,
                input_type="counts", target_sum=1e6, filter_gene_min_cpm_cell=5.0, fdr_scope="per_pert").to_pandas()
gi = {g: i for i, g in enumerate(genes)}
de["want"] = [np.sign(deltas[t][gi[f]]) for t, f in zip(de.target, de.feature)]
de["preset"] = de.target.str.split("|").str[0]; de["sig"] = de.p_adj < 0.05
de["agree"] = de.sig & (np.sign(de.log2_fold_change) == de.want)
de["sh_agree"] = de.sig & (np.sign(de.log2_fold_change) == np.sign(shared[[gi[f] for f in de.feature]]))
s = de.groupby(["preset", "target"]).agg(n_sig=("sig", "sum"), agree=("agree", "sum"), shared_sign=("sh_agree", "sum"),
                                          med_abs_lfc_sig=("abs_log2_fold_change", lambda v: np.nan))
s["med_abs_lfc_sig"] = de[de.sig].groupby(["preset", "target"]).abs_log2_fold_change.median()
print(s.groupby("preset").median().round(2).to_string())
de.to_csv(f"e2e_de_{ctx}.csv.gz", index=False)
