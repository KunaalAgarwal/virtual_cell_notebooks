"""How far can the same gene's profile in another line push cosine discrimination and the
mse ratio, when the truth is clean? Truth: RPE1 (Replogle 2022) knockdowns; source: K562 GWPS.
Both converted from z to approximate lfc = log1p(z * sd / mean) per line."""
import numpy as np, pandas as pd, anndata as ad
D = "/mnt/project-files/public-screens/data"
def load(path):
    a = ad.read_h5ad(path)
    Z = np.nan_to_num(np.asarray(a.X, dtype=np.float64), posinf=0, neginf=0).clip(-10, 10)
    tg = a.obs_names.str.split("_").str[1].to_numpy()
    keep = ~pd.Series(tg).str.contains("non-targeting").to_numpy()
    df = pd.DataFrame(Z[keep], index=tg[keep], columns=a.var["gene_name"].astype(str).values).groupby(level=0).mean()
    df = df.T.groupby(level=0).mean().T
    cv = pd.Series((a.var["std"] / a.var["mean"].clip(lower=1e-3)).to_numpy(), index=a.var["gene_name"].astype(str).values).groupby(level=0).mean()
    mean = pd.Series(a.var["mean"].to_numpy(), index=a.var["gene_name"].astype(str).values).groupby(level=0).mean()
    return df, cv, mean
K, cvK, mK = load(f"{D}/K562_gwps_normalized_bulk_01.h5ad")
R, cvR, mR = load(f"{D}/rpe1_normalized_bulk_01.h5ad")
genes = K.columns.intersection(R.columns)
common = R.index.intersection(K.index)
# a 300-target panel with effects, like vcc2026-val-1 (55% of targets have a K562 effect)
nR = (R.loc[common, genes].abs() > 2).sum(1)
panel = nR.sort_values(ascending=False).index[:300]
feat = genes.difference(panel)                      # scorer: all panel target genes excluded
def lfc(df, cv): return np.log1p(np.clip(df.to_numpy() * cv.reindex(df.columns).to_numpy(), -0.95, 20))
truth = lfc(R.loc[panel, feat], cvR); src = lfc(K.loc[panel, feat], cvK)
shared_R = lfc(R.loc[R.index.difference(panel), feat], cvR).mean(0)      # RPE1's own shared response (other KDs)
shared_K = lfc(K.loc[K.index.difference(panel), feat], cvK).mean(0)
def pds(P, T):
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12); Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    S = Pn @ Tn.T; d = np.diag(S)[:, None]
    return (1 - ((S > d).sum(1) + ((S == d).sum(1) - 1) / 2) / (len(S) - 1)).mean()
def mse_ratio(P, T): return ((P - T) ** 2).sum() / (T ** 2).sum()
def var_expl_shared(T, s):
    b = (T @ s) / (s @ s); return 1 - ((T - np.outer(b, s)) ** 2).sum() / (T ** 2).sum(), np.median(b)
rows = []
def add(name, P): rows.append((name, round(pds(P, truth), 3), round(mse_ratio(P, truth), 3)))
add("control (zero)", np.zeros_like(truth))
add("RPE1 own shared mean (oracle shared)", np.tile(shared_R, (300, 1)))
add("K562 shared mean, ported", np.tile(shared_K, (300, 1)))
for k in (1.0, 0.5, 0.3):
    add(f"K562 shared x{k}", np.tile(k * shared_K, (300, 1)))
add("K562 same gene, raw", src)
for q in (0.5, 0.3):
    add(f"K562 same gene x{q}", q * src)
# sparsify: keep each target's top-k genes by |z| in K562
zK = K.loc[panel, feat].to_numpy()
for k in (200, 1000):
    thr = -np.sort(-np.abs(zK), axis=1)[:, k - 1][:, None]
    add(f"K562 same gene, top{k} genes only", np.where(np.abs(zK) >= thr, src, 0))
add("K562 same gene x0.5 + K562 shared x0.5", 0.5 * src + 0.5 * shared_K)
add("K562 same gene x0.5 + RPE1 shared (oracle)", 0.5 * src + shared_R)
add("oracle (truth)", truth)
# what fraction of truth is the shared response?
ve, b = var_expl_shared(truth, shared_R); veK, bK = var_expl_shared(truth, shared_K)
print(f"variance of RPE1 panel responses explained by RPE1's own shared mean: {ve:.2f} (median coef {b:.2f}); by K562's shared mean: {veK:.2f} (coef {bK:.2f})")
print("cosine(shared_R, shared_K) =", round(float(shared_R @ shared_K / np.linalg.norm(shared_R) / np.linalg.norm(shared_K)), 3))
print("median norm truth / src / shared_R:", np.median(np.linalg.norm(truth, axis=1)).round(2), np.median(np.linalg.norm(src, axis=1)).round(2), np.linalg.norm(shared_R).round(2))
print(pd.DataFrame(rows, columns=["prediction", "pds_cosine", "mse_ratio"]).to_string(index=False))
