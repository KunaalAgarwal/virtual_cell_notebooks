import numpy as np, pandas as pd
exec(open("xline_pds.py").read().split("rows = []")[0])
rng = np.random.default_rng(0)
zK = K.loc[panel, feat].to_numpy()
# per-target strength: oracle coefficient b_t of RPE1 shared, and a K562-based predictor
bt = (truth @ shared_R) / (shared_R @ shared_R)
sK = np.abs(zK).sum(1); sK = sK / np.median(sK)
print("corr(oracle shared coef, K562 strength) =", round(np.corrcoef(bt, sK)[0, 1], 3), "; corr with truth norm:", round(np.corrcoef(np.linalg.norm(truth, axis=1), sK)[0, 1], 3))
# gene-wise sd of a pseudobulk delta (noise), proportional to 1/sqrt(mean expression) in K562 units
gsd = 1 / np.sqrt(np.maximum(mK.reindex(feat).to_numpy(), 1e-3)); gsd = gsd / np.linalg.norm(gsd)
rows = []
def add(name, P, noise=0.0):
    N = noise * gsd * rng.standard_normal(P.shape) if noise else 0
    rows.append((name, noise, round(pds(P + N, truth), 3), round(mse_ratio(P + N, truth), 3)))
for nz in (0, 3.6, 7):
    add("K562 same gene (x1)", src, nz)
    add("K562 same gene x3", 3 * src, nz)
    add("K562 same gene x3 + RPE1 shared x median b (oracle shared, one size)", 3 * src + np.median(bt) * shared_R, nz)
    add("K562 same gene x3 + RPE1 shared x b_t (oracle per-target size)", 3 * src + bt[:, None] * shared_R, nz)
    add("K562 same gene x3 + RPE1 shared x 2.5*strength_K562", 3 * src + 2.5 * sK[:, None] * shared_R, nz)
    add("K562 same gene x3 + K562 shared x 5*strength_K562", 3 * src + 5 * sK[:, None] * shared_K, nz)
print(pd.DataFrame(rows, columns=["prediction", "pred noise norm", "pds_cosine", "mse_ratio"]).to_string(index=False))
