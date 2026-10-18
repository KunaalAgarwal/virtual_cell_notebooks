import numpy as np, pandas as pd
exec(open("pipeline_on_rpe1.py").read().split("k562 = xb.load_k562")[0])
rng = np.random.default_rng(0)
sigK = pd.read_pickle if False else None
panel = pd.Index(rng.choice(common, 300, replace=False)); feat = genes.difference(panel)
truth = lfc(R.loc[panel, feat], cvR); src = lfc(K.loc[panel, feat], cvK); zR = R.loc[panel, feat].to_numpy(); zK = K.loc[panel, feat].to_numpy()
sig = np.abs(zR) > 2          # stand-in for the reference's significant genes
def nmae(P): 
    v = [np.abs(P[i, sig[i]] - truth[i, sig[i]]).mean() / np.abs(truth[i, sig[i]]).mean() for i in range(len(P)) if sig[i].sum() >= 10]
    return float(np.mean(v))
def reach(P):
    out = []
    for i in range(len(P)):
        g = np.flatnonzero(sig[i])
        if len(g) < 5: continue
        order = g[np.argsort(-np.abs(P[i, g]))]
        ok = np.sign(P[i, order]) == np.sign(truth[i, order]); pur = np.cumsum(ok) / np.arange(1, len(ok) + 1)
        hit = np.flatnonzero(pur >= 0.9); out.append((hit[-1] + 1) / len(g) if len(hit) else 0)
    return float(np.mean(out))
rows = []
for s in (1, 2, 3, 4, 6):
    rows.append((f"K562 direct x{s}", round(pds(s * src, truth), 3), round(mse_ratio(s * src, truth), 3), round(nmae(s * src), 3), round(reach(s * src), 3)))
# shrink per gene by K562 confidence: keep only |z|>1 genes
for s in (2, 3):
    P = np.where(np.abs(zK) > 1, s * src, 0); rows.append((f"K562 direct x{s}, |z|>1 genes only", round(pds(P, truth), 3), round(mse_ratio(P, truth), 3), round(nmae(P), 3), round(reach(P), 3)))
# wide call set with per-target deviation signs (c1_dev idea) on top of direct x3
e = 0.2 * np.ones(len(feat)); dev = src - src.mean(0)
for q in (0.5, 1.0):
    w = np.sign(dev) * e * q; w[np.abs(dev) < 1e-6] = 0
    P = np.where(np.abs(3 * src) >= np.abs(w), 3 * src, w); rows.append((f"direct x3 + wide dev-sign calls at {0.2*q:.1f} log", round(pds(P, truth), 3), round(mse_ratio(P, truth), 3), round(nmae(P), 3), round(reach(P), 3)))
rows.append(("truth", 1.0, 0.0, 0.0, round(reach(truth), 3)))
print(pd.DataFrame(rows, columns=["prediction", "pds", "mse_ratio", "nmae", "reach"]).to_string(index=False))
