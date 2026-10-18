"""Run our actual c1 pipeline (xline_benchmark.build_features/predict) with RPE1 as the target
line, truth = RPE1 knockdown profiles. Panels of different difficulty."""
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "public_screens")); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import xline_benchmark as xb, build_submission as bs
exec(open("xline_pds.py").read().split("# a 300-target panel")[0])
def lfc(df, cv): return np.log1p(np.clip(df.to_numpy() * cv.reindex(df.columns).to_numpy(), -0.95, 20))
def pds(P, T):
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12); Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    S = Pn @ Tn.T; d = np.diag(S)[:, None]
    return (1 - ((S > d).sum(1) + ((S == d).sum(1) - 1) / 2) / (len(S) - 1)).mean()
def mse_ratio(P, T): return ((P - T) ** 2).sum() / (T ** 2).sum()

common = common.intersection(genes)
k562 = xb.load_k562(f"{D}/K562_gwps_normalized_bulk_01.h5ad")
z = np.load(str(Path(__file__).resolve().parents[1] / "data/h1_pseudobulk.npz"), allow_pickle=True)
rg = np.array(genes)
pos = pd.Series(np.arange(len(z["genes"])), index=z["genes"].astype(str)).reindex(rg).to_numpy(); have = np.isfinite(pos)
Dm = np.zeros((len(z["perts"]), len(rg)), np.float32); Dm[:, have] = z["delta"][:, pos[have].astype(int)]
h1 = dict(delta=Dm, perts=z["perts"].astype(str))
meta = k562["meta"]
rng = np.random.default_rng(0)
nR = (R.loc[common, genes].abs() > 2).sum(1).sort_values(ascending=False)
sigK = (meta["energy_test_p_value"].reindex(common) < 0.05)
panels = {"strongest 300 in RPE1": nR.index[:300],
          "random 300": rng.choice(common, 300, replace=False),
          "random 300 with 55% K562-significant (like vcc2026)": np.concatenate([rng.choice(common[sigK.reindex(common).fillna(False).to_numpy()], 165, replace=False),
                                                                                 rng.choice(common[~sigK.reindex(common).fillna(False).to_numpy()], 135, replace=False)])}
rows = []
for pname, panel in panels.items():
    panel = pd.Index(panel); feat = genes.difference(panel)
    truth = lfc(R.loc[panel, feat], cvR)
    line = dict(name="RPE1", genes=rg, perts=list(panel), ctrl_mean=np.log1p(mR.reindex(rg).fillna(0).to_numpy() / 10), det=np.ones(len(rg)))
    F = xb.build_features(line, h1, k562)
    fi = np.isin(rg, feat)
    def P(p):
        return xb.predict(line, F, p)[:, fi]
    base = dict(bs.PRESETS["c1_calibrated"]); base["gate"] = 0
    variants = {"c1_calibrated (shared H1 + borrowed H1 + K562 direct)": base,
                "K562 direct only (shared none, C 0)": {**base, "shared": "none", "C": 0.0},
                "shared H1 + K562 direct (no borrowed)": {**base, "C": 0.0},
                "borrowed H1 + K562 direct (shared none)": {**base, "shared": "none"},
                "c4: shared H1 + 4x borrowed + direct": {**base, "C": 4.0}}
    for vn, p in variants.items():
        pr = P(p); rows.append((pname, vn, round(pds(pr, truth), 3), round(mse_ratio(pr, truth), 3), round(float(np.median(np.linalg.norm(pr, axis=1))), 2)))
    src = lfc(K.loc[panel, feat], cvK); rows.append((pname, "same gene K562, plain lfc conversion", round(pds(src, truth), 3), round(mse_ratio(src, truth), 3), round(float(np.median(np.linalg.norm(src, axis=1))), 2)))
    rows.append((pname, "truth norm (median)", None, None, round(float(np.median(np.linalg.norm(truth, axis=1))), 2)))
pd.set_option("display.width", 200)
print(pd.DataFrame(rows, columns=["panel", "prediction", "pds_cosine", "mse_ratio", "median norm"]).to_string(index=False))
