"""Does adding a common shared shift hurt cosine discrimination? H1 pseudobulk proxy:
pred_t = signal_t (true H1 delta + noise, tuned to pds ~0.55) + lam * w (wide shift, cross-line sign)."""
from pathlib import Path
MODELS = str(Path(__file__).resolve().parents[2])  # repo models/
import numpy as np, pandas as pd, anndata as ad
PF = "/mnt/project-files"
z = np.load(f"{PF}/logo-shap/v2/h1/h1_pseudobulk.npz", allow_pickle=True)
G = z["genes"].astype(str); D = z["delta"].astype(np.float64); cm = z["ctrl_mean"]
sr = np.load(MODELS + "/submission_2026/data/shared_response.npz", allow_pickle=True)
# leave H1 out of the shared vector is not possible from the file; rebuild quickly without H1
import importlib.util, sys
sys.argv = ["x"]
from importlib import import_module
spec = importlib.util.spec_from_file_location("srm", MODELS + "/submission_2026/shared_response.py")
srm = importlib.util.module_from_spec(spec); spec.loader.exec_module(srm)
P = f"{PF}/public-screens/data"
src = {"K562": srm.zmean(f"{P}/K562_gwps_normalized_bulk_01.h5ad", .2), "RPE1": srm.zmean(f"{P}/rpe1_normalized_bulk_01.h5ad", .2),
       "HepG2": srm.zmean(f"{P}/hepg2_bulk.h5ad"), "Jurkat": srm.zmean(f"{P}/jurkat_bulk.h5ad")}
s = pd.concat({k: v / v.abs().median() for k, v in src.items()}, axis=1).reindex(G).mean(1).fillna(0).to_numpy()
# context A CV per gene for the detectable size
c = ad.read_h5ad(f"{PF}/vcc2026-data/context_A.h5ad"); X = c.X.tocsr(); lib = np.asarray(X.sum(1)).ravel()
N = X.multiply(1e6 / lib[:, None]).tocsr(); mu = np.asarray(N.mean(0)).ravel()
sd = np.sqrt(np.maximum(np.asarray(N.multiply(N).mean(0)).ravel() - mu**2, 0))
e_all = pd.Series(6 * sd / np.maximum(mu, 1e-12) / 20, index=c.var_names).reindex(G).to_numpy()
mu_s = pd.Series(mu, index=c.var_names).reindex(G).to_numpy()
ok = np.isfinite(e_all) & (mu_s >= 5) & (e_all <= .35) & (s != 0)
pick = np.argsort(-np.where(ok, np.abs(s), -1))[:8000]
m = np.expm1(cm)  # CP10k mean
w = np.zeros(len(G)); w[pick] = np.log1p(m[pick] * (1 + np.sign(s[pick]) * e_all[pick])) - np.log1p(m[pick])
w = np.nan_to_num(w)
def pds(P, R):
    Pn = P / np.linalg.norm(P, axis=1, keepdims=True); Rn = R / np.linalg.norm(R, axis=1, keepdims=True)
    S = Pn @ Rn.T; n = len(P)
    rank = (S > np.diag(S)[:, None]).sum(1)
    return (1 - rank / (n - 1)).mean()
rng = np.random.default_rng(0)
R = D
for noise in (150.0, 300.0):
    sig = D + noise * np.linalg.norm(D, axis=1, keepdims=True) / np.sqrt(D.shape[1]) * rng.standard_normal(D.shape)
    print("norm signal median", np.median(np.linalg.norm(D, axis=1)).round(2), "norm w", np.linalg.norm(w).round(2))
    for amp in (1, 4, 8):
        for lam in (0, 0.25, 0.5, 1):
            print(f"noise {noise} specific x{amp} lam {lam}: pds {pds(amp*sig + lam*w, R):.3f}")
