"""H1 check of wide-call variants: pds cosine, fidelity, nmae (approx) on real H1 knockdowns.
Predictions for H1 targets come from the same c1 machinery (K562 similarity + H1 effects of other
genes + K562 direct), so this is the closest local stand-in for the leaderboard."""
import sys, numpy as np, pandas as pd, anndata as ad
from pathlib import Path
R = str(Path(__file__).resolve().parents[2]); PF = "/mnt/project-files"  # repo models/, project data
sys.path.insert(0, f"{R}/public_screens"); sys.path.insert(0, f"{R}/submission_2026")
import xline_benchmark as xb, build_submission as bs
z = np.load(f"{R}/submission_2026/data/h1_pseudobulk.npz", allow_pickle=True)
G = z["genes"].astype(str); P = z["perts"].astype(str); D = z["delta"].astype(np.float64); cm = z["ctrl_mean"]
k562 = xb.load_k562(f"{PF}/public-screens/data/K562_gwps_normalized_bulk_01.h5ad")
T = [p for p in P if p in k562["Z"].index and p in set(G)]
line = dict(name="H1", genes=G, perts=T, ctrl_mean=cm, det=(np.expm1(cm) > 0.01).astype(float))
h1 = dict(delta=D, perts=P)
F = xb.build_features(line, h1, k562)
pp = dict(bs.PRESETS["c1_calibrated"]); pp["gate"] = 0.0
pred = xb.predict(line, F, pp)
Rl = D[[list(P).index(t) for t in T]]
# per-gene detectable size (kappa SE of a 400-cell CPM mean), borrowed from context A CVs
c = ad.read_h5ad(f"{PF}/vcc2026-data/context_A.h5ad"); X = c.X.tocsr(); lib = np.asarray(X.sum(1)).ravel()
N = X.multiply(1e6 / lib[:, None]).tocsr(); mu = np.asarray(N.mean(0)).ravel()
sd = np.sqrt(np.maximum(np.asarray(N.multiply(N).mean(0)).ravel() - mu ** 2, 0))
cvA = pd.Series(sd / np.maximum(mu, 1e-12), index=c.var_names).reindex(G).to_numpy()
muA = pd.Series(mu, index=c.var_names).reindex(G).to_numpy()
zs = np.load(f"{R}/submission_2026/data/shared_response.npz", allow_pickle=True)
# shared without H1 (avoid leak)
src = {}
import importlib.util
spec = importlib.util.spec_from_file_location("srm", f"{R}/submission_2026/shared_response.py"); srm = importlib.util.module_from_spec(spec); spec.loader.exec_module(srm)
Dd = f"{PF}/public-screens/data"
for k, f, t in (("K562", "K562_gwps_normalized_bulk_01.h5ad", .2), ("RPE1", "rpe1_normalized_bulk_01.h5ad", .2), ("HepG2", "hepg2_bulk.h5ad", None), ("Jurkat", "jurkat_bulk.h5ad", None)):
    src[k] = srm.zmean(f"{Dd}/{f}", t)
shared = pd.concat({k: v / v.abs().median() for k, v in src.items()}, axis=1).reindex(G).mean(axis=1).fillna(0).to_numpy()
m = np.expm1(cm)
def lognorm(e):  # relative change e -> H1 log1p(CP10k) shift
    return np.log1p(m * (1 + e)) - np.log1p(m)
def wide_mask(kappa, cap, n, score):
    e = kappa * cvA / 20
    ok = np.isfinite(e) & (muA >= 5) & (e <= cap)
    return ok, np.where(ok, e, 0)
panel = np.isin(G, T)
def pds(Pm, Rm):
    keep = ~panel
    A = Pm[:, keep]; B = Rm[:, keep]
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12); Bn = B / np.linalg.norm(B, axis=1, keepdims=True)
    S = An @ Bn.T
    gt = (S > np.diag(S)[:, None]).sum(1); eq = (S == np.diag(S)[:, None]).sum(1) - 1
    return (1 - (gt + eq / 2) / (len(S) - 1)).mean()
real = pd.read_csv(f"{PF}/logo-shap/v2/h1/celleval/mean_effect/real_de.csv", usecols=["target", "feature", "log2_fold_change", "fdr"])
real = real[(real.target != real.feature) & real.target.isin(T)]
gi = {g: i for i, g in enumerate(G)}; ti = {t: i for i, t in enumerate(T)}
def de_scores(Pm, called):
    """fid with the called set (bool targets x genes), nmae over real-sig genes (log2 of mean ratio)."""
    fids, nm = [], []
    for t, g in real.groupby("target"):
        i = ti[t]; cols = np.array([gi.get(f, -1) for f in g.feature]); okc = cols >= 0
        g = g[okc]; cols = cols[okc]
        lr = g.log2_fold_change.to_numpy(); sig = (g.fdr < .05).to_numpy(); adj = np.isfinite(lr) & (lr != 0)
        dp = Pm[i, cols]; cl = called[i, cols] & adj
        k = (cl & (np.sign(dp) == np.sign(lr))).sum()
        fids.append(k / max(cl.sum(), sig.sum(), 1))
        lp = np.log2((np.expm1(cm[cols] + dp) + 1e-9) / (m[cols] + 1e-9))
        gate = sig & np.isfinite(lr) & np.isfinite(lp)
        nm.append(np.abs(lp[gate] - lr[gate]).mean() / np.abs(lr[gate]).mean())
    return np.mean(fids), np.mean(nm)
dev = pred - pred.mean(0)
res = []
base_called = np.abs(pred) > lognorm(np.full_like(cm, 0.15))  # c1 calls ~only big moves
res.append(("c1 (submitted 1st)", pds(pred, Rl), *de_scores(pred, base_called)))
for kappa in (5, 4, 3.5):
    ok, e = wide_mask(kappa, .5, 8000, shared)
    for name, sgn in (("shared signs (submitted 2nd)" , np.sign(shared)[None, :] * np.ones_like(pred)),
                      ("target-specific deviation signs", np.sign(dev)),
                      ("deviation + 0.3 shared", np.sign(dev / np.abs(dev).mean(1, keepdims=True) + 0.3 * shared / np.abs(shared).mean()))):
        sc = np.where(ok, np.abs(shared), -1)
        pick = np.zeros_like(ok); pick[np.argsort(-sc)[:8000]] = True; pick &= ok
        w = np.where(pick[None, :], sgn * lognorm(e)[None, :], 0)
        Pm = np.where(np.abs(pred) >= np.abs(w), pred, w)
        called = base_called | (pick[None, :] & (np.abs(Pm) >= np.abs(w) - 1e-12) & (w != 0))
        res.append((f"k{kappa} {name}", pds(Pm, Rl), *de_scores(Pm, called)))
    print(pd.DataFrame(res, columns=["variant", "pds_cos", "fid_raw", "nmae"]).round(3).to_string(), flush=True)
