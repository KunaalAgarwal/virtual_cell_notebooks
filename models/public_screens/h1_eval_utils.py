"""Shared H1 scoring helpers (from the GRN exploration's eval_h1.py and celleval_h1.py).

score            pseudobulk Pearson / top-100 Pearson / top-50 sign agreement, target gene excluded
build_pred       a cell-eval submission: resampled control cells shifted by a per-gene log change
control_mean_pred  the repo baseline: every cell = rounded control mean
"""

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SELF_DROP = -0.633  # H1 median change of the knocked-down gene itself, log1p(CP10k)


def score(pred, true, k):
    mask = np.ones(len(true), bool)
    mask[k] = False
    p, t = pred[mask], true[mask]
    out = {}
    out["pearson_all"] = np.corrcoef(p, t)[0, 1] if p.std() > 0 else np.nan
    top = np.argsort(-np.abs(t))
    t100 = top[:100]
    out["pearson_top100"] = np.corrcoef(p[t100], t[t100])[0, 1] if p[t100].std() > 0 else np.nan
    t50 = top[:50]
    out["sign_top50"] = float(np.mean(np.sign(p[t50]) == np.sign(t[t50])))
    return out


def to_lognorm(X):
    X = sp.csr_matrix(X, dtype=np.float32)
    lib = np.asarray(X.sum(1)).ravel()
    Xn = sp.diags((1e4 / np.maximum(lib, 1)).astype(np.float32)) @ X
    Xn.data = np.log1p(Xn.data)
    return Xn, lib


def shifted_counts(ctrl_log, lib, delta):
    """Apply a per-gene log change to control cells and return integer counts."""
    out = np.asarray(ctrl_log.todense(), dtype=np.float32)
    out += delta[None, :].astype(np.float32)
    np.clip(out, 0, None, out=out)
    counts = np.expm1(out) * (lib / 1e4)[:, None]
    return sp.csr_matrix(np.rint(counts).astype(np.float32))


def build_pred(real, deltas, cells_per_pert, rng):
    """deltas: dict pert -> change vector over real.var_names (zeros for 'no change')."""
    ctrl_mask = (real.obs["target"] == "non-targeting").to_numpy()
    ctrl = real[ctrl_mask]
    ctrl_log, lib = to_lognorm(ctrl.X)
    blocks, obs = [ctrl.X.copy()], [pd.DataFrame({"target": ["non-targeting"] * ctrl.n_obs})]
    for pert, d in deltas.items():
        pick = rng.choice(ctrl.n_obs, cells_per_pert, replace=cells_per_pert > ctrl.n_obs)
        blocks.append(shifted_counts(ctrl_log[pick], lib[pick], d))
        obs.append(pd.DataFrame({"target": [pert] * cells_per_pert}))
    obs = pd.concat(obs, ignore_index=True)
    obs.index = [f"p{i}" for i in range(len(obs))]
    return ad.AnnData(X=sp.vstack(blocks).tocsr(), obs=obs, var=real.var.copy())


def control_mean_pred(real, perts, cells_per_pert):
    ctrl = real[(real.obs["target"] == "non-targeting").to_numpy()]
    mean = np.rint(np.asarray(ctrl.X.mean(0)).ravel()).astype(np.float32)
    row = sp.csr_matrix(mean)
    blocks = [ctrl.X.copy()] + [sp.vstack([row] * cells_per_pert) for _ in perts]
    obs = pd.DataFrame({"target": ["non-targeting"] * ctrl.n_obs + [p for p in perts for _ in range(cells_per_pert)]})
    obs.index = [f"p{i}" for i in range(len(obs))]
    return ad.AnnData(X=sp.vstack(blocks).tocsr(), obs=obs, var=real.var.copy())
