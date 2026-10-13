"""Leave-one-group-out CV + SHAP feature importance on non-targeting control cells.

Implements the edge-learning half of grn_linear_model.pdf: model one gene as a function
of all the others, N_i = f(N_-i), on control cells from a single cell context.

  1. Partition the cells into groups (default: their non-targeting guide, so 46 groups
     of 400 for a 2026 context, 31 unequal groups for H1 2025).
  2. Leave-one-group-out CV: train on all groups but one, predict the held-out group.
  3. SHAP the held-out cells (one row of SHAP values per cell).
  4. Collapse each fold's per-cell SHAP matrix to one row per fold, write that CSV.
  5. Average the fold rows into a single row per target gene. Stacking those rows over
     target genes gives the edge matrix (targets x regulators).

Run on its own:
    python logo_shap.py --h5ad context_A.h5ad --group-col ntc_id --targets DPM1 NFYA \
        --models ridge lgbm --out results/
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

AGG_STATS = ("mean_abs_shap", "signed_shap", "native_importance")


# ---------------------------------------------------------------- data

def load_controls(h5ad_path, group_col, max_cells_per_group=None, seed=0):
    """Read a controls .h5ad, optionally capping cells per group (for quick local runs)."""
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path)
    if "target_gene" in adata.obs and (adata.obs["target_gene"] != "non-targeting").any():
        adata = adata[adata.obs["target_gene"] == "non-targeting"].copy()
    if max_cells_per_group:
        rng = np.random.default_rng(seed)
        keep = []
        for _, idx in adata.obs.groupby(group_col, observed=True).indices.items():
            keep.extend(rng.choice(idx, min(len(idx), max_cells_per_group), replace=False))
        adata = adata[np.sort(keep)].copy()
    return adata


def normalize(X, target_sum=1e4):
    """Library-size normalize to target_sum counts per cell, then log1p."""
    X = sp.csr_matrix(X, dtype=np.float32)
    libsize = np.asarray(X.sum(axis=1)).ravel()
    libsize[libsize == 0] = 1
    X = sp.diags((target_sum / libsize).astype(np.float32)) @ X
    X.data = np.log1p(X.data)
    return X


def select_features(X, gene_names, n_top=2000, min_frac=0.05, always=()):
    """Highly variable genes among those detected in >= min_frac of cells.

    n_top=None keeps every gene that passes min_frac (the full GENIE3-style run).
    Genes in `always` (e.g. the targets) are kept regardless.
    """
    detected = np.asarray((X > 0).mean(axis=0)).ravel()
    mean = np.asarray(X.mean(axis=0)).ravel()
    var = np.asarray(X.multiply(X).mean(axis=0)).ravel() - mean**2
    ok = detected >= min_frac
    order = np.argsort(-np.where(ok, var, -np.inf))
    n_ok = int(ok.sum())
    chosen = order[: n_ok if n_top is None else min(n_top, n_ok)]
    names = np.asarray(gene_names)
    extra = [np.flatnonzero(names == g)[0] for g in always if g in set(names)]
    chosen = np.unique(np.concatenate([chosen, extra]).astype(int))
    return chosen


# ---------------------------------------------------------------- models

def make_model(name, seed=0):
    if name == "ridge":
        from sklearn.linear_model import RidgeCV
        return RidgeCV(alphas=np.logspace(-1, 4, 11))
    if name == "lasso":
        from sklearn.linear_model import LassoCV
        return LassoCV(n_alphas=20, cv=3, max_iter=2000, random_state=seed, n_jobs=-1)
    if name == "lgbm":
        from lightgbm import LGBMRegressor
        return LGBMRegressor(
            n_estimators=200, learning_rate=0.05, num_leaves=31, subsample=0.8,
            subsample_freq=1, colsample_bytree=0.3, min_child_samples=20,
            random_state=seed, n_jobs=-1, verbose=-1,
        )
    raise ValueError(f"unknown model {name!r}")


def native_importance(model):
    """The model's own importance: |coef| for linear models, total gain for trees."""
    if hasattr(model, "coef_"):
        return np.abs(np.ravel(model.coef_))
    if hasattr(model, "booster_"):
        return model.booster_.feature_importance(importance_type="gain").astype(float)
    return np.asarray(model.feature_importances_, dtype=float)


def shap_values(model, X_train, X_test, background=200, seed=0):
    import shap

    if hasattr(model, "coef_"):
        rng = np.random.default_rng(seed)
        bg = X_train[rng.choice(len(X_train), min(background, len(X_train)), replace=False)]
        explainer = shap.LinearExplainer(model, bg)
    else:
        explainer = shap.TreeExplainer(model)
    return np.asarray(explainer.shap_values(X_test))


def fold_row(shap_mat, X_test, native):
    """Collapse one fold's (cells x features) SHAP matrix into per-feature summaries.

    mean_abs_shap      mean over held-out cells of |SHAP|: how much the gene moves the
                       prediction, the edge weight.
    signed_shap        mean_abs_shap * sign(corr(feature value, SHAP value)): + means
                       higher expression of the regulator pushes the target up
                       (activation), - means down (repression). Plain mean SHAP is ~0
                       by construction, so it cannot carry the sign.
    native_importance  the model's own importance (|coef| or gain), for comparison.
    """
    mean_abs = np.abs(shap_mat).mean(axis=0)
    xc = X_test - X_test.mean(axis=0)
    sc = shap_mat - shap_mat.mean(axis=0)
    cov = (xc * sc).sum(axis=0)
    sign = np.sign(cov)
    return {"mean_abs_shap": mean_abs, "signed_shap": mean_abs * sign, "native_importance": native}


# ---------------------------------------------------------------- LOGO CV

def logo_for_target(Xn, gene_names, feat_idx, target, groups, model_name, seed=0, log=print):
    """Run LOGO CV + SHAP for one target gene with one model.

    Returns (per_fold DataFrame: one row per (fold, stat), metrics DataFrame).
    """
    from sklearn.model_selection import LeaveOneGroupOut
    from scipy.stats import pearsonr

    names = np.asarray(gene_names)
    t_idx = int(np.flatnonzero(names == target)[0])
    f_idx = feat_idx[feat_idx != t_idx]  # never use the target to predict itself
    X = Xn[:, f_idx].toarray()
    y = Xn[:, t_idx].toarray().ravel()
    feats = names[f_idx]

    rows, metrics = [], []
    for k, (tr, te) in enumerate(LeaveOneGroupOut().split(X, y, groups)):
        t0 = time.perf_counter()
        model = make_model(model_name, seed)
        model.fit(X[tr], y[tr])
        pred = model.predict(X[te])
        shap_mat = shap_values(model, X[tr], X[te], seed=seed)
        summary = fold_row(shap_mat, X[te], native_importance(model))
        held = str(groups[te][0])
        for stat, vec in summary.items():
            rows.append({"target": target, "model": model_name, "fold": k,
                         "held_out_group": held, "stat": stat, **dict(zip(feats, vec))})
        ss_res = ((y[te] - pred) ** 2).sum()
        ss_tot = ((y[te] - y[tr].mean()) ** 2).sum()
        r = pearsonr(y[te], pred)[0] if np.std(pred) > 0 and np.std(y[te]) > 0 else np.nan
        metrics.append({"target": target, "model": model_name, "fold": k, "held_out_group": held,
                        "n_test": len(te), "r2": 1 - ss_res / ss_tot, "pearson": r,
                        "seconds": time.perf_counter() - t0})
    m = pd.DataFrame(metrics)
    log(f"  {target:>12s} {model_name:6s} folds={len(m)}  R2={m.r2.mean():.3f}  "
        f"pearson={m.pearson.mean():.3f}  ({m.seconds.sum():.0f}s)")
    return pd.DataFrame(rows), m


def average_folds(per_fold, weights=None):
    """Average fold rows into one row per (target, model, stat).

    Unweighted mean over folds by default; pass weights={group: n_cells} to weight by
    held-out group size (matters for H1, whose guides have unequal cell counts).
    """
    meta = ["target", "model", "fold", "held_out_group", "stat"]
    feats = [c for c in per_fold.columns if c not in meta]
    out = []
    for (t, m, s), df in per_fold.groupby(["target", "model", "stat"], sort=False):
        w = None if weights is None else df["held_out_group"].map(weights).to_numpy(float)
        vals = np.average(df[feats].to_numpy(float), axis=0, weights=w)
        out.append({"target": t, "model": m, "stat": s, "n_folds": len(df), **dict(zip(feats, vals))})
    return pd.DataFrame(out)


def write_outputs(out_dir, per_fold, metrics, models, weights=None):
    per_fold = pd.concat(per_fold, ignore_index=True)
    single = average_folds(per_fold, weights)
    per_fold.to_csv(out_dir / "shap_per_fold.csv", index=False)
    pd.concat(metrics, ignore_index=True).to_csv(out_dir / "cv_metrics.csv", index=False)
    single.to_csv(out_dir / "shap_single_row.csv", index=False)
    # Edge matrix view: rows are targets, columns are regulators.
    for m in models:
        for stat in ("mean_abs_shap", "signed_shap"):
            em = single[(single.model == m) & (single.stat == stat)].drop(
                columns=["model", "stat", "n_folds"]).set_index("target")
            em.to_csv(out_dir / f"edge_matrix_{m}_{stat}.csv")
    return single


def run(h5ad, group_col, targets, models=("ridge", "lgbm"), out_dir="results",
        n_top_features=2000, min_frac=0.05, max_cells_per_group=None, weight_by_size=False,
        seed=0, log=print):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata = load_controls(h5ad, group_col, max_cells_per_group, seed)
    genes = adata.var_names.to_numpy()
    groups = adata.obs[group_col].astype(str).to_numpy()
    log(f"{adata.n_obs} cells x {adata.n_vars} genes, {len(np.unique(groups))} groups by {group_col!r}")

    Xn = normalize(adata.X)
    missing = [t for t in targets if t not in set(genes)]
    if missing:
        log(f"skipping targets not in the data: {missing}")
    targets = [t for t in targets if t in set(genes)]
    feat_idx = select_features(Xn, genes, n_top_features, min_frac, always=targets)
    Xn = Xn.tocsc()
    log(f"{len(feat_idx)} candidate regulator genes; {len(targets)} targets; models {list(models)}")

    weights = pd.Series(groups).value_counts().to_dict() if weight_by_size else None
    per_fold, metrics = [], []
    for t in targets:
        for m in models:
            pf, mt = logo_for_target(Xn, genes, feat_idx, t, groups, m, seed, log)
            per_fold.append(pf)
            metrics.append(mt)
        # Rewrite the outputs after every target so a long run that stops early keeps its work.
        single = write_outputs(out_dir, per_fold, metrics, models, weights)
    log(f"wrote {out_dir}/shap_per_fold.csv, shap_single_row.csv, cv_metrics.csv, edge_matrix_*.csv")
    per_fold = pd.concat(per_fold, ignore_index=True)
    metrics = pd.concat(metrics, ignore_index=True)
    return per_fold, single, metrics


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True, help="controls .h5ad (one cell context)")
    p.add_argument("--group-col", default="ntc_id", help="obs column defining LOGO groups (ntc_id for 2026, guide_id for H1)")
    p.add_argument("--targets", nargs="+", required=True, help="genes to model as N_i = f(N_-i)")
    p.add_argument("--models", nargs="+", default=["ridge", "lgbm"], choices=["ridge", "lasso", "lgbm"])
    p.add_argument("--out", default="results")
    p.add_argument("--n-top-features", type=int, default=2000, help="HVG regulators to use; 0 = all detected genes")
    p.add_argument("--min-frac", type=float, default=0.05, help="min fraction of cells a regulator is detected in")
    p.add_argument("--max-cells-per-group", type=int, default=None, help="subsample for quick runs")
    p.add_argument("--weight-by-size", action="store_true", help="weight the fold average by held-out group size")
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    run(a.h5ad, a.group_col, a.targets, a.models, a.out, a.n_top_features or None, a.min_frac,
        a.max_cells_per_group, a.weight_by_size, a.seed)


if __name__ == "__main__":
    main()
