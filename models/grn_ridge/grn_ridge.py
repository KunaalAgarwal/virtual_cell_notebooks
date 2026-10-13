"""All-genes ridge GRN from control cells, with grouped CV and stability stats.

The fast path for stage 1 of grn_linear_model.pdf. Instead of fitting one regression per
target gene, every "gene ~ all other genes" ridge comes out of one matrix inverse:
with standardized expression Z and P = (Z'Z/n + lam*I)^-1, the coefficients predicting
gene i are  -P[:, i] / P[i, i]  (exactly what sklearn's Ridge gives per gene). One
eigendecomposition per fold makes every lam essentially free.

Per fold (GroupKFold over guides, or full leave-one-group-out):
  coef      regulator j -> target i, in log1p(CP10k) units (the slope, for propagation)
  shap      mean |SHAP| on held-out cells = |coef| * mean|x_j - mean_j| (exact for linear)
  r2        held-out R^2 for every target
Across folds:
  coef_mean, coef_std, sign_agree (fraction of folds whose sign matches the mean's),
  shap_share (each fold's |SHAP| column scaled to sum 1 per target, then averaged), r2.

Cell-state removal (--n-pcs, --covariates) regresses log library size, mitochondrial
fraction and the top principal components out of every gene before fitting, so edges
reflect residual co-variation rather than shared programs (cell cycle, quality,
differentiation).

    python grn_ridge.py --h5ad context_A.h5ad --group-col ntc_id --out grn_A.npz
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

LAMBDAS = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)


def load(h5ad, group_col, min_frac=0.05):
    """log1p(CP10k) expression (dense, genes detected in >= min_frac of cells), plus covariates."""
    import anndata as ad

    a = ad.read_h5ad(h5ad)
    if "target_gene" in a.obs and (a.obs["target_gene"] != "non-targeting").any():
        a = a[(a.obs["target_gene"] == "non-targeting").to_numpy()].copy()
    X = sp.csr_matrix(a.X, dtype=np.float32)
    lib = np.asarray(X.sum(1)).ravel()
    mt = np.asarray(X[:, np.flatnonzero(np.char.startswith(a.var_names.to_numpy().astype(str), "MT-"))].sum(1)).ravel() / lib
    X = sp.diags((1e4 / lib).astype(np.float32)) @ X
    X.data = np.log1p(X.data)
    X.eliminate_zeros()
    keep = np.bincount(X.indices, minlength=X.shape[1]) / X.shape[0] >= min_frac
    covars = np.column_stack([np.log(lib), mt])
    return X[:, keep].toarray(), a.var_names.to_numpy()[keep], a.obs[group_col].astype(str).to_numpy(), covars


def remove_cell_state(X, covars, n_pcs=10, log=print):
    """Residualize every gene on covariates and the top n_pcs principal components."""
    parts = [np.ones((len(X), 1))]
    if covars is not None:
        parts.append(covars)
    if n_pcs:
        Xs = (X - X.mean(0)) / (X.std(0) + 1e-8)
        from sklearn.utils.extmath import randomized_svd
        u, s, _ = randomized_svd(Xs, n_components=n_pcs, random_state=0)
        parts.append(u * s)
        log(f"  removed {n_pcs} PCs (variance share {(s**2).sum() / (Xs**2).sum():.1%})")
    C = np.column_stack(parts)
    beta, *_ = np.linalg.lstsq(C, X, rcond=None)
    return (X - C @ beta).astype(np.float32)


def _edges_for_lambdas(G_eig, lam):
    L, V = G_eig
    P = (V / (L + lam)) @ V.T
    d = np.diag(P).copy()
    B = -P / d[None, :]
    np.fill_diagonal(B, 0.0)
    return B  # B[j, i]: standardized coef of regulator j predicting target i


def _r2(X_te, mu, sd, B_raw):
    pred = (X_te - mu) @ B_raw + mu
    ss_res = ((X_te - pred) ** 2).sum(0)
    ss_tot = ((X_te - mu) ** 2).sum(0)
    return 1 - ss_res / np.maximum(ss_tot, 1e-12)


def fit(X, groups, n_folds=5, lambdas=LAMBDAS, lam=None, log=print):
    from sklearn.model_selection import GroupKFold, LeaveOneGroupOut

    n_genes = X.shape[1]
    splitter = LeaveOneGroupOut() if n_folds in (None, 0, "logo") else GroupKFold(n_splits=n_folds)
    coef_sum = np.zeros((n_genes, n_genes), dtype=np.float32)
    coef_sq = np.zeros((n_genes, n_genes), dtype=np.float32)
    pos = np.zeros((n_genes, n_genes), dtype=np.int32)
    share_sum = np.zeros((n_genes, n_genes), dtype=np.float32)
    shap_sum = np.zeros((n_genes, n_genes), dtype=np.float32)
    r2s, k = [], 0
    for k, (tr, te) in enumerate(splitter.split(X, groups=groups), start=1):
        t0 = time.perf_counter()
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-8
        Z = (X[tr] - mu) / sd
        G = (Z.T @ Z / len(tr)).astype(np.float64)
        del Z
        L, V = np.linalg.eigh(G)
        del G
        eig = (L.astype(np.float32), V.astype(np.float32))
        del L, V
        if lam is None:  # pick lam once, on the first fold, by mean held-out R^2
            scores = {}
            for l in lambdas:
                B = _edges_for_lambdas(eig, l) * (sd[None, :] / sd[:, None])
                scores[l] = float(np.mean(_r2(X[te], mu, sd, B)))
            lam = max(scores, key=scores.get)
            log("  lambda search (mean held-out R2): " + ", ".join(f"{l:g}={s:.3f}" for l, s in scores.items()))
        B = _edges_for_lambdas(eig, lam) * (sd[None, :] / sd[:, None])  # raw-unit slopes
        r2 = _r2(X[te], mu, sd, B)
        shap = np.abs(B) * np.abs(X[te] - mu).mean(0)[:, None]
        coef_sum += B
        coef_sq += B**2
        pos += B > 0
        shap_sum += shap
        share_sum += shap / np.maximum(shap.sum(0, keepdims=True), 1e-12)
        r2s.append(r2)
        log(f"  fold {k}: {len(te)} held-out cells, mean R2 {r2.mean():.3f}, median {np.median(r2):.3f} "
            f"({time.perf_counter() - t0:.0f}s)")
    coef_mean = coef_sum / k
    sign_agree = np.where(coef_mean > 0, pos / k, 1 - pos / k)
    np.fill_diagonal(sign_agree, 0)
    return {
        "coef_mean": coef_mean.astype(np.float32),
        "coef_std": np.sqrt(np.maximum(coef_sq / k - coef_mean**2, 0)).astype(np.float32),
        "sign_agree": sign_agree.astype(np.float32),
        "shap_mean": (shap_sum / k).astype(np.float32),
        "shap_share": (share_sum / k).astype(np.float32),
        "r2": np.mean(r2s, 0).astype(np.float32),
        "lam": lam,
        "n_folds": k,
    }


def top_edges(res, genes, per_target=20, min_sign_agree=0.8):
    """Long-format table of each target's strongest stable regulators."""
    rows = []
    share, coef, agree = res["shap_share"], res["coef_mean"], res["sign_agree"]
    for i, t in enumerate(genes):
        order = np.argsort(-share[:, i])
        order = [j for j in order if agree[j, i] >= min_sign_agree][:per_target]
        for rank, j in enumerate(order, 1):
            rows.append({"target": t, "regulator": genes[j], "rank": rank, "shap_share": share[j, i],
                         "coef": coef[j, i], "sign_agree": agree[j, i], "target_r2": res["r2"][i]})
    return pd.DataFrame(rows)


def run(h5ad, group_col, out, n_folds=5, n_pcs=10, covariates=True, min_frac=0.05, lam=None, log=print):
    t0 = time.perf_counter()
    X, genes, groups, covars = load(h5ad, group_col, min_frac)
    log(f"{X.shape[0]} cells x {X.shape[1]} genes, {len(np.unique(groups))} groups by {group_col!r}")
    if n_pcs or covariates:
        X = remove_cell_state(X, covars if covariates else None, n_pcs, log)
    res = fit(X, groups, n_folds, lam=lam, log=log)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, genes=genes.astype(str), **{k: v for k, v in res.items()})
    top_edges(res, genes).to_csv(out.with_suffix(".top_edges.csv"), index=False)
    log(f"lambda {res['lam']:g}, {res['n_folds']} folds, mean R2 {res['r2'].mean():.3f}; "
        f"wrote {out} and {out.with_suffix('.top_edges.csv').name} ({time.perf_counter() - t0:.0f}s)")
    return genes, res


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--group-col", default="ntc_id")
    p.add_argument("--out", required=True, help="output .npz")
    p.add_argument("--folds", default="5", help="number of grouped folds, or 'logo' for leave-one-group-out")
    p.add_argument("--n-pcs", type=int, default=10, help="principal components to regress out (0 = none)")
    p.add_argument("--no-covariates", action="store_true", help="don't regress out library size / MT fraction")
    p.add_argument("--min-frac", type=float, default=0.05)
    p.add_argument("--lam", type=float, default=None, help="fix lambda instead of searching")
    a = p.parse_args()
    folds = 0 if a.folds == "logo" else int(a.folds)
    run(a.h5ad, a.group_col, a.out, folds, a.n_pcs, not a.no_covariates, a.min_frac, a.lam)


if __name__ == "__main__":
    main()
