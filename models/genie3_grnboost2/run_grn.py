"""Infer a GRNBoost2 or GENIE3 network with arboreto's own per-target regressors.

arboreto.core.infer_partial_network is the exact routine arboreto's grnboost2() / genie3()
call for each target gene (same hyperparameters, early stopping for GRNBoost2); we call it
through joblib instead of dask, which arboreto 0.1.6 no longer runs on with current dask.
Every gene is a candidate regulator of every other gene (no TF list), so knocked-down
genes that are not TFs still get out-edges.

    python run_grn.py controls.npz grnboost2|genie3 out.npz [n_cells] [n_trees] [n_jobs]
If the input npz has a `targets` array, only those genes are fit as targets (other columns
of imp stay 0). Output: imp (regulators x targets importance, regulator == target set to 0), genes, rho
(Spearman sign between regulator and target, as pySCENIC uses to sign edges).
"""
import sys
import time
import numpy as np
from joblib import Parallel, delayed
from arboreto import core
from scipy.stats import spearmanr


def one(kind, kw, X, genes, t, seed):
    reg = np.delete(np.arange(X.shape[1]), t)
    df = core.infer_partial_network(kind, kw, X[:, reg], list(genes[reg]), genes[t], X[:, t],
                                    early_stop_window_length=core.EARLY_STOP_WINDOW_LENGTH, seed=seed)
    col = np.zeros(X.shape[1], np.float32)
    gi = {g: i for i, g in enumerate(genes)}
    col[[gi[g] for g in df["TF"]]] = df["importance"].to_numpy()
    return col


def main(src, method, dst, n_cells=0, n_trees=0, n_jobs=4, max_targets=0):
    d = np.load(src, allow_pickle=True)
    X, genes = d["X"].astype(np.float32), d["genes"]
    if n_cells and n_cells < len(X):
        X = X[np.random.default_rng(1).choice(len(X), n_cells, replace=False)]
    if method == "grnboost2":
        kind, kw = "GBM", dict(core.SGBM_KWARGS)
    else:
        kind, kw = "RF", dict(core.RF_KWARGS)
        if n_trees:
            kw["n_estimators"] = n_trees
    targets = range(X.shape[1]) if not max_targets else range(max_targets)
    if "targets" in d.files:                       # optional subset of target genes (columns of imp)
        targets = [list(genes).index(g) for g in d["targets"]]
    t0 = time.time()
    cols = Parallel(n_jobs=n_jobs, verbose=5)(delayed(one)(kind, kw, X, genes, t, 666) for t in targets)
    imp = np.zeros((X.shape[1], X.shape[1]), np.float32)
    imp[:, list(targets)] = np.stack(cols, 1)
    print(f"{method}: {imp.shape} in {time.time() - t0:.0f}s on {len(X)} cells, kwargs {kw}")
    rho = spearmanr(X).statistic.astype(np.float32) if not max_targets else None
    np.savez_compressed(dst, imp=imp, genes=genes, rho=rho, n_cells=len(X), kwargs=str(kw),
                        targets=genes[list(targets)])


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2], *(int(x) for x in a[3:]))
