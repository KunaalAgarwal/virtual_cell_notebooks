"""Build a GRNBoost2 input from the Replogle K562 genome-wide screen instead of control cells.

Rows are the ~11k perturbation pseudobulk profiles (z-scored vs controls). Columns are the
2,000 H1 scoring genes that vary most across K562 perturbations, plus every H1 knockdown
target that K562 reads out. A tree model fit across perturbations learns
"when gene j moves because something was knocked down, gene g moves too", i.e. it sees
interventions rather than only natural co-variation in unperturbed cells.

    python prep_k562_rows.py h1_pseudobulk.npz grn_h1_pc0.top_edges.csv K562_gwps_normalized_bulk_01.h5ad out.npz [n_var] [n_extra_targets]
"""
import sys
import anndata as ad
import numpy as np


def main(pb_npz, genes_csv, screen, dst, n_var=2000, n_extra_targets=600):
    import pandas as pd
    genes = list(dict.fromkeys(pd.read_csv(genes_csv, usecols=["target"])["target"]))
    perts = set(np.load(pb_npz, allow_pickle=True)["perts"])
    a = ad.read_h5ad(screen)
    X = np.asarray(a.X, dtype=np.float32)
    X[~np.isfinite(X)] = 0
    np.clip(X, -10, 10, out=X)
    sg = list(a.var["gene_name"])
    si = {}
    for j, g in enumerate(sg):
        si.setdefault(g, j)
    cand = [g for g in genes if g in si]
    v = X[:, [si[g] for g in cand]].var(0)
    top = {cand[j] for j in np.argsort(-v)[:n_var]}
    use = [g for g in cand if g in top or g in perts]
    X = X[:, [si[g] for g in use]]
    print(f"{X.shape[0]} K562 perturbations x {len(use)} genes ({sum(g in perts for g in use)} H1 knockdown targets)")
    # Fit only the H1 knockdown targets plus a random sample of other genes as targets
    # (run_grn.py reads `targets`); every gene stays a candidate regulator.
    kd = [g for g in use if g in perts]
    rest = [g for g in use if g not in perts]
    extra = list(np.random.default_rng(0).choice(rest, min(n_extra_targets, len(rest)), replace=False))
    np.savez_compressed(dst, X=X, genes=np.array(use), targets=np.array(kd + extra))


if __name__ == "__main__":
    main(*sys.argv[1:5], *(int(x) for x in sys.argv[5:]))
