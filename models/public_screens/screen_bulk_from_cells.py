"""Turn a single-cell Perturb-seq .h5/.h5ad (log-normalized X, obs target_gene) into a z-scored
pseudobulk file in the same layout as Replogle's *_normalized_bulk_01.h5ad, so
screen_transfer_h1.load_screen can read it.

Profile of knockdown k = (mean of k's cells - mean of controls in the same batch) / control std,
averaged over batches weighted by k's cells. Readout genes with < min_detect controls expressing
them are dropped.

    python screen_bulk_from_cells.py hepg2.h5 hepg2_bulk.h5ad [--batch batch_var]
"""

import argparse

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--batch", default="batch_var")
    ap.add_argument("--control", default="non-targeting")
    ap.add_argument("--min-detect", type=float, default=0.05)
    a = ap.parse_args()

    x = ad.read_h5ad(a.src)
    X = x.X.toarray() if sp.issparse(x.X) else np.asarray(x.X)
    X = X.astype(np.float32)
    tg = x.obs["target_gene"].astype(str).to_numpy()
    bt = x.obs[a.batch].astype(str).to_numpy() if a.batch in x.obs else np.zeros(len(tg), str)
    ctrl = tg == a.control
    keep_g = (X[ctrl] > 0).mean(0) >= a.min_detect
    X = X[:, keep_g]; genes = x.var_names.to_numpy()[keep_g]
    sd = X[ctrl].std(0) + 1e-3
    gm = X[ctrl].mean(0)
    bmean = {b: X[ctrl & (bt == b)].mean(0) if (ctrl & (bt == b)).sum() >= 20 else gm for b in np.unique(bt)}
    rows, names = [], []
    for k in sorted(set(tg) - {a.control}):
        idx = np.flatnonzero(tg == k)
        ref = np.stack([bmean[b] for b in bt[idx]]).mean(0)
        rows.append((X[idx].mean(0) - ref) / sd)
        names.append(f"{len(names)}_{k}_bulk")
    out = ad.AnnData(X=np.stack(rows).astype(np.float32), obs=pd.DataFrame(index=names),
                     var=pd.DataFrame({"gene_name": genes}, index=genes))
    out.write_h5ad(a.dst)
    print(f"{out.n_obs} knockdowns x {out.n_vars} readout genes -> {a.dst}")


if __name__ == "__main__":
    main()
