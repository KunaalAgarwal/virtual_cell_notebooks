"""Cross-line shared knockdown response, one signed score per gene of the 2026 gene axis.

Most CRISPRi knockdowns move a common set of genes in the same direction (a generic response on
top of each target's specific effect; Replogle et al. 2022, Cell; Nadig et al. 2025). The 2026
scorer's `de_direction_fidelity_yield_raw` counts correct-direction DE calls divided by
max(genes called, genes the real knockdown changes), so a prediction that calls few genes is
charged for every real DE gene it misses. Calling the shared-response genes with their
cross-line sign fixes that. This script averages the mean knockdown response of five lines,
each scaled by its median |value| so no line dominates:

  * K562 genome-wide and RPE1 essential screens (Replogle et al. 2022), z-scored pseudobulk,
    restricted to their strongest 20% of knockdowns (the 2026 panel was chosen for effects),
  * HepG2 and Jurkat (Nadig et al. 2025, 2025 VCC support set), z-scored pseudobulk,
  * H1 (2025 VCC training data), log1p(CP10k) differences.

Checked leave-one-line-out (README, "Shared response"): the average of the other lines gets the
direction of a held-out line's individual knockdowns right on 53-56% of its top 4,000-8,000
genes, about as well as the best single source and never below chance.

    python shared_response.py --data vcc2026-data/ --h1 data/h1_pseudobulk.npz \
        --k562 K562_gwps_normalized_bulk_01.h5ad --rpe1 rpe1_normalized_bulk_01.h5ad \
        --hepg2 hepg2_bulk.h5ad --jurkat jurkat_bulk.h5ad --out data/shared_response.npz
"""

from __future__ import annotations

import argparse

import anndata as ad
import numpy as np
import pandas as pd


def zmean(path, top_frac=None):
    """Mean z-scored knockdown profile of a Replogle-layout pseudobulk file, per gene symbol."""
    a = ad.read_h5ad(path)
    X = np.nan_to_num(np.asarray(a.X, dtype=np.float32), posinf=0, neginf=0).clip(-10, 10)
    tg = pd.Index(a.obs_names.str.split("_").str[1])
    X = X[~tg.str.contains("non-targeting")]
    if top_frac:
        n = (np.abs(X) > 1).sum(1)
        X = X[n >= np.quantile(n, 1 - top_frac)]
    return pd.Series(X.mean(0), index=a.var["gene_name"].astype(str).values).groupby(level=0).mean()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("data", "h1", "k562", "rpe1", "hepg2", "jurkat", "out"):
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args()
    genes = pd.read_csv(f"{a.data}/gene_names.csv")["gene_name"].astype(str).tolist()
    z = np.load(a.h1, allow_pickle=True)
    src = {"K562": zmean(a.k562, 0.2), "RPE1": zmean(a.rpe1, 0.2),
           "HepG2": zmean(a.hepg2), "Jurkat": zmean(a.jurkat),
           "H1": pd.Series(z["delta"].mean(0), index=z["genes"].astype(str)).groupby(level=0).mean()}
    M = pd.concat({k: v / v.abs().median() for k, v in src.items()}, axis=1).reindex(genes)
    score = M.mean(axis=1, skipna=True).fillna(0.0).to_numpy(np.float32)
    n_src = M.notna().sum(axis=1).to_numpy(np.int8)
    np.savez_compressed(a.out, genes=np.array(genes), score=score, n_sources=n_src,
                        sources=np.array(list(src)))
    print(f"{len(genes)} genes, {int((n_src > 0).sum())} with a score; sources per gene: "
          f"{pd.Series(n_src).value_counts().sort_index().to_dict()}")


if __name__ == "__main__":
    main()
