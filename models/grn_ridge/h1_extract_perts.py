"""Pull a test set of real H1 knockdown cells (raw counts) for scoring with cell-eval.

Picks `--n-perts` knockdowns stratified by effect size (strong / medium / weak, from the
pseudobulk file), keeps up to `--cells-per-pert` cells for each, plus `--n-controls`
non-targeting cells. Writes an .h5ad with obs column `target` (cell-eval's default).

    python h1_extract_perts.py adata_Training.h5ad h1_pseudobulk.npz h1_real_test.h5ad
"""

import argparse

import anndata as ad
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("src")
    p.add_argument("pseudobulk")
    p.add_argument("dst")
    p.add_argument("--n-perts", type=int, default=40)
    p.add_argument("--cells-per-pert", type=int, default=400)
    p.add_argument("--n-controls", type=int, default=2000)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()

    pb = np.load(a.pseudobulk, allow_pickle=True)
    perts, delta = pb["perts"], pb["delta"]
    eff = np.abs(delta).mean(1)
    order = np.argsort(-eff)
    n = a.n_perts
    pick = np.concatenate([order[: n * 3 // 8], order[len(order) // 2 - n * 3 // 16: len(order) // 2 + n * 3 // 16],
                           order[-(n - n * 3 // 8 - 2 * (n * 3 // 16)):]])
    chosen = sorted(set(perts[pick]))
    print(f"{len(chosen)} knockdowns:", ", ".join(chosen))

    rng = np.random.default_rng(a.seed)
    full = ad.read_h5ad(a.src, backed="r")
    tg = full.obs["target_gene"].astype(str).to_numpy()
    keep = []
    for g in chosen + ["non-targeting"]:
        idx = np.flatnonzero(tg == g)
        cap = a.n_controls if g == "non-targeting" else a.cells_per_pert
        keep.extend(rng.choice(idx, min(cap, len(idx)), replace=False))
    keep = np.sort(np.array(keep))
    sub = full[keep].to_memory()
    sub.obs["target"] = sub.obs["target_gene"].astype(str)
    sub.write_h5ad(a.dst, compression="gzip")
    print(f"{sub.n_obs} cells x {sub.n_vars} genes -> {a.dst}")


if __name__ == "__main__":
    main()
