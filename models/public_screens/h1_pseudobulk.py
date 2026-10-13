"""Batch-matched pseudobulk knockdown effects for the 2025 H1 training set.

Streams the 15 GB CSR .h5ad in row chunks, normalizes each cell to log1p(CP10k), and
accumulates per-(target_gene, batch) sums. The true effect of knocking down gene k is

    delta_k = mean over k's cells of (x_cell - control mean of the cell's own batch)

so batch differences don't leak into the effect. Writes an .npz with delta (perts x
genes), control mean, gene names, perturbation names and cell counts.

    python h1_pseudobulk.py adata_Training.h5ad h1_pseudobulk.npz
"""

import sys

import h5py
import numpy as np
import scipy.sparse as sp

CHUNK = 4000


def _col(f, name):
    g = f["obs"][name]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]], dtype=object)
    return cats[g["codes"][:]]


def main(src, dst):
    with h5py.File(src, "r") as f:
        n, n_genes = (int(v) for v in f["X"].attrs["shape"])
        var = f["var"]["_index"][:]
        genes = np.array([g.decode() if isinstance(g, bytes) else g for g in var], dtype=object)
        target, batch = _col(f, "target_gene"), _col(f, "batch")
        perts, p_idx = np.unique(target, return_inverse=True)
        batches, b_idx = np.unique(batch, return_inverse=True)
        key = p_idx * len(batches) + b_idx
        sums = np.zeros((len(perts) * len(batches), n_genes), dtype=np.float64)
        counts = np.bincount(key, minlength=sums.shape[0])

        indptr = f["X/indptr"][:]
        data_ds, idx_ds = f["X/data"], f["X/indices"]
        # The public 2025 support-set competition_train.h5 stores log1p(counts); undo that.
        log1p_counts = np.modf(data_ds[:100000])[0].max() > 1e-3
        print("input is log1p(counts):", log1p_counts)
        for start in range(0, n, CHUNK):
            stop = min(start + CHUNK, n)
            lo, hi = indptr[start], indptr[stop]
            data = data_ds[lo:hi].astype(np.float64)
            if log1p_counts:
                data = np.rint(np.expm1(data))
            cols = idx_ds[lo:hi]
            row_ptr = indptr[start:stop + 1] - lo
            rows = np.repeat(np.arange(stop - start), np.diff(row_ptr))
            lib = np.bincount(rows, weights=data, minlength=stop - start)
            vals = np.log1p(data * (1e4 / lib[rows]))
            X = sp.csr_matrix((vals, cols, row_ptr), shape=(stop - start, n_genes))
            # Indicator (groups x cells) @ cells sums each cell's row into its (pert, batch) group.
            k = key[start:stop]
            ind = sp.csr_matrix((np.ones(len(k)), (k, np.arange(len(k)))), shape=(sums.shape[0], len(k)))
            g = np.unique(k)
            sums[g] += (ind[g] @ X).toarray()
            if (start // CHUNK) % 10 == 0:
                print(f"{stop}/{n}", flush=True)

    sums = sums.reshape(len(perts), len(batches), n_genes)
    counts = counts.reshape(len(perts), len(batches))
    ctrl = int(np.flatnonzero(perts == "non-targeting")[0])
    ctrl_batch_mean = sums[ctrl] / np.maximum(counts[ctrl], 1)[:, None]           # batches x genes
    ctrl_mean = sums[ctrl].sum(0) / counts[ctrl].sum()
    # Each perturbed cell is compared with its own batch's controls.
    pert_sum = sums.sum(1)
    matched_ctrl = counts @ ctrl_batch_mean                                        # perts x genes
    n_cells = counts.sum(1)
    delta = (pert_sum - matched_ctrl) / n_cells[:, None]
    keep = np.arange(len(perts)) != ctrl
    np.savez_compressed(dst, delta=delta[keep].astype(np.float32), ctrl_mean=ctrl_mean.astype(np.float32),
                        genes=genes.astype(str), perts=perts[keep].astype(str), n_cells=n_cells[keep])
    print(f"{keep.sum()} perturbations x {n_genes} genes -> {dst}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
