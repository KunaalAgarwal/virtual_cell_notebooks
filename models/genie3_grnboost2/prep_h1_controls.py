"""Sample H1 control cells and pick the GRN gene set (2,000 HVGs + every H1 knockdown target).

Reads the 2025 support-set competition_train.h5 (X = log1p(counts), CSR), keeps a stratified
sample of non-targeting cells, re-normalizes to log1p(CP10k), and writes a dense npz.

    python prep_h1_controls.py h1.h5 h1_pseudobulk.npz grn_h1_pc0.top_edges.csv out.npz [n_cells] [n_hvg]
"""
import sys
import h5py
import numpy as np
import pandas as pd


def col(f, name):
    g = f["obs"][name]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]], dtype=object)
    return cats[g["codes"][:]]


def main(src, pb_path, genes_csv, dst, n_cells=6000, n_hvg=2000, seed=0):
    rng = np.random.default_rng(seed)
    pb = np.load(pb_path, allow_pickle=True)
    perts = list(pb["perts"])
    score_genes = list(dict.fromkeys(pd.read_csv(genes_csv, usecols=["target"])["target"]))
    with h5py.File(src, "r") as f:
        genes = np.array([g.decode() for g in f["var"]["_index"][:]])
        tg, batch = col(f, "target_gene"), col(f, "batch")
        ctrl = np.flatnonzero(tg == "non-targeting")
        pick = np.sort(rng.choice(ctrl, size=min(n_cells, len(ctrl)), replace=False))
        indptr = f["X/indptr"][:]
        n_genes = len(genes)
        X = np.zeros((len(pick), n_genes), np.float32)
        data_ds, idx_ds = f["X/data"], f["X/indices"]
        for r, c in enumerate(pick):
            lo, hi = indptr[c], indptr[c + 1]
            X[r, idx_ds[lo:hi]] = data_ds[lo:hi]
    counts = np.expm1(X)
    lib = counts.sum(1, keepdims=True)
    X = np.log1p(counts * (1e4 / lib)).astype(np.float32)
    gi = {g: i for i, g in enumerate(genes)}
    cand = np.array([gi[g] for g in score_genes])
    det = (X[:, cand] > 0).mean(0)
    mu, var = X[:, cand].mean(0), X[:, cand].var(0)
    disp = np.where(det >= 0.05, var / (mu + 1e-6), -np.inf)
    hvg = cand[np.argsort(-disp)[:n_hvg]]
    kd = [gi[p] for p in perts if p in set(score_genes)]
    sel = np.array(sorted(set(hvg) | set(kd)))
    print(f"{len(pick)} control cells from {len(set(batch[pick]))} batches; {len(sel)} genes "
          f"({len(kd)} knockdown targets, {len(set(kd) - set(hvg))} added beyond the HVGs)")
    np.savez_compressed(dst, X=X[:, sel], genes=genes[sel], batch=batch[pick].astype(str),
                        kd_targets=np.array([genes[k] for k in kd]))


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], a[2], a[3], *(int(x) for x in a[4:]))
