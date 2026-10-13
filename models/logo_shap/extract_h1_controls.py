"""Pull the non-targeting control cells out of the 2025 H1 training .h5ad.

Writes an .h5ad shaped like a 2026 context file (raw counts, obs: target_gene, ntc_id,
batch), so logo_shap.py runs on it the same way. ntc_id is copied from guide_id.

    python extract_h1_controls.py data/atlas_cache/vcc2025/train/adata_Training.h5ad h1_controls.h5ad
"""

import sys

import anndata as ad


def main(src, dst):
    adata = ad.read_h5ad(src, backed="r")
    mask = (adata.obs["target_gene"] == "non-targeting").to_numpy()
    ctrl = adata[mask].to_memory()
    ctrl.obs["ntc_id"] = ctrl.obs["guide_id"].astype(str)
    ctrl.write_h5ad(dst, compression="gzip")
    print(f"{ctrl.n_obs} control cells x {ctrl.n_vars} genes, "
          f"{ctrl.obs['ntc_id'].nunique()} guides, {ctrl.obs['batch'].nunique()} batches -> {dst}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
