# models/

Model code as plain `.py` files. **No notebooks here.**

Every model script must:

1. Expose `predict(context_paths, gene_names_path, perts_path, out_path)` so it can be imported and called from a notebook.
2. Be runnable standalone via argparse, so it runs headless without a Jupyter kernel.

```bash
uv run python models/baseline_control_mean.py \
  --controls-dir data/validation_2026 \
  -o pred.h5ad
```

Exploration and plots belong in `notebooks/`.

## Exploration folders

`logo_shap/` and `grn_ridge/` hold exploration pipelines, not submission models, so they do not expose `predict()`:

- `logo_shap/`: leave-one-group-out CV + SHAP on control cells (ridge / LightGBM per target gene), averaged to one row per gene. See `logo_shap/README.md`, `REVIEW.md` and `PRESENTER_GUIDE.md`. The Colab runner is `notebooks/logo_shap_colab.ipynb`.
- `grn_ridge/`: fast closed-form all-genes ridge network, propagation, and the H1 knockdown evaluation (Pearson scorer, leave-one-out effect transfer, cell-eval harness). Findings are in `grn_ridge/RESULTS.md`; small result tables are in `grn_ridge/results/`.
