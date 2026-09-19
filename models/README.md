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
