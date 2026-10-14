# VCC 2026 validation submission (contexts A, B, C)

Applies the public-screen effect-borrowing method (`models/public_screens/`, RESULTS.md Round 4)
to the 300 targets of panel `vcc2026-val-1` and packages the result in the official 2026 format.

## Format (from the vcc CLI 0.2.2 and cell_eval2 0.18.0 docs)

- One `.vcc` for all three contexts, made from one `.h5ad`: `obs.target_gene` (gene symbol),
  `obs.context` (A/B/C exactly as in the control files), all 18,533 genes in `gene_names.csv` order.
- Exactly 400 cells per target per context (360,000 cells), no non-targeting cells.
- Raw integer counts, at most 1e6 per cell, sparse, at most 4.75e9 stored entries.
- Scoring (`cell_eval2` `vcc2026` profile, counts space): six metrics, each rescaled so 0 is the
  context's mean perturbation response and 1 is a split-half replicate; Overall is their mean.
  `pds_cosine` is scale-free; `expr_mse_unbiased_capped_norm` and `de_wilcoxon_lfc_nmae` penalize
  effects that are too large.

## How cells are made

Each predicted cell is a real control cell of the same context, with counts scaled gene by gene
by exp(predicted log change): binomial thinning for decreases, added Poisson counts for
increases. Counts stay integers, depth and noise stay realistic, and stored entries stay at the
controls' density (~2.06e9 in total).

## Presets (`build_submission.py --preset`)

| preset | what | median effect norm vs real H1 (2.0) |
|---|---|---|
| `c1_calibrated` (default) | Round 4 "tuned for Pearson": C 1, m 20, K562 term, 1% gate | 2.0 |
| `c8_strength` | Round 4 recommended default: C 8 x strength^2, m 10, K562 term, 1% gate | 4.8 (15x at p90) |

`c8_strength` was chosen for cell-eval 0.8's L1 discrimination, which rewards large effects.
The 2026 discrimination metric is cosine, so amplitude no longer helps it, while two of the six
metrics charge for overshoot. Submitting both is the direct test.

## Files

- `build_submission.py`: predictions and cells, streamed to one `.h5ad` (~3 GB RAM).
- `package_vcc.py`: runs the official `vcc prep --dry-run` on each context and the CLI's own
  label checks on the combined file, then writes the `.vcc` the way `vcc prep` does. Needed
  because `vcc prep` holds the whole ~16 GB matrix (about 25 GB peak); on a High-RAM machine,
  plain `vcc prep pred.h5ad -g gene_names.csv --perts pert_counts.csv -o pred.vcc` does the same.
- `sanity_checks.py`: label/count/integer checks, read-back of written vs intended changes, and
  the 13 targets also knocked down in H1.
- `data/h1_pseudobulk.npz`: log1p(CP10k) effects of the 150 H1 knockdowns (from the public 2025
  VCC training data, `models/public_screens/h1_pseudobulk.py`).
- `notebooks/submission_2026_colab.ipynb`: end to end on Colab, asks for the VCC key at runtime.

## Submit

```bash
uv tool install vcc-cli
vcc login --token-stdin          # or: export VCC_TOKEN=... on a headless machine
vcc submit vcc2026_val_c1_calibrated.vcc -m "screen-borrow c1 calibrated"
vcc status <entry-id> --wait
```
Two scored submissions per team per day; one in flight at a time.

## References

Replogle et al. 2022, Cell (K562 genome-wide Perturb-seq); Nadig et al. 2025 / Arc VCC 2025 H1
data; Arc Institute vcc-cli and cell_eval2 (vcc2026 metric specification, 2026-08-19).
