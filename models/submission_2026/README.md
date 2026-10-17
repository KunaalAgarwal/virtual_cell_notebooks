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
| `c1_calibrated` | Round 4 "tuned for Pearson": C 1, m 20, K562 term, 1% gate | 2.0 |
| `c8_strength` | Round 4 recommended default: C 8 x strength^2, m 10, K562 term, 1% gate | 4.8 (15x at p90) |
| `c1_wide` (default) | `c1_calibrated` + cross-line shared response called on ~7,000-8,000 genes (see below) | per-target part 2.0 |
| `c4_wide` | `c1_wide` with the per-target part x4 (protects pds against the common shift) | per-target part ~8 |

`c8_strength` was chosen for cell-eval 0.8's L1 discrimination, which rewards large effects.
The 2026 discrimination metric is cosine, so amplitude no longer helps it, while two of the six
metrics charge for overshoot. Submitting both is the direct test.

## Leaderboard 1 (2026-10-14) and the `_wide` presets

`c1_calibrated` scored overall -0.1945 (rank 1283). The leaderboard shows each metric's scaled
score (0 = the context's real mean response, 1 = a split-half replicate) before its raw value:

| metric | raw | scaled |
|---|---|---|
| pds_cosine | 0.549 | +0.106 |
| expr_mse_unbiased_capped_norm | 1.156 | 0.000 (clamped to [0, 1]) |
| de sig_jaccard | 0.004 | -0.071 |
| de lfc_nmae | 0.984 | +0.028 |
| de direction_fidelity_yield_raw | 0.135 | **-1.262** |
| de direction_reach_raw | 0.107 | +0.031 |

Direction fidelity alone costs about 0.21 of overall. In cell_eval2 0.18.0
(`metrics/direction.py`) it is `k / max(n_pred, N_conf)` per target: correct-direction calls
among the genes the prediction calls DE, divided by the larger of the number called and the
number the real knockdown changes. Calling random directions on at least `N_conf` genes reads
about 0.5, which is where the mean-response baseline sits. The official DE test on our cells
(wilcoxon on CPM+log1p, 400 cells vs the 18,400 controls, per-target BH, genes >= 5 CPM) calls a
median of **5 genes** per target for `c1_calibrated`, so it was charged for almost every real
DE gene. The same formula on our H1 cell-eval tables gives 0.146 for the mean-effect method
(15 genes called, 72% right), matching the leaderboard.

The `_wide` presets add the cross-line shared knockdown response (`shared_response.py`: the
mean knockdown profile of K562 and RPE1 (Replogle et al. 2022), HepG2 and Jurkat (Nadig et al.
2025) and H1, each scaled by its median) on the 8,000 genes with the largest |score|. Each gene
moves by 5 standard errors of a 400-cell CPM mean (at most 50%), the smallest change the
official test reliably calls; the target-specific prediction is kept wherever it moves a gene
more. Increases add Poisson counts at the gene's depth-scaled mean rate, because a rank test
cannot see an increase that only scales the non-zero counts of a mostly-zero gene (measured:
0-10% of such increases were called before, 96-100% after).

Checks (scripts in `diagnostics/`, run on the shared project data):

* H1, real DE tables: calling the top N genes of a shared response with its sign gives raw
  fidelity 0.38 (N 1,000), 0.43-0.48 (4,000) and 0.50-0.53 (8,000) from any other line or
  their average, against 0.15 for what we submitted; H1's own mean reaches 0.57 (`sim_topn.py`).
* Leave one line out over H1, K562, RPE1, HepG2, Jurkat: the average of the other lines gets
  the direction of a held-out line's individual knockdowns right on 53-66% of its top 4,000
  genes, as good as the best single source (`xline_sign.py`).
* Official DE (cell_eval2 0.18.0, pdex backend) on 12 real targets per context: `c1_wide`
  calls a median of ~7,000 genes in context B with 99.5% in the intended direction, median
  |log2FC| 0.34 (`e2e.py`). Stored counts grow ~2%.
* Cosine discrimination, H1 pseudobulk proxy: the common shift costs about 0.016 raw pds at the
  current per-target size and about 0.003 at x4 (`pds_sim.py`). That is the reason for `c4_wide`.

Expected effect, not measured on A/B/C (the leaderboard is the only truth there): fidelity near
or slightly above the baseline (scaled about -0.1 to +0.2 instead of -1.26), so overall around 0,
with jaccard likely up and pds slightly down.

## Files

- `build_submission.py`: predictions and cells, streamed to one `.h5ad` (~3 GB RAM).
- `package_vcc.py`: runs the official `vcc prep --dry-run` on each context and the CLI's own
  label checks on the combined file, then writes the `.vcc` the way `vcc prep` does. Needed
  because `vcc prep` holds the whole ~16 GB matrix (about 25 GB peak); on a High-RAM machine,
  plain `vcc prep pred.h5ad -g gene_names.csv --perts pert_counts.csv -o pred.vcc` does the same.
- `shared_response.py`: builds `data/shared_response.npz` (signed cross-line shared response per gene).
- `diagnostics/`: the checks behind the `_wide` presets (paths point at the project's shared data folder).
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

Arc Institute cell_eval2 0.18.0 source (`metrics/direction.py`, `catalog.py`, `configs/vcc2026.yaml`) for the 2026 metric definitions.

Replogle et al. 2022, Cell (K562 genome-wide Perturb-seq); Nadig et al. 2025 / Arc VCC 2025 H1
data; Arc Institute vcc-cli and cell_eval2 (vcc2026 metric specification, 2026-08-19).
