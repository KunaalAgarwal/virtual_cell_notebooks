# VCC 2026 submission (final phase D, E, F; validation A, B, C)

Applies the public-screen effect-borrowing method (`models/public_screens/`, RESULTS.md Round 4)
to the targets of a controls bundle and packages the result in the official 2026 format.

## Final phase

The final phase scores new contexts D, E, F against a new perturbation list that is published
only inside its own bundle (`vcc datasets download controls-final`, `vcc` 0.3.0 or newer). A file
built from the validation bundle uploads and then fails with "Submission is missing predictions
for context(s): D, E, F". Nothing here assumes A/B/C any more:

- `bundle.py` reads the contexts from the bundle's `manifest.json` (as `vcc prep` does), the
  targets from `pert_counts.csv` (per context if it has a `context` column) and finds each
  context's control `.h5ad` (`context_<X>.h5ad`, else by its obs `context` label).
- `package_vcc.py` takes the contexts from the same manifest through `vcc.panel`, and writes the
  panel into `meta.json` (`contexts`, `perts_fingerprint`, `panel_id`) exactly as `vcc prep` 0.3.0
  does, so `vcc submit` refuses a file made for the other phase before uploading it.
- The whole K562 screen ships (`data/k562_gwps_all/`), because the final targets are unknown
  until the bundle is downloaded.

Tested 2026-10-22 on a mock D/E/F bundle (validation controls relabelled, 29 K562 targets not in
the validation panel, panel id `mock-final`): both presets build, every context passes
`vcc prep --dry-run` 0.3.0, and the packaged `meta.json` matches the one plain `vcc prep` writes.
In the final phase only a team's last final entry counts, and no scores are shown.

## Format (from the vcc CLI 0.2.2 and cell_eval2 0.18.0 docs)

- One `.vcc` for all three contexts, made from one `.h5ad`: `obs.target_gene` (gene symbol),
  `obs.context` (A/B/C exactly as in the control files), all 18,533 genes in `gene_names.csv` order.
- Exactly 400 cells per target per context (360,000 cells), no non-targeting cells.
- Raw integer counts, at most 1e6 per cell, sparse, at most 4.75e9 stored entries.
- Scoring (`cell_eval2` `vcc2026` profile, counts space): six metrics, each rescaled so 0 is the
  context's mean perturbation response and 1 is a split-half replicate; Overall is their mean.
  `pds_cosine` is scale-free; `expr_mse_unbiased_capped_norm` and `de_wilcoxon_lfc_nmae` penalize
  effects that are too large.

## Inputs shipped in `data/`

- `h1_pseudobulk.npz`: log1p(CP10k) effects of the 2025 H1 knockdowns.
- `shared_response.npz`: cross-line shared knockdown response (see "Shared response").
- `k562_gwps_all/` (111 MB in 5 .npz files): all 11,258 rows of the K562 screen as int16
  (z x 1000 after the loader's clip to +-10 and NaN -> 0), byte-shuffled and zlib-compressed.
  `python k562_pack.py unpack --parts data/k562_gwps_all --out K562_gwps_rebuilt.h5ad` rebuilds
  the full `.h5ad` (original obs, var, row order) and copies in the exact rows of the subset below.
  Validation-panel folds are bit-identical to the original file (`c1_dev`, `d3_dev`, context A).
  On 300 other K562 targets the core prediction agrees to max |log change| 0.004 (`d3_dev`
  without the wide calls) and 0.03 (`c1_calibrated`); with the wide calls a few genes whose
  per-target deviation is near zero can flip sign (row correlation >= 0.996 for `d3_dev`, >= 0.91
  for `c1_dev`). With the old subset alone, 299 of those 300 targets would have had no K562 profile.
- `k562_gwps_bulk_subset.h5ad` (31 MB): the rows of Replogle et al. 2022
  `K562_gwps_normalized_bulk_01.h5ad` (figshare article 20029387, file 35773217) that the build
  reads: the 300 panel targets, the 2025 H1 knockdowns and the non-targeting controls (996 of
  11,258 rows, all 8,248 readout genes). `build_submission.py` only uses those rows (correlation
  similarity, the target's own profile, strength metadata), so predictions are bit-identical to
  the full 375 MB file (checked on context A for `c1_dev` and `d3_dev`, max |fold difference| 0).
  Added 2026-10-19 after figshare started answering downloads with HTTP 500 JSON error pages,
  which `wget -q -O` saved as the .h5ad and the build then failed on with
  "file signature not found". Since the final phase the notebook uses `k562_gwps_all/` instead
  (this file still supplies its exact rows) and has no figshare download at all.

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
| `c1_wide` | `c1_calibrated` + cross-line shared response called on ~7,000-8,000 genes (see below) | per-target part 2.0 |
| `c1_dev` (default) | `c1_wide` genes, but each target gets its own sign (its predicted change minus the panel average), 4 SE | per-target part 2.0 |
| `d3_dev` | the target's own K562 profile x3, no H1 terms, plus `c1_dev`'s per-target wide calls (see Leaderboard 2, "Where the top teams gain") | per-target part ~4.5 |
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

## Leaderboard 2 (2026-10-18) and the `c1_dev` preset

`c1_wide` scored overall -0.0295 (rank 1139, from -0.1945 and 1283):

| metric | c1_calibrated raw / scaled | c1_wide raw / scaled |
|---|---|---|
| pds_cosine | 0.549 / +0.106 | **0.500 / 0.000** |
| mse (capped norm) | 1.156 / 0.000 | 17.683 / 0.000 (clamped) |
| sig_jaccard | 0.004 / -0.071 | 0.033 / +0.007 |
| lfc_nmae | 0.984 / +0.028 | **1.091 / -0.152** |
| direction_fidelity_yield_raw | 0.135 / **-1.262** | 0.511 / -0.005 |
| direction_reach_raw | 0.107 / +0.031 | 0.055 / -0.027 |

Fidelity was fixed as intended. The cell builder works as designed: the official DE test calls
the intended direction on 99.5% of genes. The losses come from the shared sign pattern itself.
The same ~8,000-gene shift was added to every target, and it is much larger in norm than the
per-target part. So every predicted change pointed the same way, and cosine discrimination
collapsed to a common ranking, which reads exactly 0.5. The shared directions were only 51%
right in A/B/C, so the extra error on real DE genes raised nmae. The shared genes also crowd the
top of the significance ranking, which lowers reach.

`c1_dev` keeps the same genes at a slightly smaller callable size (4 SE instead of 5). Each target's sign
on a gene is the sign of its own predicted change minus the panel-average predicted change.
Fidelity needs only enough calls with better-than-random directions, while per-target signs keep
the targets apart. On H1 with the same c1 machinery (`diagnostics/h1_variants.py`: pds over 136 H1
knockdowns with K562 profiles, fidelity and nmae over the 40 with DE tables, shared response without H1):

| variant | pds_cosine | fidelity raw | nmae (approx.) |
|---|---|---|---|
| c1_calibrated | 0.627 | 0.486 | 0.856 |
| c1_wide (shared signs, 5 SE) | 0.531 | 0.558 | 0.969 |
| **c1_dev** (per-target signs, 4 SE) | **0.751** | 0.569 | 0.933 |
| per-target + 0.3 x shared, 4 SE | 0.740 | 0.574 | 0.931 |

The H1 proxy reproduces both leaderboard moves of c1_wide (pds down to chance, nmae +0.11).
With the official DE test on 12 real targets in context C, `c1_dev` gets 7,538 genes called per
target, 99.6% in the intended direction, median |log2FC| 0.29 (c1_wide: 8,010 genes, 0.32), and
only 42% agree with the shared sign. Expected on A/B/C, not measured: fidelity near 0 as now, pds
back above c1's +0.11, nmae a little better than c1_wide, so overall slightly above 0.

### Where the top teams gain (2026-10-18)

Top-10 entries (overall 0.31-0.47) score pds 0.75-0.88 (raw 0.84-0.90), mse 0.44-0.79 (raw
0.24-0.55), reach 0.25-0.46 and nmae 0.23-0.44, while their fidelity (0.03-0.17) and jaccard
(about 0) sit near the baseline. With RPE1 (Replogle et al. 2022) as a clean stand-in for a new
cell line and the K562 genome-wide screen as the source (`diagnostics/xline_pds.py`,
`pipeline_on_rpe1.py`, `rpe1_scale.py`; lfc = log1p(z x sd / mean)):

| prediction, RPE1 panel of 300 | strongest KDs | random | random, 55% K562-significant (vcc-like) |
|---|---|---|---|
| same gene's K562 profile only | 0.88 | 0.80 | 0.78 |
| our c1 recipe (H1 mean + borrowed H1 + K562 profile) | 0.84 | 0.69 | 0.65 |
| H1 mean + K562 profile | 0.86 | 0.71 | 0.66 |
| borrowed H1 + K562 profile | 0.87 | 0.75 | 0.72 |

The same gene's profile in another line reaches the top teams' pds. Both H1 terms lower it.
Scaling the K562 profile x3 leaves pds unchanged and improves the nmae proxy from 0.92 to 0.78
(x4: 0.75); it raises the raw mse ratio, which the leaderboard clamps at 1 for every entry above
the baseline, ours included. The mse gains of the top teams need the context's own shared
response: RPE1's own mean explains 43% of its response variance, but K562's mean ported over
explains 6% (cosine 0.28), and adding any large common vector costs pds (0.88 -> 0.69-0.72 even
with the oracle shared response). How the leaders get mse 0.24-0.55 in a line without
perturbation data is open; private data or leaderboard tuning over many submissions would both
do it.

## Files

- `build_submission.py`: predictions and cells, streamed to one `.h5ad` (~3 GB RAM).
- `bundle.py`: contexts, targets and control files of an unzipped controls bundle.
- `k562_pack.py`: packs / rebuilds the full K562 screen shipped in `data/k562_gwps_all/`.
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
