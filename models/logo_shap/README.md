# LOGO CV + SHAP exploration

## Summary

We built and tested the edge-learning step from `grn_linear_model.pdf`. For each target gene we fit a regression model that predicts its expression from the other genes, using **control cells only**. Leave-one-guide-out CV checks that the model generalizes. SHAP on the held-out cells measures how much each other gene contributes, and the fold results are averaged into **one row per target gene**. Stacking those rows over many targets gives the N x N edge matrix (the GRN) that the PDF then uses to simulate knockdowns.

It runs end to end. The test covered 3 genes on the H1 2025 controls with ridge and LightGBM, 31 folds each.

## Files

| File | What it is |
|---|---|
| `logo_shap.py` | The pipeline as a plain module, also runnable from the command line. |
| `logo_shap_colab.ipynb` | Colab notebook. Choose `DATA = "vcc2026"` (asks for your vcc key at runtime) or `"h1"` (public data). |
| `extract_h1_controls.py` | Pulls the 38,176 non-targeting cells out of the 15.5 GB H1 training file. |
| `h1_test_run/shap_per_fold.csv` | One row per (target, model, fold, stat). The columns are the ~2,000 regulator genes. |
| `h1_test_run/shap_single_row.csv` | The single averaged row per (target, model, stat). |
| `h1_test_run/edge_matrix_*.csv` | The same rows reshaped to targets x regulators. |
| `h1_test_run/cv_metrics.csv` | Held-out R² and Pearson for each fold. |

## The method, step by step

1. **One set of cells.** We use the control cells of a single context: H1 here, or one of 2026 contexts A, B or C.
2. **Normalize.** Counts are scaled to 10,000 per cell, then log1p is applied.
3. **Partition into groups by guide.** Every control cell received one of several non-targeting guides. In 2026 that is 46 guides x 400 cells, which is exactly the 46 groups in the PDF. In H1 it is 31 guides with 483 to 2,749 cells each. Grouping by guide asks whether the learned relationships hold for cells the model never saw from a different guide, which approximates a different well or batch.
4. **Regulators.** To keep runs fast we use the 2,000 most variable genes that are detected in at least 5% of cells. The target is always removed from its own inputs. Setting `N_TOP_FEATURES=None` uses every detected gene, about 10,000 in 2026.
5. **LOGO CV.** We train on all guides but one, predict the held-out guide's cells, and record R² and Pearson.
6. **Models.**
   - Ridge: the linear model the PDF describes. It uses exact linear SHAP.
   - LightGBM: tree boosting in the style of GRNBoost2/GENIE3. It uses TreeSHAP.
   - Lasso is also available.
7. **SHAP.** For each fold this gives a cells x regulators matrix of SHAP values.
8. **Collapse each fold to one row.** Three statistics are computed per regulator:
   - `mean_abs_shap`: the mean of |SHAP| over the held-out cells. This is the edge weight.
   - `signed_shap`: `mean_abs_shap` multiplied by the sign of corr(regulator expression, its SHAP). A positive value means more of the regulator pushes the target up (activation). A negative value means it pushes the target down (repression). We need this because the plain mean of SHAP over cells is about 0 by construction, so it cannot carry the sign.
   - `native_importance`: |coefficient| for ridge, total split gain for LightGBM.
9. **Average the folds into a single row.** We take the unweighted mean over folds of each statistic. The `WEIGHT_BY_SIZE` option weights folds by held-out cell count instead, which matters for H1 because its guides differ in size.

## Test results (H1, 200 cells per guide, 3 targets)

| Target | Ridge R² / r | LightGBM R² / r | Ridge vs LGBM rank agreement (Spearman) |
|---|---|---|---|
| TMSB4X | 0.31 / 0.57 | 0.26 / 0.52 | 0.64 |
| PRCP | 0.06 / 0.25 | 0.04 / 0.22 | 0.61 |
| HIRA | 0.06 / 0.25 | 0.04 / 0.21 | 0.58 |

- **Highly expressed genes are predictable and lowly expressed ones mostly are not.** TMSB4X is a very abundant transcript, so it is well predicted. PRCP and HIRA are mostly sparse counts, so the signal is largely noise.
- **Ridge is very stable across folds and LightGBM much less so.** For ridge, a typical fold's importance ranking correlates about 0.95 with the averaged row, and about 18 of the top 20 genes are the same. For LightGBM the figures are about 0.65 and 15 of 20. Averaging over folds is what turns LightGBM's noisy importances into a usable row.
- **The two models only partly agree on which genes matter** (Spearman about 0.6).
- **Many top regulators reflect cell state rather than direct regulation.** They include histones (HIST1H1B, cell cycle), mitochondrial genes (MT-ND4, cell quality), and VIM/COL1A1/TAGLN (a differentiation axis). This is expected: in unperturbed cells, co-variation is driven mostly by shared programs. That is the main limitation of learning edges from controls alone.

## The data

| | H1 2025 (training) | 2026 contexts A, B, C (what we're scored on) |
|---|---|---|
| Cells | 221,273 in total, of which 38,176 are non-targeting controls; the rest are knockdowns of 150 genes | 18,400 per context, all controls |
| Guides | 31 non-targeting guides of unequal size | 46 non-targeting guides x 400 cells |
| Genes | 18,080 (18,077 shared with the 2026 list) | 18,533 |
| Depth | median about 54,600 UMI per cell | median about 20,000 UMI per cell, about 6,000 genes detected |
| Access | public GCS bucket, no key | `vcc datasets download controls` (needs your key); or `data/validation_2026/` once downloaded |

- The **2026 task** is to predict, for each of 3 unseen cell lines, 400 cells x 18,533 genes for each of 300 gene knockdowns: 360,000 cells in total. Arc measures the real knockdowns and we never see them.
- **H1 is the only data where we observe both control cells and knocked-down cells** from the same experiment. It is therefore where we can check whether an edge matrix learned from controls predicts real knockdown effects.
- H1 is sequenced at about 2.7x the depth of 2026. Library-size normalization handles most of that difference, but sparse genes will be noisier in 2026.

## How this relates to the task

The PDF's plan is GRN_control → knock out node → GRN_perturb → rescale control cells by the ratio. This exploration builds the **edge matrix** part of that plan for one context. To finish the loop we still need to:

1. Run all ~300 target genes for a context to fill in the rows. This is about 6 minutes per gene per model on this CPU with 2,000 regulators, so it means Colab chunks.
2. Run the PDF's propagation (remove the node, re-solve) to get a predicted shift for every gene.
3. **Validate on H1.** Build the edge matrix from H1 controls, simulate the 150 H1 knockdowns, and score against the real H1 knockdowns with `cell-eval` (already in the repo's dependencies). That is the only honest check before we spend submissions.

## Other methods to run or refine

**Cheaper and better edge learning**
- *Multi-output ridge.* One solve gives every target at once, and for linear models SHAP is exactly coef x (x − mean). That should bring the full 10k x 10k matrix down from days to minutes (an estimate, not yet timed), making it the fastest route to a full GRN.
- *GRNBoost2 / GENIE3 (arboreto).* These are purpose-built versions of the tree approach and are in the background papers. Restricting regulators to transcription factors (about 1,600 genes) cuts cost and removes many cell-state edges.
- *Remove the cell-state axes first.* Regress out library size, cell-cycle scores and mitochondrial fraction, or the top PCs. Without this, edges are dominated by shared programs.
- *Graphical lasso / partial correlation.* This is a direct estimate of conditional dependencies and helps with the A → B → C problem noted in `kunaal_ideas.md`.

**Better validation**
- Group by **batch** instead of guide for H1 (48 GEM batches). Batch effects are likely stronger than guide effects (not yet measured), so this should be a harder and more realistic split.
- Repeat across contexts A/B/C and compare the edge matrices. Edges that are stable across contexts are the ones to trust for zero-shot prediction.

**Using perturbation data, not just controls**
- Edges learned from controls are correlational. Knockdown datasets show real causal effects: H1 (150 genes), Replogle K562/RPE1, Nadig HepG2/Jurkat, X-Atlas HCT116/HEK293T, and Feng iPSC (the closest to H1). These can calibrate or replace control-learned edges, for example by fitting a model that maps (edge-matrix features of the target gene, context control profile) → knockdown effect.
- Add prior networks (STRING, GO, TRRUST) as features or as a mask on allowed edges, as in GEARS/TxPert.

**Baselines to beat first**
- The repo's control-mean baseline (`models/baseline_control_mean.py`).
- "Knocked-down gene goes down, nothing else changes."
- The mean knockdown effect across all H1 perturbations, applied to each context. These are often hard to beat on the challenge metrics (see `notes/scoring.md`).
