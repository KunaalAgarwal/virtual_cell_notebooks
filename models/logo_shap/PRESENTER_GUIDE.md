# Presenter's guide: Control-Cell GRN Review

A companion to the [review page](https://claude.ai/artifact/1zm4X79ZrwXpUazKmA2JRe). It follows the page section by section. For each section it gives the point to make, how to explain it to someone with little background, and what the code does underneath. Part B is a thorough "how can we improve this". A glossary and likely audience questions close it out.

Suggested timing for a 20-minute talk: Part A sections 0–2 in 4 minutes, 3–7 in 7 minutes, 8 (the test) in 4 minutes, 9–10 and Part B in 5 minutes.

---

## The one-paragraph version

> Cells are governed by genes that switch each other on and off. The challenge asks us to predict how three unfamiliar cell types respond when one gene is switched down, given only their normal, unperturbed cells. We tried the natural first idea: learn a "wiring diagram" of which genes influence which from those normal cells, then simulate pulling one wire. We built the pipeline, made it about 100,000× faster, fixed a flaw in the maths, and tested it on the only dataset where we can check the answer. The simulated effects don't match the real ones at all (correlation 0.00), while a trivial "every knockdown looks like the average knockdown" guess scores 0.19. That negative result is the finding: the signal we need lives in real knockdown experiments, not in the co-variation of unperturbed cells.

---

# Part A: walking through the page

## 0. Background to give first

Say this before showing the page.

- **Gene expression.** Every cell carries about 20,000 genes. How strongly each is "on" is its expression level, measured by counting RNA molecules.
- **Single-cell RNA-seq.** Measures expression for thousands of genes in each of thousands of individual cells. The result is a table: one row per cell, one column per gene, each value an integer count (a **UMI** count). Most entries are zero, because a single cell has few molecules of most genes. That sparsity matters later.
- **CRISPRi knockdown.** A guide RNA steers a dead Cas9 protein to one gene and turns it down. Other genes then change because genes regulate each other. The full set of changes is the knockdown's **effect**.
- **Non-targeting controls.** Cells that went through the identical procedure with a guide that matches nothing in the genome. They show the cell's normal state, with the same handling and stress as the real knockdowns.
- **The challenge.** For three new cell types (contexts A, B and C) we receive only control cells. We must output 400 cells × 18,533 genes for each of 300 knockdowns per context, 360,000 cells in total. Arc has run the real experiments and scores us against them. "Zero-shot" means we never see knockdown data from these cells.

---

## 1. Header and the three numbers

**The point:** lead with the result.

| Number | What it means |
|---|---|
| **0.00** | Correlation between what our network predicts a knockdown does and what it really does, averaged over 147 knockdowns in H1. Zero means no relationship. |
| **0.19** | The same score for "predict that this knockdown does what the average knockdown does". A low bar, but a real one. |
| **36 s** | Time to learn a 10,468 × 10,468 edge matrix with the all-genes shortcut, compared with about 40 days one gene at a time. |

Explain Pearson correlation as: +1 means predictions move exactly with reality, 0 means no relationship, −1 means systematically backwards.

---

## 2. "The plan, and where we are"

**The point:** we followed a specific proposal, the `grn_linear_model.pdf` notes, which have four stages.

1. **Learn the edge matrix.** Genes are nodes in a graph. An edge from gene A to gene B means "A's level helps predict B's level", with a weight for strength and a sign for direction (activation vs repression). The full set of weights is an N × N matrix.
2. **Simulate a knockdown.** Lower one node and let the change travel along the edges.
3. **Check against reality.** Compare the simulated changes with real knockdown experiments.
4. **Build a submission.** Apply the predicted changes to each context's control cells and package the file.

The status pills: stages 1 and 2 are built, stage 3 was run and found no signal, stage 4 was never reached because stage 3 failed.

---

## 3. "The data"

**The point:** two datasets with different jobs.

- **2026 contexts A/B/C** are what we're scored on. Each has 18,400 control cells, every one carrying one of 46 non-targeting guides, 400 cells per guide. About 10,500 of the 18,533 genes are detected in at least 5% of cells; the rest are too sparse to model.
- **H1 (2025)** is a human embryonic stem cell line from last year's challenge. It is the only dataset we have with both control cells (38,176) and real knockdowns (150 genes, 183,000 cells). That makes it our test bench.

**What the code does.**
- `extract_h1_controls.py` opens the 15.5 GB H1 file in "backed" mode (reading from disk rather than loading it), keeps the rows where `target_gene == "non-targeting"`, and writes a 0.8 GB controls-only file.
- **Normalization**, applied everywhere before modelling: divide each cell's counts by that cell's total and multiply by 10,000 (so a deeply sequenced cell doesn't look like it has more of everything), then take `log(1 + x)` (so a handful of extremely abundant genes don't dominate). Every "change" on the page is in these log units. A typical knocked-down gene drops by 0.63 in these units, which is roughly a halving.

---

## 4. "Stage 1: LOGO cross-validation with SHAP"

**The point:** to learn the edges, ask for each gene, "can I predict this gene from all the others, and which others matter most?" The answer to "which matter most" is that gene's row of the edge matrix.

Walk the diagram left to right.

1. **Control cells, one context.** The expression table for 18,400 cells.
2. **46 groups.** Cells are grouped by which non-targeting guide they got. The groups are the unit for cross-validation.
3. **Each fold: leave one group out (LOGO).** Train on 45 groups, test on the held-out one, and repeat 46 times so each group is held out once. Testing on cells the model never saw tells us whether it learned something real or just memorized noise.
   - The model is **regression**: target gene ≈ weighted sum of the other genes.
   - **Ridge regression** is linear regression plus a penalty on large weights. With thousands of inputs and only thousands of cells, unpenalized regression would overfit badly; the penalty shrinks weights toward zero and stabilizes them.
   - **LightGBM** is gradient-boosted decision trees. It can capture non-linear, "if gene A is high AND gene B is low" patterns. GRNBoost2 and GENIE3, the standard GRN-inference tools, use this family.
4. **SHAP on the held-out cells.** SHAP values (Shapley values, from cooperative game theory) split each prediction fairly among the inputs: "of this cell's predicted level, +0.3 came from gene X and −0.1 from gene Y." For 400 held-out cells and 2,000 input genes that gives a 400 × 2,000 table.
5. **One row per fold.** The 400 rows are collapsed into one row per gene with three statistics:
   - `mean_abs_shap`: the average magnitude of each gene's contribution. This is the **edge weight**.
   - `signed_shap`: the same magnitude with a sign. The sign is the sign of the correlation between the input gene's value and its SHAP contribution: positive means "more of this gene pushes the target up" (activation), negative means repression. Raw SHAP values average to about zero by construction, so a plain average can't carry the sign.
   - `native_importance`: the model's own measure (|coefficient| for ridge, split gain for trees), kept for comparison.
6. **One row.** The 46 fold rows are averaged. This is the "average the CSV into a single row" step from the original ask.
7. **N × N.** Repeat for every target gene and stack the rows into the edge matrix.

**What the code does** (`logo_shap.py`).
- `load_controls` reads the file and can subsample cells per group for quick runs.
- `normalize` does the depth scaling and log transform.
- `select_features` picks the 2,000 most variable genes detected in at least 5% of cells as inputs, for speed. The target is always removed from its own inputs.
- `make_model` builds `RidgeCV` (which picks its own penalty strength from a grid), `LassoCV` or `LGBMRegressor`.
- `logo_for_target` runs the fold loop for one target using `LeaveOneGroupOut`, and records held-out R² and Pearson per fold.
- `shap_values` uses `shap.LinearExplainer` for ridge (exact) and `shap.TreeExplainer` for LightGBM (exact for trees).
- `fold_row` computes the three statistics for a fold.
- `average_folds` averages the fold rows, optionally weighting by group size.
- `write_outputs` writes `shap_per_fold.csv` (one row per fold), `shap_single_row.csv` (the averaged rows), `edge_matrix_*.csv` (targets × regulators) and `cv_metrics.csv`. It rewrites them after every target so a long run that stops early keeps its work.

**The results table.**
- R² is the share of a gene's variation the model explains. 1 is perfect, 0 is no better than guessing the mean.
- TMSB4X, one of the most abundant transcripts in the cell, is partly predictable (R² 0.31). PRCP and HIRA, which are mostly zeros at the single-cell level, barely are (0.06).
- Ridge's importance rankings are nearly identical across folds (0.95 or more). LightGBM's vary a lot (0.6 to 0.7), so averaging over folds is doing real work for the tree model.
- Ridge and LightGBM agree only moderately on which genes matter (0.6).
- The top "regulators" were histone genes (cell cycle), mitochondrial genes (cell health) and VIM / COL1A1 (differentiation). That is the first hint of the problem in section 6.

---

## 5. "The fast version: every gene in one solve"

**The point:** a piece of linear algebra turns weeks into seconds.

- **The naive way.** Fit one regression per target gene. At about 6 minutes per gene and 10,000 genes, that is about 40 days per network.
- **The shortcut.** All 10,000 ridge regressions share the same inputs, so they are all hidden in one matrix. Standardize the expression table, compute the gene × gene covariance matrix, add the ridge penalty to its diagonal, and invert it once. Column *i* of the inverse, divided by its diagonal entry and negated, is exactly the set of weights for predicting gene *i* from all the others. This is a standard identity for regressing each variable on the rest.
- **Checked.** On context A with 10,468 genes the full matrix took 36 seconds and agreed with scikit-learn's per-gene ridge to 14 decimal places.
- **SHAP for free.** For a linear model SHAP has a closed form: a gene's mean |SHAP| is |weight| × the average distance of that gene from its mean. So the edge weights from section 4 come out of the same computation.

**What the code does** (`v2/grn_ridge.py`).
- `load` reads a context, normalizes it, keeps genes detected in at least 5% of cells, and also returns per-cell library size and mitochondrial fraction.
- `fit` splits cells with `GroupKFold` (5 grouped folds) or full leave-one-group-out, standardizes using training-fold statistics only, and does one **eigendecomposition** of the covariance matrix per fold. In the eigen-basis, adding a penalty λ to the diagonal is trivial, so trying seven λ values costs almost nothing. It picks the λ with the best held-out R² on the first fold (λ = 10 for both H1 and context A).
- Per edge it accumulates across folds: `coef_mean` (the slope: change in target per unit change in regulator, which stage 2 needs), `coef_std` (fold-to-fold spread), `sign_agree` (the fraction of folds that give the edge the same sign as the mean), `shap_mean` and `shap_share` (each fold's |SHAP| column normalized to sum to 1 per target, then averaged, so targets are comparable).
- `top_edges` writes each gene's 20 strongest regulators with sign agreement ≥ 0.8 to a CSV.
- Runtime on H1: about 2 minutes per fold for 10,993 genes and 38,176 cells, about 16 minutes for a 5-fold network. The 1.8 GB `.npz` output holds the full matrices.

---

## 6. "Most co-variation is cell state, not regulation"

**The point:** two genes rising and falling together across cells does not mean one controls the other.

- Cells in the same dish differ in cell-cycle phase, size, how much RNA was captured, and health. Each of those moves thousands of genes at once. A gene that predicts another may just be reading the same clock.
- `remove_cell_state` regresses every gene on a set of nuisance variables and keeps the residual: log library size, mitochondrial fraction, and the top 10 **principal components** (the ten largest shared patterns of variation, found with a randomized SVD). What remains is co-variation not explained by those shared programs.
- **Result:** average held-out R² drops from 0.08 to 0.02 in both H1 and context A. About three-quarters of the predictability was shared cell state. What is left is weak.

This is the first of two reasons the approach fails. Even the "network" that remains after cleaning is mostly noise.

---

## 7. "Stage 2: making the knockdown simulation well-posed"

**The point:** the notes' equation had a flaw that would make every prediction zero, and we replaced it.

- **The notes.** Each gene's expression equals a weighted sum of its neighbours' expression. Written as a matrix equation that is (edges − I) × x = 0.
  - If that matrix is invertible, the only solution is x = 0: every gene switched off.
  - If it is not invertible, there are infinitely many solutions and no way to choose.
  - Either way it doesn't produce a usable expression vector. Yagnik raised the invertibility worry early and it is a real issue, just not the only one.
- **The fix.** Solve for the *change* δ caused by the knockdown instead of the absolute level.
  - Fix the knocked-down gene's change δ_k (for example −0.63).
  - Every other gene's change is the weighted sum of its regulators' changes: δ = e_k · δ_k + Wᵀ δ.
  - Iterate, keeping gene k pinned at δ_k. Each round is one more hop along the edges: direct targets first, then their targets, and so on.
- **Why it is safe.** The iteration converges if the largest eigenvalue magnitude of W is below 1, meaning each hop shrinks the ripple. Measured on context A: 0.92 for the raw network and 0.57 with cell state removed. So (I − W) is invertible and well-conditioned (condition numbers 15 and 2.6). The ridge penalty is what keeps the eigenvalues small; an unpenalized fit could blow up.

**What the code does.** `propagate` in `v2/eval_h1.py` does this for all 147 knockdowns at once: a genes × knockdowns matrix D, updated as D ← D₀ + Wᵀ D with the knocked-down rows re-clamped each step. It was run with 1, 3 and 20 steps.

---

## 8. "The test: do simulated knockdowns match real ones?"

**The point:** this is the honest check, and the answer is no. Spend the most time here.

**Step 1: the ground truth** (`v2/h1_pseudobulk.py`).
- For each of H1's 150 knockdowns, average the normalized expression of its cells and subtract the average of control cells **from the same sequencing batch**. H1 was processed in 48 batches, and batches differ technically; subtracting the matched batch's controls prevents those differences from masquerading as knockdown effects. The result is one **pseudobulk** effect vector per knockdown.
- The 15 GB file never fits in memory, so the script streams it 4,000 cells at a time, normalizes each chunk, and adds each cell into its (knockdown, batch) accumulator with a sparse indicator-matrix multiply.

**Step 2: the simulation.** The network is learned from H1 control cells only, mirroring the zero-shot setup. Each knockdown is pushed through it with `propagate`.

**Step 3: the score.** Pearson correlation between predicted and real change across the network's 10,993 genes, excluding the knocked-down gene itself (its drop is an input, not a prediction). Also reported: correlation on the 100 genes that really changed most, and the fraction of the 50 most-changed genes whose direction is predicted correctly.

**Step 4: the comparisons.**
- **Average knockdown.** Predict that this knockdown does what the other 149 do on average (leave-one-out, so it never sees itself). It contains zero knowledge about the specific gene.
- **Co-expression.** Each gene moves in proportion to its plain correlation with the knocked-down gene across control cells. This is the network idea without the multi-gene modelling.
- **Real vs fixed δ_k.** "Real" uses each knockdown's observed self-drop (an oracle we wouldn't have in practice). "Fixed" uses the H1 median of −0.63, which is what we would assume zero-shot.

**Reading the bar chart.** The average-knockdown baseline is at 0.19. Every network variant is at 0.00 or slightly negative, with and without cell-state removal, for 1, 3 or 20 steps, and for real or fixed δ_k. Co-expression is also 0.00. Adding the network's shift on top of the baseline makes the baseline slightly worse (0.190 → 0.182–0.186).

**Reading the dot plot.** One dot per knockdown; larger dots are knockdowns with bigger real effects. The orange (baseline) dots sit mostly to the right of zero with a median of 0.19. The blue (network) dots are centred on zero. Hovering shows each gene and its three scores.

**The takeaway.** There *is* a real shared response to knockdowns, a generic stress program that the average captures (RNF20, PRDM14 and KDM1A score 0.5 or more under the baseline). But how genes co-vary across unperturbed cells is not the route a knockdown's effect follows. This is the second and decisive reason the approach fails.

---

## 9. "Contexts A, B and C are very different cells"

**The point:** even a working network would not transfer between contexts.

A parallel run of the stage-1 pipeline on each 2026 context found:
- only 541 of each context's 2,000 most variable genes are shared by all three;
- on those shared genes, the regulator rankings barely correlate across contexts (Spearman −0.07 to 0.18), and no top-10 regulator is shared by all three for any target;
- marker genes suggest A is lymphoid (IL32, IKZF1, IL2RG), B epithelial or mesenchymal (MUC16, MMP2), and C epithelial (KRT13, KRT17) with an interferon signal;
- predictability follows the same pattern as H1: TMSB4X is well predicted (R² 0.35–0.64), the other tested genes reach 0.02–0.12.

So any effect we learn elsewhere must be adapted per context, and a network built from H1 stem cells can't be reused as-is.

---

## 10. "What this means, and what to do next"

1. Make the average-knockdown baseline a real submission to get a leaderboard reference.
2. Get the causal signal from real perturbation screens (section B.4).
3. Use each context's control cells to adapt an effect to that context, not to supply the effect.
4. Keep the tooling: the fast fitter, the batch-matched effects and the scorer are reusable.

The caveats on the page are worth saying out loud: the test covered the linear network only; scores are correlations on pseudobulk, not the challenge's official metrics; and the context A/B/C runs were subsampled.

---

# Part B: how can we improve this?

Three levels: things we could do to make the control-cell network better (B.1–B.3), the strategic change the result points to (B.4–B.6), and engineering and process improvements (B.7–B.8). The short answer to "should we keep polishing the network?" is no: co-expression scoring zero means there is little signal to polish. But B.1–B.3 are cheap and would make the negative result airtight.

## B.1 Make the negative result bulletproof

The page tests one network design. Before closing the door, cover the obvious objections.

| Objection | What to run | Cost |
|---|---|---|
| "A non-linear model would capture what ridge misses." | Put the LightGBM edge matrix through `eval_h1.py`. The all-genes shortcut doesn't apply to trees, so use the 2,000-gene subset and 5 folds. | Hours on Colab |
| "Single cells are too noisy; the network is real but buried." | Build **metacells** (average each cell with its 20–50 nearest neighbours) before fitting. This raises per-gene R² substantially in most datasets. | Easy; scanpy has neighbour graphs |
| "The score hides a few genes that do work." | Report the distribution per knockdown and per target gene, not just the mean. Check whether any subset (transcription factors, highly expressed genes) does better than zero. | Minutes |
| "Correlation on pseudobulk isn't the challenge's metric." | Run Arc's `cell-eval` (already a repo dependency) on simulated vs real H1 cells. | Hours |
| "The network was learned on the wrong scale." | Try edges on standardized expression, on raw counts with a Poisson model, or on metacells. | Hours |

If all of these still score about zero, the conclusion is settled.

## B.2 Improvements to stage 1 (edge learning), if we keep it

These come from `REVIEW.md`, updated with what we now know.

- **Use slopes, not mean |SHAP|, as propagation weights.** Mean |SHAP| equals |slope| × how much the regulator varies, which is right for ranking but wrong for simulating. `grn_ridge.py` already outputs slopes (`coef_mean`); `logo_shap.py` doesn't yet. For trees, add `shap_slope = cov(x, SHAP) / var(x)`.
- **Restrict regulators to transcription factors.** About 1,600 human TFs. This is what GRNBoost2/SCENIC do; it removes most cell-state edges by construction and cuts compute 6×.
- **Choose features inside each fold.** `select_features` currently picks the 2,000 variable genes on all cells, a mild leak of test information into training.
- **Standardize before ridge.** The per-gene loop in `logo_shap.py` doesn't, so the penalty is lighter on high-variance genes. `grn_ridge.py` already does.
- **Keep only stable edges.** Threshold on `sign_agree` (e.g. ≥ 0.8) and on `coef_std` relative to `coef_mean`. Unstable edges add noise to propagation.
- **Normalize each fold's row before averaging** (as GENIE3 does) so one fold with a larger scale doesn't dominate. `shap_share` in `grn_ridge.py` does this; `average_folds` in `logo_shap.py` doesn't.
- **Fewer, more meaningful folds.** In the 2026 data all 46 guides are non-targeting and come from one dish, so LOGO-by-guide is close to random K-fold. Five grouped folds give almost the same averaged row at a ninth of the compute. For H1, grouping by sequencing batch (what `grn_ridge.py` does) is the harder, more honest split.
- **Try a sparse model.** Graphical lasso or lasso-per-gene gives a sparse matrix of *conditional* dependencies, which is closer to "direct regulation" than ridge's dense one, and addresses the A → B → C problem noted in `kunaal_ideas.md`.
- **Add prior knowledge.** Mask or weight edges using STRING, TRRUST (TF → target) or GO co-annotation, as GEARS and TxPert do. This is where the control data could contribute: as evidence for or against a known edge rather than as the sole source of edges.

## B.3 Improvements to stage 2 (propagation)

- **Damping and depth.** The clamped iteration already converges. A damping factor (δ ← δ₀ + αWᵀδ with α < 1) or a fixed small number of steps makes it robust to a badly conditioned W from an unpenalized or larger fit.
- **Sign-aware weights.** Propagation should use signed slopes; `signed_shap` from stage 1 is magnitude × sign and is not a slope. Section B.2's first point fixes this.
- **Apply the change per cell, not per mean.** The notes propose scaling each control cell by the ratio of perturbed to control GRN solution. In log space that is "add δ to every control cell". Doing it per cell preserves the realistic cell-to-cell spread that the challenge metrics reward.
- **Handle genes the network doesn't cover.** Sparse genes (about 7,500 of the 18,533) have no edges. They need a default: no change, or the average-knockdown change.

## B.4 The strategic change: learn effects from real knockdown screens

This is the main improvement. The H1 test says control cells don't carry the causal signal. Real knockdown experiments do, and several large public ones exist (listed in the repo's `data_access_guide.md`):

| Dataset | Cells lines | Scale | Why it matters |
|---|---|---|---|
| Replogle 2022 | K562, RPE1 | ~2,000–10,000 knockdowns each | The classic genome-wide CRISPRi screens; h5ad, ready to use |
| Nadig 2025 | HepG2, Jurkat | ~2,400 essential genes each | Liver and T-cell contexts |
| X-Atlas / Orion 2025 | HCT116, HEK293T | 18,903 knockdowns, 8M cells | Largest available; deep sequencing |
| Feng 2024 | 34 iPSC lines | 7,226 knockdowns | Closest biology to H1 stem cells; shows cross-background variation |
| H1 (VCC 2025) | H1 hESC | 150 knockdowns | Already in hand; our test bench |

Concrete approaches, roughly in order of effort:

1. **Shared-response baseline, done properly.** Compute the average knockdown effect in each public screen and in H1. Check how much of it is universal (the same in every cell line) versus line-specific. A universal stress response applied to each 2026 context is the first thing to submit.
2. **Per-gene effect transfer.** For each of the 300 target genes, collect its measured effect in every screen where it was knocked down. Predict the 2026 effect as a weighted combination, with weights set by how similar each screen's control profile is to the target context. The control-cell data is used here for similarity, which is a job it can do.
3. **Gene-feature regression.** Train a model that maps (features of the target gene, features of the context) → effect vector, where gene features come from the screens (its mean effect, which pathways it hits), prior networks (STRING, GO) and expression in the context. This is the GEARS/TxPert family and can be done with plain gradient boosting before touching graph neural networks.
4. **Fine-tune or adapt a published model.** STATE (Arc's own transformer, in `virtual_cell_background`) and others were trained for exactly this task. Running one as-is on the 2026 controls is a strong baseline before any of our own modelling.

Validation for all of these stays the same: hold out H1 (or hold out one cell line from the public screens), predict its knockdowns zero-shot, score with `cell-eval` against baselines.

## B.5 Use each context's control cells for what they are good at

The control data is not useless; it is being asked the wrong question. It is good for:

- **Which genes are expressed.** A gene that isn't expressed in context B can't respond. Mask predicted effects by expression.
- **Expression scale.** Effects learned in log space transfer better when rescaled to each context's own levels and library size (A/B/C have ~20,000 UMIs per cell, H1 has ~55,000).
- **Cell-state composition.** The PCs we regressed out are real biology: cycle phase, differentiation state. The proportion of cells in each state differs between contexts and changes under knockdown. Modelling how a knockdown shifts those proportions may matter more than gene-by-gene edges.
- **Context similarity.** Comparing each context's control profile to the public cell lines picks which screens to borrow from (section B.4, approach 2).
- **Realistic cell-to-cell variation.** The 400 output cells per knockdown should look like real single cells, not 400 copies of a mean. Sampling control cells and shifting them keeps the right sparsity and spread. The repo's `baseline_control_mean.py` already produces the right file format.

## B.6 Scoring: match the challenge before trusting a number

- All results so far are Pearson correlation on pseudobulk effect vectors. The challenge scores on its own metrics (differential-expression overlap, per-gene discrimination, and a mean-absolute-error term; see `notes/scoring.md` and `vcc2026-metrics.pdf`).
- `cell-eval` is in the repo's dependencies. Wrap the H1 test so it writes simulated cells in the submission format and runs `cell-eval` against the real H1 cells. Then every method, baseline included, has a number that means the same thing as the leaderboard.
- Submit the control-mean baseline and the average-knockdown baseline early to learn what the leaderboard scale looks like.

## B.7 Engineering improvements to the code

- **Memory.** `grn_ridge.py` holds several 11k × 11k float32 matrices (about 500 MB each) and peaked near 9 GB on H1. Accumulating in float32, deleting intermediates and chunking the SHAP computation would halve that and allow all ~12,700 genes detected in ≥5% of cells rather than the 10,993 at ≥20%.
- **Column-wise selection in `logo_shap.py`.** Converting to dense per target (`Xn[:, f_idx].toarray()`) is fine for 2,000 genes but not for all genes; the fast path should replace it, with the per-gene loop kept for trees only.
- **Reproducibility.** Pin seeds everywhere (LightGBM, subsampling, randomized SVD all take one; the ridge path is deterministic). Record the λ chosen, the fold assignment and the gene list in each output `.npz`.
- **Config file, not flags.** The two scripts take overlapping command-line flags. A small YAML or dataclass config shared by both would prevent drift.
- **Streaming for the 2026 contexts.** Each context is only ~215 MB, so loading is fine. For the public screens (1.5–9 GB), reuse the chunked reader from `h1_pseudobulk.py`.
- **Tests.** Three cheap unit tests would catch most regressions: the closed-form ridge against scikit-learn on a small random matrix (already done once by hand), `propagate` on a 3-gene toy network with known answer, and `fold_row`'s sign logic on a synthetic positive and negative edge.
- **Put it in the repo.** Nothing is in `virtual_cell_notebooks` yet. Per the repo's conventions, `grn_ridge.py` and `eval_h1.py` belong in `models/` (with a `predict(...)` entry point for the submission path), the notebook in `notebooks/`, and the H1 pseudobulk file tracked with git LFS. This needs a PR, which waits for Yagnik's go-ahead.

## B.8 Process improvements

- **Run on Colab for anything over a few minutes.** The cloud container here has 4 CPUs and 15 GB RAM; the per-gene LightGBM loop took 20+ minutes per gene on full-size 2026 data. The Colab notebook is ready and asks for the vcc key at runtime without storing it.
- **Validate before building.** The biggest saving in this project was running the H1 check before building a submission pipeline. Keep that order for every new idea: baseline score first, then the method, then the submission.
- **Keep the baselines in the table.** Every result table should include the control-mean and average-knockdown rows so improvements are measured against something.
- **Write down what was ruled out.** `RESULTS.md` and this guide do that for the control-cell network. Doing the same for each future method prevents re-running dead ends.

---

## Glossary

| Term | Plain meaning |
|---|---|
| Expression | How strongly a gene is "on" in a cell |
| UMI | One counted RNA molecule; the raw data is UMI counts |
| CRISPRi knockdown | Switching one gene down on purpose with a guide RNA |
| Non-targeting / control cell | Same procedure, guide hits nothing; shows the normal state |
| Pseudobulk | Averaging many single cells into one profile |
| log1p(CP10k) | Our normalization: scale to 10,000 counts per cell, then log(1 + x) |
| GRN / edge matrix | A table of "gene A influences gene B with strength w and sign ±" |
| Regression | Predicting one number from others with a weighted formula |
| Ridge | Regression with a penalty that keeps weights small and stable |
| λ (lambda) | The strength of that penalty; chosen by held-out performance |
| LightGBM | Gradient-boosted decision trees; captures non-linear rules |
| Cross-validation / LOGO | Train on most groups, test on the one left out, repeat |
| Fold | One train/test split in cross-validation |
| R² | Share of variation explained (0 = none, 1 = all) |
| Pearson r / Spearman | Correlation of values / of rankings (−1 to 1) |
| SHAP | Splits a prediction fairly into each input's contribution |
| Principal components (PCs) | The largest shared patterns of variation across all genes |
| Residual | What is left of a gene after removing what the nuisance variables explain |
| Eigenvalue / spectral radius | How much a matrix stretches vectors; below 1 means ripples shrink |
| Condition number | How numerically touchy a matrix inverse is; small is good |
| Metacell | A small group of similar cells averaged to reduce noise |
| Transcription factor (TF) | A gene whose protein directly turns other genes on or off |
| cell-eval | Arc's official scoring package for the challenge |

## Likely questions

- **"Isn't 0.19 also bad?"** It is modest, but it is the bar. Mean-effect baselines are often hard to beat on perturbation benchmarks, which is why the challenge exists.
- **"Maybe a non-linear model would work?"** Plain co-expression also scored about zero, so there is little signal in control cells for any model to find. We haven't run LightGBM through the knockdown test, and B.1 says to do it for completeness.
- **"Why group by guide or batch instead of random splits?"** Cells from the same group share technical quirks. Holding out whole groups tests whether results generalize to new cells rather than to near-duplicates.
- **"Is the matrix invertible?"** The PDF's equation needed W − I invertible, and even then only gave the zero solution. The change-based version needs I − W invertible, which we measured: it is, comfortably.
- **"Is this the official metric?"** No. It is correlation on averaged profiles. The challenge uses Arc's `cell-eval` metrics; B.6 says to switch to them before trusting any number.
- **"Was any knockdown data used to build the network?"** No. Only control cells, mirroring the zero-shot setting. Knockdown data was used only to score it.
- **"So was this wasted?"** No. We now have a fast network fitter, batch-corrected ground truth for 150 knockdowns, a scorer with baselines, and a clear answer about where the signal is not. That directs the next month of work.
- **"What would you submit tomorrow?"** The average-knockdown response from H1, applied to each context's control cells with expression masking, in the repo's existing submission format. Then iterate from B.4.
