# Where this sits, and what to improve

## Where it sits in the plan

Your PDF has three stages:

| Stage | PDF | Status |
|---|---|---|
| 1. Learn the edge matrix | Model N_i = f(N_-i). Split 18,400 cells into 46 guide groups and run LOGO-CV. SHAP each held-out group (400 rows), aggregate the 46 folds into one N-dim vector, and repeat for every gene to get an N x N matrix. | **Built, for one gene at a time.** The test used 3 genes, about 2,000 regulators rather than 18,533, and H1 controls. The same code runs on 2026 context A/B/C. |
| 2. Simulate a knockdown | Solve for GRN_control, remove the knocked-down node and its edges, then re-solve for GRN_perturb. | Not started |
| 3. Apply it to cells | perturb_i = control_i x (GRN_perturb / GRN_control), giving 400 cells x 300 targets x 3 contexts. Package with `vcc prep`. | Not started (the repo's control-mean baseline already handles the output format) |

How it fits the challenge blurb: the approach only needs each context's own control cells. That makes it zero-shot by construction, because we learn a separate edge matrix for A, B and C and never need their knockdown data. The open question is whether edges learned from unperturbed cells predict knockdown effects. H1 is the one dataset where we can check that, since it has controls and 150 real knockdowns.

## Improvements, most important first

**1. Fit every gene at once for ridge (measured: about 40 days becomes about 30 minutes).**
Ridge regressions of each gene on all the others can be read directly off one matrix inverse. With P = (XᵀX + λI)⁻¹, the coefficients for gene i are −P[:, i] / P[i, i].
- I ran this on 2026 context A with all 10,468 genes detected in at least 5% of cells. The full 10,468 x 10,468 edge matrix took **36 s** and matched sklearn's per-gene Ridge to about 1e-14.
- The current loop takes about 6 minutes per gene, so 10k genes would take about 40 days.
- Linear SHAP is exact (|β| x mean|x − μ|), so mean |SHAP| for every edge comes almost for free. With 46 LOGO folds the whole stage-1 matrix takes about 30 minutes.
- LightGBM stays useful as a check on a subset of genes.

**2. The PDF's solve step needs a small fix.**
"(Edge matrix) x = 0" has only the trivial solution x = 0 unless the matrix is singular. So GRN_init and GRN_perturb would both come out as zero. The working version propagates the *change* instead:
- set the knocked-down gene's change δ_k (for example, its observed CRISPRi log fold change),
- compute δ = (I − Wᵀ)⁻¹ δ_k, or iterate δ ← Wᵀδ for a few damped steps,
- add δ to each control cell in log space.

This is the same idea as removing the node and re-solving, but it is well posed.

**3. Use slopes, not mean |SHAP|, as propagation weights.**
- Mean |SHAP| ≈ |β| x (how much the regulator varies). That is good for *ranking* edges, but it is not "change in target per unit change in regulator", which is what propagation needs.
- For ridge, use β directly. For trees, add a `shap_slope` stat = cov(x, SHAP) / var(x).
- Keep mean |SHAP| for ranking and for choosing which edges to keep.

**4. Make the single row more trustworthy.**
- Add fold stability stats alongside the mean: the standard deviation, and the fraction of folds that agree on the sign. Then keep only edges with consistent signs.
- Normalize each fold's row to sum to 1 before averaging, as GENIE3 does, so one fold with a large scale doesn't dominate.
- Report each edge as a share of that target's total importance, so rows are comparable across targets.

**5. Fewer folds; guide is a weak grouping.**
- In 2026 the 46 guides are all non-targeting and come from the same cell line, so LOGO-by-guide behaves almost like random K-fold. It is more a stability device than a generalization test.
- About 5 to 10 grouped folds give nearly the same averaged row at 5 to 9x less compute.
- For H1, grouping by **batch** (48 GEM batches) is the more meaningful split.

**6. Strip shared cell-state signal before learning edges.**
- The top "regulators" in the test were histones (cell cycle), MT- genes (cell quality) and VIM/COL1A1 (differentiation). These move together because of shared programs, not regulation.
- Regress out cell-cycle score, mitochondrial fraction, library size and/or the top ~10 PCs before fitting.
- Optionally restrict regulators to transcription factors (about 1,600 genes), as GRNBoost2/SCENIC do.

**7. Denoise sparse genes.**
- Sparse genes had R² ≈ 0.06 at the single-cell level.
- Averaging each cell with its ~20 nearest neighbours (metacells) before regressing usually raises signal a lot.
- Alternatively, weight edges by the target's CV R², so unpredictable targets contribute little.

**8. Smaller code fixes.**
- Standardize regulators before ridge. Right now the penalty hits high-variance genes less, which biases edges toward them.
- Choose highly variable genes inside each fold's training set. Today they are chosen on all cells, which is a mild leak.
- Exclude genes that could be the same transcript as the target (paralogs, or readthrough loci such as HIST1H1B with HIST1H1C) if they look trivially predictive.

**9. Validate before submitting.**
Build the edge matrix from H1 controls, simulate H1's 150 knockdowns with step 2, and score with `cell-eval` against the real H1 knockdowns. Compare against the control-mean baseline and the "mean knockdown effect" baseline. This is the only signal we have on whether the GRN adds anything.

## Suggested order

1. Add the all-genes ridge fast path, the slope stat and stability stats (items 1, 3 and 4).
2. Run it on H1 controls with cell-state signal removed (item 6).
3. Run the propagation (item 2) on the 150 H1 knockdowns and score with cell-eval (item 9).
4. If it beats the baselines, run it on A/B/C and package with `vcc prep`.
