# Improvements audit: public-screen effect borrowing

*Written 13 October 2026, with the noise-check results added the same day. Companion to [`research_summary.md`](research_summary.md) and [`models/public_screens/RESULTS.md`](../models/public_screens/RESULTS.md). Every idea below cites the work it comes from; the full references are at the end, keyed by [AuthorYear].*

## Where we stand

The current best method predicts a knockdown's change as a weighted average of other knockdowns' real H1 changes, weighted by how alike the two target genes look in public Perturb-seq screens [Replogle2022; Nadig2025]. The knockdown-specific part is then amplified C times. On Arc's `cell-eval` [cell-eval] over 40 H1 knockdowns, C = 8 gives overlap_at_N 0.238 and discrimination 0.661, against 0.146 / 0.520 for the control-mean baseline and 0.548 / 0.904 for the oracle. A shuffled null scores 0.163 / 0.543, so the discrimination gain comes from the screens. The noise check (bootstrap CIs and held-out tuning, in [`models/public_screens/RESULTS.md`](../models/public_screens/RESULTS.md)) confirms both gains are beyond noise for ×8.

Three facts checked for this audit change what matters most next (computed 2026-10-13 from the 2026 validation files and the public screen data):

| Fact | Number | Consequence |
|---|---|---|
| 2026 target genes also knocked down in H1 | **13 of 300** | We cannot borrow "the same gene's" H1 effect. Borrowing by similarity, as now, is the only route from H1. |
| 2026 targets in the Replogle K562 genome-wide screen | **272 of 300** | Similarity is available for almost every target, but only from one screen. RPE1, HepG2, Jurkat and the K562 essential screen cover **0** of them. |
| Those 272 K562 knockdowns with a significant transcriptional effect (Replogle's energy test, p < 0.05) | **55%** (48% screen-wide); 96% knocked the target below half | Nearly half the targets have weak or no measurable effect in K562. A single global amplification C will over-predict them. |
| 2026 targets detected in ≥ 5% of control cells | 299 / 297 / 298 in A / B / C | Every target is expressed in each context, so knockdowns are not trivially null. |

The H1 test so far borrowed effects **within the same cell line** (other H1 knockdowns, scored on H1). The 2026 task borrows **across cell lines** (from H1 and K562 into A, B and C). The biggest risk is therefore not the similarity step but the transfer step, which no experiment has tested yet. The audit is ordered accordingly.

## Summary table

Effort: S = hours, M = a day or two, L = a week or more. "Disc." = discrimination_score_l1, "Overlap" = overlap_at_N.

| # | Improvement | Main metric it targets | Effort | Grounding |
|---|---|---|---|---|
| 1 | A validation set that mimics 2026: new cell line, unseen targets | All (trustworthy numbers) | M | [Roohani2025; Adduri2025; AhlmannEltze2025; Nadig2025; BenDavid2010] |
| 2 | Nested tuning + bootstrap CIs (done on H1; keep as standard) | All | S | [VarmaSimon2006; CawleyTalbot2010; Efron1979] |
| 3 | Per-knockdown amplitude from screen strength (empirical-Bayes shrinkage instead of one global C) | Overlap, MAE | S–M | [JamesStein1961; EfronMorris1975; Stephens2017; Zhu2019; Love2014; Guo2017; Mejia2025] |
| 4 | Better gene similarity: low-rank embeddings of the screen, fused with prior knowledge for uncovered genes | Disc. | M | [AhlmannEltze2025; Littman2025; Roohani2024; Szklarczyk2023; ChenZou2025; Wang2014; CandesRecht2009] |
| 5 | Shared-plus-specific decomposition across cell lines; weight sources by context similarity | Disc., MAE | M | [Daume2007; BenDavid2010; SunFengSaenko2016; Lotfollahi2019] |
| 6 | Context-aware rescaling and expression gating of predicted changes | MAE, Overlap | S | [Hafemeister2019; Svensson2020; AhlmannEltze2025] |
| 7 | Realistic single-cell generation (count noise, incomplete knockdown) | Overlap, Disc. | M | [Svensson2020; Love2014; Replogle2022; Horlbeck2016] |
| 8 | Run Arc's State model zero-shot as a baseline and ensemble member | All | M–L | [Adduri2025; Roohani2025] |
| 9 | Ensemble / stack several predictors | All | S–M | [Wolpert1992; Breiman1996; Lakshminarayanan2017] |
| 10 | Report perturbation-specific metrics and guard against metric gaming | Trust | S | [VinasTorne2025; ArcWrapUp2025; Marrakchi2026; Csendes2025] |

Recommended order: 1 first (2 is done on H1 and should be reused there), then 3 and 6 (cheap, likely wins), then 4, 5 and 9, with 7 and 8 in parallel on Colab.

---

## 1. Build a validation set that looks like 2026

**What changes.** Today the method is tuned and scored on H1 knockdowns, borrowing from other H1 knockdowns. Add a second benchmark that matches the challenge: hold out a whole cell line, give the model only its control cells, and ask it to predict knockdowns of genes it has not seen in that line. Two ready-made options:

- **HepG2 and Jurkat from the 2025 support set** [Nadig2025] already have single cells for controls and 68 knockdowns each. Predict them from H1 + K562 only, and score with `cell-eval`. Jurkat is a T-cell line, close to context A (lymphoid-like); HepG2 is epithelial, closer to B and C.
- **K562 → RPE1**, the split used by Ahlmann-Eltze et al. [AhlmannEltze2025] and by Arc's State paper [Adduri2025] for cross-context generalization.

**Why it should help.** Nadig et al. report only moderate consistency of knockdown effects across K562, RPE1, HepG2 and Jurkat (in their preprint, average pairwise correlation ~0.62, from ~0.40 for K562/RPE1 to ~0.75 for HepG2/RPE1) [Nadig2025]. Cross-line transfer will therefore lose signal, and we need to measure how much. Ben-David et al.'s domain-adaptation bound [BenDavid2010] says error on a new domain depends on how different that domain is from the training one, which a within-H1 test never measures. Arc framed the challenge itself around exactly this generalization gap [Roohani2025]. Tuning C, m and the shrinkage on a cross-line split also stops us picking settings that only suit H1 [VarmaSimon2006].

## 2. Tune without peeking, and put error bars on every number

**Status: done for H1 (2026-10-13; see "Noise checks" in [`models/public_screens/RESULTS.md`](../models/public_screens/RESULTS.md)).** Held-out tuning and bootstrap CIs over knockdowns [Efron1979] confirm the ×8 gains: overlap +0.092 (CI +0.042 to +0.141) and discrimination +0.141 (CI +0.059 to +0.229) over the control-mean baseline, and C tuned on 20 knockdowns scores 0.260 / 0.648 on the other 20. ×4's overlap gain is not significant (CI −0.024 to +0.110).

**What changes now.** Keep the same protocol for every new idea, and carry it to the cross-line benchmark in idea 1. Two lessons from the check shape that protocol:
- With 40 knockdowns, chance discrimination ranges from about 0.41 to 0.59. Differences of a few hundredths between methods are noise; compare methods on paired bootstraps, not on point estimates.
- The H1 truth is itself noisy (split-half correlation median 0.58; 8 of 40 knockdowns below 0.12). Score on more knockdowns where possible (the 147 with pseudobulk effects for Pearson), and report results separately for reliable and unreliable knockdowns so a method isn't judged on knockdowns nobody can predict.

**Why it matters.** Selecting hyperparameters on the same data used to report performance gives optimistic numbers [VarmaSimon2006], and the bias grows with the number of settings tried [CawleyTalbot2010].

## 3. Give each knockdown its own effect size

**What changes.** Replace `mean + C × (similarity-weighted − mean)` with a per-knockdown multiplier C_k that grows with the evidence for a real effect. Candidate signals, all already in the K562 file's metadata: Replogle's energy-test p-value, the number of Anderson–Darling DE genes, and how strongly the target was knocked down (`fold_expr`) [Replogle2022]. Fit the mapping from these signals to C_k on held-out knockdowns. Two principled versions:

- **Empirical-Bayes shrinkage.** Treat each knockdown's specific effect as a noisy estimate and shrink it toward the shared response in proportion to its noise. This is the James–Stein / Efron–Morris idea [JamesStein1961; EfronMorris1975], used in genomics as fold-change shrinkage in DESeq2 and apeglm [Love2014; Zhu2019] and adaptive shrinkage (ashr) [Stephens2017]. Heavy-tailed priors (apeglm) are the right fit here: most knockdowns small, a few large.
- **Post-hoc calibration.** Fit one or two scalars on held-out data, as temperature scaling does for classifiers [Guo2017]. Our global C is already a crude version of this; making it depend on screen strength is the next step.

**Why it should help.** About 45% of the 2026 targets show no significant effect in K562, and in H1 a fifth of knockdowns (8 of 40) have effects too weak to reproduce between halves of their own cells. Amplifying their predicted change 8× invents DE genes, which costs MAE and adds false positives to overlap. Strong knockdowns, meanwhile, are probably under-amplified. Averaging ten neighbours and blending 50/50 with the mean effect shrinks every prediction toward the mean; amplification undoes that shrinkage, but uniformly. Per-knockdown shrinkage is the textbook fix. Arc's 2025 wrap-up notes that the discrimination score rises with predicted effect magnitude, and that top teams traded MAE for discrimination and DE overlap [ArcWrapUp2025]; Mejia et al. describe the opposite failure, "mode collapse" toward the mean, which scale-aware metrics penalize [Mejia2025]. A per-knockdown scale sits between the two: large where the evidence is strong, small where it is not.

## 4. Learn a better gene similarity

**What changes.** Similarity is now the plain correlation of two genes' z-scored K562 profiles over ~8,000 readout genes. Options, in order of effort:

1. **Low-rank embedding first.** Compute the top ~50 principal components of the screen's knockdown × gene matrix and measure similarity in that space. Ahlmann-Eltze et al. [AhlmannEltze2025] showed a linear model whose gene embeddings come from a PCA of *perturbation data from another cell line* matched or beat GEARS and scGPT, and that pre-training on K562 perturbations helped predict RPE1. Denoising by low rank is the matrix-completion idea [CandesRecht2009]. PRESAGE likewise found the knowledge source matters more than the architecture, with Perturb-seq data from other cell systems the most useful source [Littman2025].
2. **Prior-knowledge similarity for genes the screen lacks** (28 targets). GEARS builds a gene graph from Gene Ontology co-annotation [Roohani2024]; STRING gives protein-association scores [Szklarczyk2023]; GenePT embeds gene descriptions with a language model [ChenZou2025]. Any of these can stand in where K562 has no profile, or give a second opinion where the K562 profile is weak.
3. **Fuse the similarities.** Similarity network fusion [Wang2014] combines several similarity matrices into one; a simpler learned weighted average, fit on validation, is a good first step.
4. **Learn the similarity end to end** (metric learning): fit a small model that maps two genes' screen profiles to the correlation of their true effects, training on pairs where both are known. The oracle similarity reaches Pearson 0.46 against our 0.25, so the gap is large.

**Why it should help.** Discrimination rewards predictions that are distinct for each knockdown and close to the right one. Better neighbours give more distinct and more correct knockdown-specific parts. Raw correlation over 8,000 noisy readouts is dominated by noise for weak knockdowns; a low-rank space keeps the shared biological programs and drops the noise.

## 5. Separate the shared response from the gene-specific one, across cell lines

**What changes.** Write each knockdown's effect in each cell line as `shared response (cell line) + gene-specific effect (gene) + interaction`. Estimate the shared part per source line (H1, K562, HepG2, Jurkat) and predict it for A/B/C from how similar their control profiles are to each source line. Borrow the gene-specific part through the similarity step. This is Daumé's "frustratingly easy" domain adaptation (shared + domain-specific copies of each feature) [Daume2007]. Weighting sources by context similarity follows the domain-adaptation bound [BenDavid2010]; aligning their covariances before transfer is CORAL [SunFengSaenko2016]. scGen's latent-space arithmetic, which moves a perturbation effect learned in one cell type onto another [Lotfollahi2019], is the deep-learning version of the same idea.

**Why it should help.** The mean effect is most of the pseudobulk Pearson (0.19 of 0.25) but scores at chance on discrimination, because it is identical for every knockdown. Getting the shared part right per context mainly helps MAE; getting the specific part right helps discrimination. Treating them separately lets each be tuned for its own metric, and avoids carrying H1's stem-cell-specific shared response into lymphoid or epithelial cells.

## 6. Rescale effects to each context and gate by expression

**What changes.** Before applying a borrowed change to context A, B or C:

- **Gate:** set the predicted change to zero for genes barely detected in that context's controls (e.g. < 1% of cells). A gene that is off cannot go down, and a large log change on it is noise.
- **Rescale:** scale the change by how expressed the gene is in the target context relative to the source, since log-scale effects on low-count genes are both smaller and noisier. The negative-binomial mean–variance relationship behind sctransform [Hafemeister2019; Svensson2020] gives the scale. Ahlmann-Eltze et al. found K562→RPE1 predictions were more accurate for genes whose expression was similar in the two lines [AhlmannEltze2025], which suggests down-weighting genes whose expression differs most between source and target.

**Why it should help.** A/B/C have ~20,000 UMIs per cell against H1's ~55,000, and their expressed genes differ a lot (only 541 of the top 2,000 variable genes are shared by all three). Changes predicted on genes that are silent in a context cost MAE and can create false DE calls in the overlap metric.

## 7. Generate cells the way real cells look

**What changes.** Submissions currently resample 400 control cells per knockdown and shift them by the predicted change. Two refinements:

- **Count noise.** Generate counts from a negative-binomial with the predicted mean per cell and the gene's dispersion estimated from controls [Love2014; Svensson2020], instead of shifting log values.
- **Incomplete knockdown.** Real CRISPRi populations contain cells where the guide worked weakly or not at all; knockdown strength varies by guide [Horlbeck2016], and Replogle et al. filtered and measured this per perturbation [Replogle2022]. Model each knockdown's cells as a mixture: a fraction fully affected, the rest near control.

**Why it should help.** overlap_at_N comes from a differential-expression test run on our cells against controls, keeping genes at false discovery rate < 0.05 [BenjaminiHochberg1995], ranking them by |log2 fold change| and measuring recall of the true top-N genes [cell-eval]. How many genes pass depends on the cell-to-cell spread, not only the mean change. Matching real spread makes our DE calls behave like the real ones, so a correct mean change turns into a correct DE set.

## 8. Run Arc's State model zero-shot

**What changes.** State [Adduri2025] is Arc's model trained to predict perturbation responses across contexts from control cells, and the reference model for the challenge [Roohani2025]. Its paper reports over 30% better discrimination than baselines and better DE-gene recovery. Run the released weights on the 2026 controls for the 300 targets, score it on the validation set from idea 1, and use it as both a baseline and an ensemble member.

**Why it should help.** It is the closest published model to the exact task. If it beats our method on the cross-line validation, we should build on it; if not, it still adds diversity to an ensemble. Note the caution from benchmarks that large models often fail to beat simple baselines on perturbation prediction [AhlmannEltze2025; Csendes2025; Wu2025]; the validation set decides.

## 9. Ensemble and stack predictors

**What changes.** Combine several predictions of each knockdown's change (similarity borrowing with different neighbour sets or similarity sources, the K562 profile mapped directly, State, a GEARS-style model) with weights fit on validation knockdowns (stacking) [Wolpert1992]. A cheap version is bagging over random subsets of neighbours [Breiman1996].

**Why it should help.** Averaging diverse predictors reduces variance and tends to improve calibration [Breiman1996; Lakshminarayanan2017]. Ours are noisy for weak knockdowns, where averaging helps most. Stacking also lets the weights differ by metric: favour the most distinct predictor for discrimination and the most accurate for MAE.

## 10. Measure what is specific, and don't game the metric

**What changes.** Alongside `cell-eval`, always report (a) the shuffled-null control (already done), (b) Pearson on the change *after* subtracting the mean knockdown effect, and (c) the systematic-variation check from Systema [VinasTorne2025].

**Why it should help.** Several benchmarks found that much of what perturbation models score comes from systematic, perturbation-independent shifts that a mean baseline captures for free [VinasTorne2025; AhlmannEltze2025; Csendes2025]. Amplification raises overlap even with no knowledge of the gene (our "mean effect ×8" row gets 0.203 overlap with chance discrimination). Arc's own 2025 wrap-up reports that almost all entries scored worse than baseline on MAE and that discrimination rises with effect magnitude [ArcWrapUp2025]; Marrakchi et al. argue raw-pseudobulk discrimination can misrank perturbations [Marrakchi2026]. Keeping metrics that only reward gene-specific signal stops us from optimizing towards a high but empty score.

---

## Ideas considered and not recommended now

- **Polishing the control-cell network** (more layers, priors, sparse models). The H1 tests found no causal signal in control-cell co-variation (Pearson ~0.00), consistent with a comparison across 11 perturbation datasets that found expression-forecasting methods rarely beat simple baselines [Kernfeld2025].
- **Training a large foundation model ourselves** (scGPT [Cui2024], Geneformer [Theodoris2023]). Their published perturbation results do not beat linear baselines [AhlmannEltze2025; Csendes2025], and the compute would compete with ideas 1–9. Their gene embeddings remain useful as similarity inputs (idea 4).
- **CPA / CellOT-style distribution models** [Lotfollahi2023; Bunne2023] need perturbed cells from the target context for training, which the zero-shot setting does not provide. Their idea of mapping control cells to perturbed cells per cell (rather than per mean) is folded into idea 7.

---

## References

Bibliographic details were checked against Crossref, Europe PMC, publisher pages or proceedings sites on 2026-10-13. "Online" gives the first online date where it differs from the issue year.

### Perturbation data and the challenge
- **[Replogle2022]** Replogle JM, Saunders RA, Pogson AN, et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell* 185(14):2559–2575.e28 (2022). https://doi.org/10.1016/j.cell.2022.05.013 — the K562 genome-wide and RPE1 screens we use for gene similarity.
- **[Nadig2025]** Nadig A, Replogle JM, Pogson AN, et al. Transcriptome-wide analysis of differential expression in perturbation atlases. *Nature Genetics* 57(5):1228–1237 (2025). https://doi.org/10.1038/s41588-025-02169-3. Preprint: Transcriptome-wide characterization of genetic perturbations, bioRxiv (2024), https://doi.org/10.1101/2024.07.03.601903 — the HepG2 and Jurkat screens; cross-cell-line consistency figures quoted above come from the preprint.
- **[Roohani2025]** Roohani YH, Hua TJ, Tung PY, et al. Virtual Cell Challenge: Toward a Turing test for the virtual cell. *Cell* 188(13):3370–3374 (2025). https://doi.org/10.1016/j.cell.2025.06.008
- **[Adduri2025]** Adduri AK, Gautam D, Bevilacqua B, et al. Predicting cellular responses to perturbation across diverse contexts with State. *Cell* 189(19):5914–5931.e20 (2026). https://doi.org/10.1016/j.cell.2026.07.052. Preprint: bioRxiv (2025), https://doi.org/10.1101/2025.06.26.661135 — the State model; also introduces cell-eval.
- **[cell-eval]** Arc Institute. cell-eval: evaluation suite for perturbation prediction. https://github.com/ArcInstitute/cell-eval (checked at commit 6928cf8). Cite via [Adduri2025], as the repository asks.
- **[ArcWrapUp2025]** Arc Institute. Virtual Cell Challenge 2025 Wrap-Up: Winners and Reflections (6 Dec 2025). https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up

### Perturbation-prediction models and benchmarks
- **[AhlmannEltze2025]** Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nature Methods* 22(8):1657–1661 (2025). https://doi.org/10.1038/s41592-025-02772-6. Preprint: https://doi.org/10.1101/2024.09.16.613342
- **[Roohani2024]** Roohani Y, Huang K, Leskovec J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nature Biotechnology* 42(6):927–935 (2024; online 2023). https://doi.org/10.1038/s41587-023-01905-6
- **[Cui2024]** Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods* 21(8):1470–1480 (2024). https://doi.org/10.1038/s41592-024-02201-0
- **[Theodoris2023]** Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature* 618(7965):616–624 (2023). https://doi.org/10.1038/s41586-023-06139-9 — Geneformer.
- **[Lotfollahi2019]** Lotfollahi M, Wolf FA, Theis FJ. scGen predicts single-cell perturbation responses. *Nature Methods* 16(8):715–721 (2019). https://doi.org/10.1038/s41592-019-0494-8
- **[Lotfollahi2023]** Lotfollahi M, Klimovskaia Susmelj A, De Donno C, et al. Predicting cellular responses to complex perturbations in high-throughput screens. *Molecular Systems Biology* 19(6):e11517 (2023). https://doi.org/10.15252/msb.202211517 — CPA.
- **[Bunne2023]** Bunne C, Stark SG, Gut G, et al. Learning single-cell perturbation responses using neural optimal transport. *Nature Methods* 20(11):1759–1768 (2023). https://doi.org/10.1038/s41592-023-01969-x — CellOT.
- **[Littman2025]** Littman R, Levine J, Maleki S, et al. Gene-embedding-based prediction and functional evaluation of perturbation expression responses with PRESAGE. bioRxiv (2025). https://doi.org/10.1101/2025.06.03.657653
- **[Kernfeld2025]** Kernfeld E, Yang Y, Weinstock JS, Battle A, Cahan P. A comparison of computational methods for expression forecasting. *Genome Biology* 26:388 (2025). https://doi.org/10.1186/s13059-025-03840-y
- **[Csendes2025]** Csendes G, Sanz G, Szalay KZ, Szalai B. Benchmarking foundation cell models for post-perturbation RNA-seq prediction. *BMC Genomics* 26:393 (2025). https://doi.org/10.1186/s12864-025-11600-2
- **[Wu2025]** Wu Y, Wershof E, Schmon SM, et al. PerturBench: Benchmarking machine learning models for cellular perturbation analysis. *NeurIPS 38* (2025), pp. 106937–106977. https://arxiv.org/abs/2408.10609
- **[VinasTorne2025]** Viñas Torné R, Wiatrak M, Piran Z, et al. Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation. *Nature Biotechnology* 44(6):1050–1059 (2026; online 2025). https://doi.org/10.1038/s41587-025-02777-8
- **[Mejia2025]** Mejia GM, Miller HE, Leblanc FJA, et al. Diversity by design: addressing mode collapse improves scRNA-seq perturbation modeling on well-calibrated metrics. arXiv:2506.22641 (2025). https://arxiv.org/abs/2506.22641
- **[Marrakchi2026]** Marrakchi Y, D'Ascenzo D, Cultrera di Montesano S. Score distributions, not cells: evaluating single-cell perturbations under class overlap. arXiv:2607.04595 (2026). https://arxiv.org/abs/2607.04595

### Gene knowledge and embeddings
- **[Szklarczyk2023]** Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Research* 51(D1):D638–D646 (2023). https://doi.org/10.1093/nar/gkac1000
- **[ChenZou2025]** Chen Y, Zou J. Simple and effective embedding model for single-cell biology built from ChatGPT. *Nature Biomedical Engineering* 9(4):483–493 (2025). https://doi.org/10.1038/s41551-024-01284-6. Preprint (GenePT): https://doi.org/10.1101/2023.10.16.562533
- **[Horlbeck2016]** Horlbeck MA, Gilbert LA, Villalta JE, et al. Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation. *eLife* 5:e19760 (2016). https://doi.org/10.7554/eLife.19760

### Single-cell statistics
- **[Love2014]** Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15(12):550 (2014). https://doi.org/10.1186/s13059-014-0550-8
- **[Zhu2019]** Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35(12):2084–2092 (2019; online 2018). https://doi.org/10.1093/bioinformatics/bty895 — apeglm.
- **[Stephens2017]** Stephens M. False discovery rates: a new deal. *Biostatistics* 18(2):275–294 (2017). https://doi.org/10.1093/biostatistics/kxw041 — ashr.
- **[Hafemeister2019]** Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology* 20(1):296 (2019). https://doi.org/10.1186/s13059-019-1874-1 — sctransform.
- **[Svensson2020]** Svensson V. Droplet scRNA-seq is not zero-inflated. *Nature Biotechnology* 38(2):147–150 (2020). https://doi.org/10.1038/s41587-019-0379-5
- **[Wang2014]** Wang B, Mezlini AM, Demir F, et al. Similarity network fusion for aggregating data types on a genomic scale. *Nature Methods* 11(3):333–337 (2014). https://doi.org/10.1038/nmeth.2810

### Statistics and machine learning
- **[JamesStein1961]** James W, Stein C. Estimation with quadratic loss. *Proc. Fourth Berkeley Symposium on Mathematical Statistics and Probability* 1:361–379 (1961). https://projecteuclid.org/euclid.bsmsp/1200512173
- **[EfronMorris1975]** Efron B, Morris C. Data analysis using Stein's estimator and its generalizations. *JASA* 70(350):311–319 (1975). https://doi.org/10.1080/01621459.1975.10479864
- **[Efron1979]** Efron B. Bootstrap methods: another look at the jackknife. *Annals of Statistics* 7(1):1–26 (1979). https://doi.org/10.1214/aos/1176344552
- **[BenjaminiHochberg1995]** Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS-B* 57(1):289–300 (1995). https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- **[VarmaSimon2006]** Varma S, Simon R. Bias in error estimation when using cross-validation for model selection. *BMC Bioinformatics* 7:91 (2006). https://doi.org/10.1186/1471-2105-7-91
- **[CawleyTalbot2010]** Cawley GC, Talbot NLC. On over-fitting in model selection and subsequent selection bias in performance evaluation. *JMLR* 11:2079–2107 (2010). https://www.jmlr.org/papers/v11/cawley10a.html
- **[Guo2017]** Guo C, Pleiss G, Sun Y, Weinberger KQ. On calibration of modern neural networks. *ICML*, PMLR 70:1321–1330 (2017). https://proceedings.mlr.press/v70/guo17a.html
- **[CandesRecht2009]** Candès EJ, Recht B. Exact matrix completion via convex optimization. *Foundations of Computational Mathematics* 9(6):717–772 (2009). https://doi.org/10.1007/s10208-009-9045-5
- **[Wolpert1992]** Wolpert DH. Stacked generalization. *Neural Networks* 5(2):241–259 (1992). https://doi.org/10.1016/S0893-6080(05)80023-1
- **[Breiman1996]** Breiman L. Bagging predictors. *Machine Learning* 24(2):123–140 (1996). https://doi.org/10.1007/BF00058655
- **[Lakshminarayanan2017]** Lakshminarayanan B, Pritzel A, Blundell C. Simple and scalable predictive uncertainty estimation using deep ensembles. *NeurIPS 30* (2017). https://papers.nips.cc/paper_files/paper/2017/hash/9ef2ed4b7fd2c810847ffa5fa85bce38-Abstract.html
- **[Daume2007]** Daumé III H. Frustratingly easy domain adaptation. *Proc. 45th ACL*, pp. 256–263 (2007). https://aclanthology.org/P07-1033/
- **[BenDavid2010]** Ben-David S, Blitzer J, Crammer K, et al. A theory of learning from different domains. *Machine Learning* 79(1–2):151–175 (2010). https://doi.org/10.1007/s10994-009-5152-4
- **[SunFengSaenko2016]** Sun B, Feng J, Saenko K. Return of frustratingly easy domain adaptation. *Proc. AAAI* 30(1) (2016). https://doi.org/10.1609/aaai.v30i1.10306 — CORAL.
