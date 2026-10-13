# From control-cell networks to public screens

*Research summary for the 2026 Arc Virtual Cell Challenge, 13 October 2026. Code: [`models/grn_ridge`](../models/grn_ridge), [`models/logo_shap`](../models/logo_shap), [`models/public_screens`](../models/public_screens). Detailed numbers: [`models/public_screens/RESULTS.md`](../models/public_screens/RESULTS.md). Full improvement audit with citations: [`IMPROVEMENTS_AUDIT.md`](IMPROVEMENTS_AUDIT.md). A styled version of this page is [`index.html`](index.html).*

---

## In one paragraph

> We are asked to predict how three unfamiliar cell types respond when one gene is switched down, given only their normal cells. Our first idea, learning a wiring diagram of genes from normal cells and simulating a knockdown through it, does not work: on last year's stem-cell data its predictions correlate with reality at 0.00. What does work is borrowing from real knockdown experiments. Large public screens in other cell lines tell us which genes behave alike when knocked down; predicting a knockdown as a blend of the effects of genes that behave like it lifts correlation from 0.19 to 0.25, and on the challenge's official scorer it moves us from chance (0.52) to 0.66 at telling knockdowns apart. A shuffled control confirms the gain comes from the screens, and confidence intervals with held-out tuning confirm it is not noise. The next step is to test the same idea across cell lines, which is what the challenge actually scores.

---

## 1. Background

- **Gene expression.** Each cell carries about 20,000 genes. How strongly each one is "on" is measured by counting its RNA molecules.
- **Single-cell RNA-seq.** Counts RNA for thousands of genes in each of thousands of cells. The result is a table: one row per cell, one column per gene, mostly zeros.
- **CRISPRi knockdown.** A guide RNA steers a disabled Cas9 to one gene and turns it down. Other genes then change because genes regulate each other. The full set of changes is the knockdown's **effect**.
- **Non-targeting controls.** Cells given a guide that matches nothing. They show the cell's normal state.
- **Perturb-seq.** Thousands of different knockdowns done at once in one dish, with single-cell RNA-seq reading out both which guide each cell got and its expression [Replogle2022].

## 2. The task

- 2026 Arc Virtual Cell Challenge [Roohani2025]. For each of three new cell lines (contexts A, B, C) we get 18,400 control cells and a list of 300 target genes. We must output 400 predicted cells for each target in each context.
- **Zero-shot:** we never see knockdown data from A, B or C. Arc ran those experiments and scores us against them.
- The contexts differ a lot: A looks lymphoid (immune), B and C epithelial. The final phase uses three more unseen contexts (D, E, F).
- **What we may train on:** last year's H1 stem-cell data (150 knockdowns with real effects) and any public data.
- **H1 is our test bench.** It is the only dataset in hand where we can hide a knockdown, predict it, and check.

## 3. Rounds 1–2: a wiring diagram from control cells does not work

**The idea** (from [`notes/grn_linear_model.pdf`](../notes/grn_linear_model.pdf)): learn which genes predict which across normal cells, giving an N × N "edge matrix", then lower one gene and let the change ripple through the edges.

**What was built.** A cross-validated regression for every gene against all others, with SHAP to rank regulators, then a fast version that fits all 11,000 genes in one matrix solve (about 100,000× faster than a gene-by-gene loop).

**The test.** For each of H1's 147 usable knockdowns, compare the simulated change with the real one.

| Method | Pearson with the real change |
|---|---|
| Network simulation | 0.00 |
| Plain co-expression | 0.00 |
| "Every knockdown does what the average knockdown does" | **0.19** |
| Smoothed (metacell) network | −0.005 |

**Takeaway.** How genes co-vary in normal cells is not the route a knockdown's effect follows. This matches a broad benchmark of network-based forecasting methods [Kernfeld2025] and the finding that simple baselines are hard to beat on perturbation prediction [AhlmannEltze2025]. The useful leftover: an **oracle** test showed that if we knew which other knockdowns a gene resembles, borrowing their effects would reach 0.46. So the job is to learn gene similarity from somewhere else.

## 4. How the challenge scores, and why the average fails it

Arc scores with its open-source `cell-eval` package [cell-eval]. Three numbers matter here.

| Metric | What it asks | Better | Chance / baseline |
|---|---|---|---|
| **overlap_at_N** | Run a differential-expression test (FDR < 0.05) on real and on predicted cells. Of the real top-N changed genes, how many are also in our top N? | Higher | 0.146 for the repo's control-mean baseline |
| **discrimination_score_l1** | Rank all real knockdown effects by L1 distance to our predicted effect (target gene excluded). Score = 1 − rank / number of knockdowns: 1.0 = own knockdown ranked first, 0.5 = random | Higher | ~0.52 |
| **mae** | Average absolute error of predicted vs real mean expression per gene | Lower | 0.049 for unchanged controls |

The leaderboard score rescales each metric against a baseline (0 = baseline, 1 = perfect, clipped at 0) and averages the three [cell-eval]. Arc's 2025 wrap-up noted that almost all entries scored worse than baseline on MAE, and that discrimination rises with the size of predicted effects [ArcWrapUp2025].

**Key point.** The average-knockdown prediction wins on correlation but is identical for every knockdown, so it cannot tell them apart: discrimination stays at chance (0.52). To score, a prediction must be *specific* to each gene.

## 5. Round 3: borrow effects from genes that behave alike in public screens

**The data.** Public Perturb-seq screens: Replogle et al. 2022 genome-wide K562 (a leukaemia line, ~10,000 knockdowns) and RPE1 (retinal epithelium) [Replogle2022]; Nadig et al. 2025 HepG2 (liver) and Jurkat (T cell) [Nadig2025].

**The method, in four steps.**
1. In the screens, describe each gene by its knockdown profile: how every other gene moves when it is knocked down.
2. Two genes are **similar** if their profiles correlate.
3. To predict a held-out H1 knockdown, take the 10 most similar *other* H1 knockdowns and average their real H1 effects, weighted by similarity. Blend 50/50 with the average knockdown.
4. For the challenge scorer, **amplify** the knockdown-specific part: prediction = average + C × (similarity-based − average). C = 8 worked best.

**Why amplify?** Averaging ten neighbours and blending with the mean produces a cautious, washed-out prediction. The scorer rewards a confident, distinct change per knockdown. Amplifying undoes the washing-out. This is a calibration step, akin to temperature scaling in ML [Guo2017].

**Results on H1.** Correlation (147 knockdowns, leave-one-out):

| Method | Pearson |
|---|---|
| Average knockdown | 0.189 |
| K562 similarity | 0.243 |
| All four screens | **0.249** |
| Oracle similarity | 0.461 |

Better on 77 of 133 covered knockdowns (Wilcoxon p = 4 × 10⁻⁴).

Official scorer (40 H1 knockdowns, 400 cells each):

| Method | overlap_at_N ↑ | discrimination ↑ | mae ↓ |
|---|---|---|---|
| Control mean (repo baseline) | 0.146 | 0.520 | 0.140 |
| Average knockdown | 0.044 | 0.521 | 0.049 |
| Average knockdown ×8 | 0.203 | 0.509 | 0.067 |
| Screen similarity ×4 | 0.190 | 0.628 | 0.050 |
| **Screen similarity ×8** | **0.238** | **0.661** | 0.069 |
| Oracle (true effect) | 0.548 | 0.904 | 0.031 |

## 6. The null control: is the gain from the screens or from amplifying?

**The test.** Give each knockdown *another* knockdown's similarity-based specific part, amplified the same 8×. Same magnitude, wrong gene.

| | overlap | discrimination |
|---|---|---|
| Screen similarity ×8 | 0.238 | 0.661 |
| Shuffled ×8 (null) | 0.163 | 0.543 |
| Average ×8 (no gene info) | 0.203 | 0.509 |

**Takeaway.** Amplifying anything raises overlap a little, because larger changes push more genes past the significance test. Only the right gene's screen information raises discrimination. That is the evidence the method is learning something real.

## 7. Noise check: are the gains real?

We re-ran everything with error bars (see "Noise checks" in [`models/public_screens/RESULTS.md`](../models/public_screens/RESULTS.md)). 95% confidence intervals come from bootstrapping over knockdowns [Efron1979].

| Check | Result |
|---|---|
| Pearson gain over the average knockdown (147 knockdowns) | +0.060, CI +0.036 to +0.086 |
| Settings tuned on half the knockdowns, scored on the other half (2,000 splits) | 0.244 vs 0.188; better in every split |
| cell-eval ×8 minus control-mean baseline (40 knockdowns) | overlap +0.092 (CI +0.042 to +0.141); discrimination +0.141 (CI +0.059 to +0.229); mae 0.071 lower |
| cell-eval ×8 minus shuffled null ×8 | overlap +0.075 (CI +0.033 to +0.120); discrimination +0.118 (CI +0.027 to +0.211) |
| Amplification tuned on 20 knockdowns, scored on the other 20 | overlap 0.260 vs 0.147, discrimination 0.648 vs 0.519; beats the baseline in 99.8% of splits |
| ×4 minus control-mean baseline | overlap +0.043 (CI −0.024 to +0.110, **not significant**); discrimination +0.108 (CI +0.039 to +0.184) |
| Three random seeds | ×8 overlap 0.236–0.241, discrimination 0.658–0.661 |
| Chance discrimination with 40 knockdowns | anywhere from about 0.41 to 0.59 |
| How noisy is the truth? Correlation of each knockdown's effect between two halves of its cells | median 0.58; 8 of 40 knockdowns below 0.12, which no method can predict |

**Takeaway.** The ×8 result holds up: it beats both the baseline and the shuffled null beyond noise, and tuning C on held-out knockdowns gives the same answer, so we did not fool ourselves by tuning on the test set [VarmaSimon2006; CawleyTalbot2010]. Two cautions. First, ×4's overlap gain is within noise; only its discrimination gain is solid. Second, with 40 knockdowns, chance discrimination can land anywhere between 0.41 and 0.59, so the null's 0.54 is chance, and differences of a few hundredths between methods mean nothing. The truth itself is noisy: a fifth of the H1 knockdowns barely change anything measurable, which caps any method well below the oracle.

## 8. What changes for 2026

Three facts, checked against the 2026 validation files:

- **Only 13 of the 300 targets were knocked down in H1.** We cannot copy a target's own H1 effect; similarity-based borrowing is the route.
- **272 of 300 targets are in the K562 genome-wide screen, and none in RPE1, HepG2 or Jurkat.** Similarity will come from one screen.
- **55% of those K562 knockdowns had a significant effect** (Replogle's energy test). The rest look weak, so one global amplification will over-predict them.

Nadig et al. found knockdown effects only moderately consistent across four cell lines (in their preprint, average correlation ~0.62; K562 vs RPE1 ~0.40) [Nadig2025], so some signal will be lost in transfer.

And the main untested step: on H1 we borrowed *within* one cell line. In 2026 we borrow *across* cell lines, from H1 and K562 into A, B and C. That transfer has not been validated yet.

## 9. Proposed improvements (details and full citations in [`IMPROVEMENTS_AUDIT.md`](IMPROVEMENTS_AUDIT.md))

| # | Idea | In plain words | Grounded in |
|---|---|---|---|
| 1 | Cross-line validation | Test as the challenge does: hide a whole cell line (HepG2 or Jurkat), predict it from H1 + K562 | Roohani 2025; Adduri 2025; Ahlmann-Eltze 2025 |
| 2 | Nested tuning + error bars (done on H1, keep doing) | Pick settings on some knockdowns, score on others; bootstrap CIs | Varma & Simon 2006; Cawley & Talbot 2010; Efron 1979 |
| 3 | Per-knockdown effect size | Amplify strong knockdowns more, weak ones less (empirical-Bayes shrinkage) | James & Stein 1961; Efron & Morris 1975; Stephens 2017; Zhu 2019; Love 2014; Mejia 2025 |
| 4 | Better gene similarity | Compare genes in a denoised low-rank space; use GO/STRING/GenePT where screens are silent; fuse sources | Ahlmann-Eltze 2025; Littman 2025 (PRESAGE); Roohani 2024 (GEARS); Szklarczyk 2023; Chen & Zou 2025; Wang 2014 |
| 5 | Shared + specific split across lines | Model each effect as a context-wide response plus a gene-specific part; adapt each separately | Daumé 2007; Ben-David 2010; Sun 2016; Lotfollahi 2019 (scGen) |
| 6 | Context rescaling and gating | Don't predict changes for genes silent in the context; scale by expression | Hafemeister & Satija 2019; Svensson 2020 |
| 7 | Realistic cells | Count noise and partial knockdowns so the DE test behaves like real data | Svensson 2020; Love 2014; Horlbeck 2016; Replogle 2022 |
| 8 | Arc's State model | Run the published model zero-shot as a baseline and ensemble member | Adduri 2025 (Cell 2026) |
| 9 | Ensembles / stacking | Combine several predictors with weights fit on held-out knockdowns | Wolpert 1992; Breiman 1996; Lakshminarayanan 2017 |
| 10 | Honest metrics | Always report the null, mean-removed Pearson and systematic-variation checks | Viñas Torné 2025 (Systema); Arc 2025 wrap-up; Marrakchi 2026; Csendes 2025 |

Recommended order: 1 first (2 is now done on H1), then 3 and 6, then 4, 5 and 9, with 7 and 8 in parallel.

## Glossary

| Term | Plain meaning |
|---|---|
| Knockdown / CRISPRi | Switching one gene down on purpose with a guide RNA |
| Control (non-targeting) cell | Same procedure, guide hits nothing |
| Effect | How every gene changes when one gene is knocked down |
| Pseudobulk | Averaging many cells into one profile |
| Perturb-seq / screen | Thousands of knockdowns read out by single-cell RNA-seq in one experiment |
| Context | A cell line (A, B, C in 2026; H1, K562 etc. elsewhere) |
| Zero-shot | Predicting for a context without any knockdown data from it |
| GRN / edge matrix | A table of "gene A influences gene B" |
| Pearson r | Correlation from −1 to 1; 0 = no relationship |
| Average knockdown / mean effect | The average change over many knockdowns; the same prediction for every gene |
| Similarity | How alike two genes' knockdown profiles are in a public screen |
| Amplification C | How many times the gene-specific part of a prediction is scaled up |
| Null control | The same pipeline with the gene information scrambled |
| Oracle | A cheat that uses the true answer, showing the ceiling |
| cell-eval | Arc's official scorer |
| overlap_at_N | Share of truly changed genes we also flag as changed |
| discrimination_score_l1 | How often our prediction is closest to its own knockdown (0.5 = chance) |
| mae | Average absolute error |
| Shrinkage / empirical Bayes | Pulling noisy estimates toward a shared value, more for noisier ones |
| Domain adaptation | Transferring a model from one setting (cell line) to another |

## Questions

- **"Is 0.66 discrimination good?"** It is well above chance (0.52) and the control-mean baseline, but far from the oracle (0.90). It shows the signal is real; most of the gap is still open.
- **"Did you tune on the test set?"** Originally yes. The noise check (section 7) re-tuned C on 20 knockdowns and scored on the other 20: overlap 0.260 and discrimination 0.648 held out, beating the baseline in 99.8% of splits.
- **"Couldn't amplification alone explain it?"** No: the shuffled null uses the same amplification and stays near chance on discrimination.
- **"Why not just use a big foundation model?"** Published benchmarks found scGPT-type models do not yet beat simple linear baselines on this task [AhlmannEltze2025; Csendes2025]. We plan to run Arc's State model [Adduri2025] as a comparison.
- **"Will H1 results hold in A, B, C?"** Unknown. H1 tested borrowing within one cell line. Idea 1 builds a cross-line test using HepG2 and Jurkat before we trust any 2026 number.
- **"Whose work is this built on?"** Replogle et al. and Nadig et al. for the screens; Arc for the challenge, the scorer and State; Ahlmann-Eltze, Huber and Anders for showing simple linear models with screen-derived gene embeddings are strong. The references list credits everyone.

---

## References

The same list as in [`IMPROVEMENTS_AUDIT.md`](IMPROVEMENTS_AUDIT.md). Bracketed keys in the text point here.

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
