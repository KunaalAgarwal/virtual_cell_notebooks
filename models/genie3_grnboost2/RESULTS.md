# GENIE3 and GRNBoost2: do they help predict knockdowns? (2026-10-18)

## Bottom line

No, not as a predictor. Both were run with arboreto's own settings and tested on the 147 real H1 knockdowns, the same leave-one-out benchmark used in `logo-shap/v2/RESULTS.md` and `public-screens/RESULTS.md`.

1. **Pushing a knockdown through the network gives Pearson about 0.00**, for GRNBoost2 and GENIE3 alike. That is the same result as our ridge and LightGBM networks.
2. **Using the network to pick which known knockdowns to borrow ties the plain average** (0.193 vs 0.189, difference inside noise). Public-screen similarity gets 0.249.
3. **Mixing network similarity into screen similarity does not help.** GRNBoost2 adds +0.003 (95% CI −0.007 to +0.014). GENIE3 and the K562-trained network make it slightly worse.
4. **GRNBoost2 trained on the K562 screen itself** (perturbation profiles as rows instead of control cells) is no better. It gives 0.179 for borrowing and −0.006 for propagation, while K562 profiles used directly give 0.196–0.249.

The literature says the same thing, and it explains the failure. These are co-expression methods. The strongest "regulators" they find in control cells are cell-cycle genes (UBE2C, CDK1, PLK1, histones) and hypoxia/glycolysis genes (SLC16A3, EGLN1, PFKFB4, PDK1, DDIT4): shared programs, not the causal paths a knockdown follows.

**Where they could still help:** describing a context's own programs (for example, which genes move together in context A vs B), or as a feature in a model trained on perturbation data. They are not a source of knockdown effects, and there is nothing to submit from them. cell-eval was not run, because nothing beat the baseline on Pearson.

## What was run

| Network | Trained on | Genes | Settings | Time (4 CPUs) |
|---|---|---|---|---|
| GRNBoost2 | 6,000 H1 control cells from all 48 batches (2025 support set) | 2,000 HVGs plus all 147 knockdown targets = 2,145 | arboreto defaults: gradient boosting, lr 0.01, up to 5,000 trees, max_features 0.1, subsample 0.9, early stopping window 25 | 38 min |
| GENIE3 | 3,000 of those cells | same 2,145 | random forest, max_features sqrt, **100 trees** (default is 1,000; cut for time) | 71 min |
| GRNBoost2 on K562 | 11,258 Replogle K562 genome-wide perturbation profiles (z-scored pseudobulk) | the 2,000 most variable across perturbations plus 120 H1 targets = 2,085; 720 fitted as targets | arboreto defaults | 55 min |

Every gene could be a regulator; there was no TF list, so non-TF knockdown targets still get out-edges. Each regressor is arboreto's `core.infer_partial_network`, the per-target routine that `grnboost2()`/`genie3()` call, driven through joblib because arboreto 0.1.6 does not run under current dask. Edges get their sign from the Spearman correlation, the same rule pySCENIC uses.

## Results (147 H1 knockdowns, leave-one-out)

Pearson is computed between predicted and real change for each knockdown, then averaged over knockdowns; the target gene itself is excluded. Full tables are in `results/summary_{gb,g3,k562}.csv` (per-knockdown scores and the network .npz files are in the project folder `genie3-grnboost2/`), and the bootstraps in `results/paired_checks.txt`.

**A. Propagate the knockdown through the network** (scored on each network's own genes)

| Method | GRNBoost2 (H1) | GENIE3 (H1) | GRNBoost2 (K562 rows) |
|---|---|---|---|
| Mean effect of the other knockdowns | 0.184 | 0.184 | 0.204 |
| Network, direct targets | −0.007 | −0.003 | −0.006 |
| Network, 3 steps | −0.008 | −0.001 | −0.006 |
| Co-expression (Spearman) | −0.001 | −0.000 | 0.004 |
| Mean effect + fitted network term | 0.184 | 0.184 | 0.203 |

**B. Borrow other knockdowns' real effects, weighted by similarity** (scored on 10,993 genes; best setting over m = 5–40 neighbours and 0/50% shrinkage)

| Similarity | Pearson | vs uniform (95% CI) |
|---|---|---|
| Uniform (mean effect) | 0.189 | |
| GRNBoost2 H1, out-edge profiles | 0.193 | +0.004 [−0.006, +0.015] |
| GENIE3 H1, out-edge profiles | 0.193 | +0.004 [−0.005, +0.014] |
| GRNBoost2 K562, out-edge profiles | 0.179 | −0.010 [−0.024, +0.005] |
| In-edge profiles (any network) | 0.187–0.191 | |
| **Public screens (K562+RPE1+HepG2+Jurkat), our current method** | **0.249** | +0.060 (see public-screens/RESULTS.md) |
| Oracle (true effect similarity) | 0.461 | |

**C. Average screen similarity with network similarity**

| Added to screen similarity | Pearson | vs screens alone (95% CI) | Better on held-out half (2,000 tune/score splits) |
|---|---|---|---|
| GRNBoost2 H1 out-edges | 0.252 | +0.003 [−0.007, +0.014] | 70% (mean +0.003) |
| GENIE3 H1 out-edges | 0.238 | −0.011 [−0.025, +0.004] | 8% |
| GRNBoost2 K562 out-edges | 0.239 | −0.010 [−0.023, +0.003] | 14% |
| any in-edges | 0.214–0.227 | −0.02 to −0.035, CI excludes 0 | 0–2% |

The +0.003 from GRNBoost2 is the best of about 16 settings tried, and its CI includes 0, so treat it as no effect.

Other observations:
- The GRNBoost2 network is sparse: early stopping leaves 4.9% of H1 edges and 2% of K562 edges non-zero. GENIE3's are all non-zero.
- The knocked-down genes rank as middling hubs (median rank 333 of 2,145 by total out-importance).

## Why this matches the literature

- **Expression-only GRN methods barely beat random at finding real regulatory edges.** In BEELINE [Pratapa2020], GENIE3 and GRNBoost2 were among the better methods. Against cell-type ChIP-seq, though, most methods' early precision was close to random (ratio ~1, never above 1.5). [ChenMar2018] and [McCalla2023] report similar results.
- **On real Perturb-seq, GRNBoost's edges are mostly not the large knockdown effects.** CausalBench [Chevalley2025] scores network inference against Replogle K562/RPE1 knockdowns [Replogle2022]. GRNBoost had high recall but low precision, and it was the worst method on the mean Wasserstein distance of its edges. That metric measures whether predicted edges are the ones a knockdown actually moves. Methods that also use the interventional data did not beat observational ones.
- **Forecasting tools built on GRNs rarely beat the mean baseline.** Across 11 perturbation datasets, "the mean or median baseline was almost always the top performer" [Kernfeld2025]. Deep models do not beat simple linear baselines either [AhlmannEltze2025].
- **SCENIC prunes GRNBoost2's raw output with motif enrichment** [Aibar2017; VandeSande2020], because the authors treat many raw edges as indirect. CellOracle [Kamimoto2023] simulates knockouts on a motif-constrained, per-cluster network and frames the output as a direction of cell-state change, not expression values. Both would need TF motif and chromatin priors for each 2026 context, and we don't have those.

## References

- [HuynhThu2010] Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. Inferring regulatory networks from expression data using tree-based methods. PLoS ONE 5:e12776 (2010). doi:10.1371/journal.pone.0012776
- [Moerman2019] Moerman T, Aibar S, Bravo González-Blas C, et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks. Bioinformatics 35:2159–2161 (2019). doi:10.1093/bioinformatics/bty916
- [Aibar2017] Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods 14:1083–1086 (2017). doi:10.1038/nmeth.4463
- [VandeSande2020] Van de Sande B, et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis. Nat Protoc 15:2247–2276 (2020). doi:10.1038/s41596-020-0336-2
- [Pratapa2020] Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. Nat Methods 17:147–154 (2020). doi:10.1038/s41592-019-0690-6
- [ChenMar2018] Chen S, Mar JC. Evaluating methods of inferring gene regulatory networks highlights their lack of performance for single cell gene expression data. BMC Bioinformatics 19:232 (2018). doi:10.1186/s12859-018-2217-z
- [McCalla2023] McCalla SG, et al. Identifying strengths and weaknesses of methods for computational network inference from single-cell RNA-seq data. G3 13:jkad004 (2023). doi:10.1093/g3journal/jkad004
- [Replogle2022] Replogle JM, et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. Cell 185:2559–2575 (2022). doi:10.1016/j.cell.2022.05.013
- [Chevalley2025] Chevalley M, Roohani YH, Mehrjou A, Leskovec J, Schwab P. A large-scale benchmark for network inference from single-cell perturbation data (CausalBench). Commun Biol 8:412 (2025). doi:10.1038/s42003-025-07764-y
- [Kernfeld2025] Kernfeld E, Yang Y, Weinstock JS, Battle A, Cahan P. A comparison of computational methods for expression forecasting. Genome Biol 26:388 (2025). doi:10.1186/s13059-025-03840-y
- [AhlmannEltze2025] Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. Nat Methods 22:1657–1661 (2025). doi:10.1038/s41592-025-02772-6
- [Kamimoto2023] Kamimoto K, et al. Dissecting cell identity via network inference and in silico gene perturbation. Nature 614:742–751 (2023). doi:10.1038/s41586-022-05688-9

The CausalBench precision and recall figures come from the arXiv v2 tables (arXiv:2210.17283).

## Files and how to rerun

```
python prep_h1_controls.py competition_train.h5 h1_pseudobulk.npz grn_h1_pc0.top_edges.csv h1_ctrl_2k.npz 6000 2000
python run_grn.py h1_ctrl_2k.npz grnboost2 h1_grnboost2.npz 6000 0 4
python run_grn.py h1_ctrl_2k.npz genie3 h1_genie3.npz 3000 100 4
python prep_k562_rows.py h1_pseudobulk.npz grn_h1_pc0.top_edges.csv K562_gwps_normalized_bulk_01.h5ad k562_rows.npz
python run_grn.py k562_rows.npz grnboost2 k562_grnboost2.npz 0 0 4
python eval_grn_h1.py --net grnboost2=h1_grnboost2.npz --out res_gb   # --pseudobulk/--genes-from/--screen-dir default to the project folder
python paired_check.py res_gb/grn_h1_scores.csv grnboost2_out grnboost2_in
```

- **Input data.** `competition_train.h5` (15.5 GB) comes from the public 2025 support-set zip. One member can be range-extracted without downloading the whole 8.7 GB file.
- **Networks.** `networks/*.npz` hold `imp` (regulator × target importance), `genes`, and `rho` (Spearman correlation).
- **Running bigger versions.** A Colab High-RAM runtime can run full GENIE3 (1,000 trees) or more genes with the same scripts; time grows roughly linearly with genes × trees. Given the results above, this is not worth doing for prediction.
