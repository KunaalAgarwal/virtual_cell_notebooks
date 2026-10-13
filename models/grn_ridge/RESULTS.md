# GRN v2: all-genes ridge, tested against H1's real knockdowns

## Bottom line

A network learned from control cells **does not predict knockdown effects** in H1. Its predicted changes correlate with the real ones at about 0.00, which is no better than chance. Simple co-expression with the knocked-down gene also scores about 0.00.

The only thing that carried signal was the trivial baseline, "every knockdown does what the average knockdown does", with a mean Pearson of 0.19 (0.33 on the strongest knockdowns). Adding the network's shift on top of that baseline made it slightly worse.

So the plan from `grn_linear_model.pdf`, run as specified (edges from controls, then propagate the knockdown), does not work on the one dataset where we can check it. The co-variation among unperturbed cells is not the causal structure a knockdown follows.

## What was built (`v2/`)

| File | What it does |
|---|---|
| `grn_ridge.py` | Fits every "gene ~ all other genes" ridge in one matrix inverse per fold. Uses grouped CV (guides, or batches for H1) and picks λ on held-out R². Outputs slopes (`coef_mean`), `coef_std`, `sign_agree` (the share of folds that agree on sign), mean \|SHAP\| (exact for linear), per-fold-normalized `shap_share`, and held-out R² for every gene. Optionally regresses out library size, MT fraction and the top PCs first. |
| `h1_pseudobulk.py` | Streams the 15 GB H1 file and computes each of the 150 knockdowns' real change versus controls **from the same batch**, in log1p(CP10k). |
| `eval_h1.py` | Pushes each knockdown through a network ((I − W)δ = push, iterated; converges because the spectral radius is below 1) and scores it against the real change. |
| `h1/grn_h1_pc0.npz`, `h1/grn_h1_pc10.npz` | The H1 networks: 10,993 x 10,993, fitted on 38,176 control cells with 5 folds over the 48 batches. One is raw; the other has cell state removed (10 PCs plus covariates). |
| `h1/*.top_edges.csv` | The top 20 sign-stable regulators per gene. |
| `h1/h1_knockdown_scores.csv`, `.summary.csv` | Per-knockdown scores for every method. |

Speed: about 2 minutes per fold for the full 11k-gene network. The per-gene loop would have taken weeks.

## Numbers

**How well the networks predict held-out control cells (R² per gene, averaged):**

| Network | H1 | 2026 context A (6k genes) |
|---|---|---|
| Raw | 0.08 | 0.09 |
| Cell state removed | 0.02 | 0.02 |

Most of the predictability in control cells comes from shared programs (cell cycle, quality, differentiation). Once those are removed, little is left.

**How well each method predicts real H1 knockdowns** (147 knockdowns, mean over knockdowns, excluding the knocked-down gene itself):

| Method | Pearson (all genes) | Pearson (100 most-changed genes) | Direction right (50 most-changed genes) | Pearson, strongest 25% of knockdowns |
|---|---|---|---|---|
| Mean effect of the other knockdowns | **0.19** | **0.28** | **0.65** | **0.33** |
| Mean effect + network (cell state removed) | 0.19 | 0.27 | 0.65 | 0.33 |
| Mean effect + network (raw) | 0.18 | 0.27 | 0.65 | 0.32 |
| Co-expression | 0.00 | 0.00 | 0.51 | 0.01 |
| Network (raw), 20 steps | 0.00 | 0.01 | 0.51 | 0.00 |
| Network (cell state removed), 20 steps | −0.02 | −0.02 | 0.50 | −0.04 |

Neither the number of propagation steps (1, 3 or 20) nor the assumed size of the knockdown (the real value vs the H1 median of −0.63) changes the picture.

## What this means for the approach

1. **The causal signal has to come from perturbation data, not control cells.** The public knockdown screens (Replogle K562/RPE1, Nadig HepG2/Jurkat, X-Atlas HCT116/HEK293T, Feng iPSC) directly show which genes move when a gene is knocked down.
2. **The "average knockdown" effect is real and is the bar to beat.** A large part of every knockdown's change is shared: a generic stress or knockdown response. A strong submission probably starts from a shared response like this, adapted to each context, plus gene-specific terms learned from perturbation screens.
3. **The control-cell network may still help as a feature.** For example, it could say which genes are co-regulated *in a given context* when moving a knockdown effect learned in K562 into context A. The 2026 data thread found A, B and C differ a lot, so that adaptation step matters. It just can't be the source of the effect on its own.

## Caveats

- This tests one network design (linear ridge on log-normalized single cells). The tree models from v1 were not run through the knockdown test. Given that co-expression also scores about 0, they are unlikely to change the conclusion, but that hasn't been checked.
- The comparison is on pseudobulk means, not the challenge's own `cell-eval` metrics.

## Round 2 (2026-10-13): the improvements, applied and run

Three of the proposed improvements were run on H1. Scripts: `metacell_grn.py`, `loo_transfer_h1.py`, `h1_extract_perts.py`, `celleval_h1.py`. Outputs in `h1/` (`h1_knockdown_scores_meta30.*`, `loo_transfer.csv`, `celleval/`).

**1. Metacell smoothing before fitting the network.** Control cells were averaged over their 30 nearest neighbours in PC space, then the ridge network was refitted (`h1/grn_h1_meta30.npz`). Held-out R² rises to 0.74 because smoothing removes noise, but the network still predicts real knockdowns at Pearson −0.005 (fixed drop, 20 steps). Adding it to the mean effect makes things worse (0.05 vs 0.19). Denoising does not create causal signal.

**2. Leave-one-out effect transfer.** Predict a held-out knockdown as a similarity-weighted average of the m most similar other knockdowns, where similarity comes from control cells (co-expression, or the network's regulator or target profiles). Pearson over all genes, mean of 147 knockdowns:

| Similarity used to pick effects to borrow | m=5 | m=15 | m=40 | all 146 (uniform) |
|---|---|---|---|---|
| Uniform (= mean effect) | | | | **0.189** |
| Co-expression in controls | 0.123 | 0.151 | 0.178 | |
| Network, target profile | 0.124 | 0.157 | 0.184 | |
| Network, regulator profile | 0.074 | 0.140 | 0.170 | |
| True effect similarity (oracle) | **0.461** | 0.454 | 0.443 | |

Control-cell similarity is never better than averaging everything. The oracle row shows the headroom: if gene similarity could be learned from perturbation data, this transfer scheme reaches ~0.46. The signal has to come from screens.

**3. Official `cell-eval` metrics** on 40 real H1 knockdowns (400 cells each, 2,000 controls). Each submission is built the way a real one would be: resample control cells and shift them by the predicted change (`celleval_h1.py`). mae lower is better; the other two higher.

| Method | overlap_at_N | mae | discrimination_score_l1 |
|---|---|---|---|
| Resampled controls, no change | 0.000 | 0.049 | 0.521 |
| Control mean in every cell (repo baseline) | 0.146 | 0.140 | 0.520 |
| Controls + mean effect of other knockdowns | 0.044 | 0.049 | 0.521 |
| Controls + control-cell network (metacell) | 0.108 | 0.062 | 0.527 |
| Controls + the true change (oracle) | 0.548 | 0.031 | 0.904 |

Discrimination is at chance (0.5) for everything but the oracle: a change shared by all knockdowns cannot tell knockdowns apart, so the mean effect does not score on it even though it wins on pseudobulk Pearson. The network's small overlap gain (0.11) comes from it flagging more genes as changed, with a worse mae. The repo's control-mean baseline has the best overlap of the non-oracle methods only because constant cells make every gene look differentially expressed; its mae is 3x worse.

**Not run:** learning gene-specific effects from public screens (Replogle K562/RPE1 etc.). zenodo.org, huggingface.co, GEO/NCBI, figshare and dataverse are all blocked from this environment; only the Arc atlas bucket on storage.googleapis.com is reachable. This needs Colab or a thread in an environment that allows those hosts.

**Housekeeping:** `h1/grn_h1_pc0.npz` and `h1/grn_h1_pc10.npz` (1.8 GB each) are no longer in the shared folder; only their `top_edges.csv` remain. `grn_ridge.py` regenerates them in about 10 minutes each.
