# Public-screen gene similarity: results on H1 (2026-10-13)

## Idea
Predict a held-out H1 knockdown as a weighted average of the *other* H1 knockdowns' real changes. The weights come from how alike the two target genes look in public Perturb-seq screens: correlation of their z-scored pseudobulk profiles in Replogle 2022 K562 genome-wide (133 of 147 H1 knockdowns covered) and RPE1 (48), plus Nadig 2025 HepG2 and Jurkat (48 each, from the 2025 VCC support set). Top 10 most similar, shrunk halfway to the mean effect. Knockdowns in no screen fall back to the mean effect.

## Pseudobulk Pearson (147 knockdowns, 10,993 genes, leave-one-out)
| Method | Pearson all | Pearson top100 | Sign top50 |
|---|---|---|---|
| Mean effect (uniform) | 0.189 | 0.275 | 0.651 |
| K562 similarity | 0.243 | 0.329 | 0.682 |
| K562+RPE1+HepG2+Jurkat similarity | **0.249** | **0.338** | **0.688** |
| Direct (own K562 profile, scaled) | 0.196 | 0.314 | 0.671 |
| Oracle similarity | 0.461 | 0.608 | 0.835 |

Paired on the 133 covered knockdowns: 0.262 vs 0.201 uniform, better on 77 of 133 (Wilcoxon p = 4e-4). Stable over m = 5–40 and shrinkage 0.5–1.

## cell-eval, 40 real H1 knockdowns (vcc profile)
"aC" = mean effect + C × (screen-weighted − mean effect), i.e. the knockdown-specific part amplified C times. "sS" also scales the shared part.

| Method | overlap_at_N ↑ | mae ↓ | discrimination_l1 ↑ |
|---|---|---|---|
| Resampled controls | 0.000 | 0.049 | 0.521 |
| Control mean (repo baseline) | 0.146 | 0.140 | 0.520 |
| Mean effect | 0.044 | 0.049 | 0.521 |
| Mean effect ×8 (s8 a0) | 0.203 | 0.067 | 0.509 |
| Null: another knockdown's specific part ×8 | 0.163 | 0.078 | 0.543 |
| Screen similarity | 0.084 | 0.048 | 0.543 |
| Screen similarity a4 | 0.190 | 0.050 | 0.628 |
| **Screen similarity a8** | 0.238 | 0.069 | **0.661** |
| Screen similarity s4 a8 | 0.259 | 0.077 | 0.646 |
| Screen similarity a16 | **0.262** | 0.125 | 0.657 |
| Oracle | 0.548 | 0.031 | 0.904 |

a4 beats the control-mean baseline on all three metrics; a8 is the best balance. The shuffled null shows the discrimination gain (0.54 → 0.66) and most of the overlap gain come from the screens, not from amplification alone. Amplifying alone raises overlap (more genes pass the DE test) but never discrimination.

## Noise checks (2026-10-13)
| Check | Result |
|---|---|
| Pearson gain over mean effect, bootstrap over 147 knockdowns | +0.060, 95% CI [+0.036, +0.086] |
| Pearson with settings tuned on half the knockdowns, scored on the other half (2000 splits) | 0.244 vs 0.188; better in 100% of splits |
| cell-eval ×8 minus control-mean baseline, bootstrap over 40 knockdowns | overlap +0.092 [+0.042, +0.141]; discrimination +0.141 [+0.059, +0.229]; mae −0.071 [−0.077, −0.065] |
| cell-eval ×4 minus control-mean baseline | overlap +0.043 [−0.024, +0.110] (not significant); discrimination +0.108 [+0.039, +0.184]; mae −0.090 |
| cell-eval ×8 minus shuffled null ×8 | overlap +0.075 [+0.033, +0.120]; discrimination +0.118 [+0.027, +0.211] |
| Amplification tuned on 20 knockdowns, scored on the other 20 | held-out overlap 0.260 vs 0.147, discrimination 0.648 vs 0.519; beats the baseline in 99.8% of splits |
| Re-run with 3 random seeds (control-cell sampling and null shuffle) | ×8: overlap 0.236–0.241, discrimination 0.658–0.661. Null ×8: overlap 0.163–0.196, discrimination 0.510–0.570 |
| Chance discrimination for 40 knockdowns (random ranks) | about 0.41–0.59 (95%) |
| Noise in the H1 ground truth: split-half correlation of real effects (200 vs 200 cells, 40 knockdowns) | median 0.58; 8 weak knockdowns below 0.12, which no method can predict. Rough Pearson ceiling about 0.77 |

×8 beats the control-mean baseline and the shuffled null beyond noise on overlap and discrimination; ×4's overlap gain is within noise. Per-knockdown split-half values are in `results/h1_split_half_reliability.csv` (project folder).

## Caveats
- C, m and shrinkage were picked on these same H1 knockdowns; held-out tuning (see Noise checks) gives the same gains.
- The 2025 support-set `competition_train.h5` stores log1p(counts); `celleval_screen_h1.py` converts back to counts. The earlier v2 run used the raw-count release, and its baselines reproduce here exactly.
- H1 is a stand-in: for 2026 the borrowed effects would be H1's (or the screens'), moved into contexts A/B/C.

## Files
See [README.md](README.md). Explainer page: https://claude.ai/artifact/6znsVo6PE54dxpW3fpWJLb
