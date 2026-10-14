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

## Round 4 (2026-10-14): hold out a whole cell line, apply the audit's improvements

Follows `docs/IMPROVEMENTS_AUDIT.md` in priority order. Scripts: `xline_benchmark.py` (pseudobulk scores, nested tuning, ablations), `celleval_xline.py` (cell-eval). Outputs in the project folder `public-screens/results_xline/`.

**Benchmark (audit idea 1).** Hold out HepG2 or Jurkat from the 2025 support set [Nadig2025], give the methods only its control cells, and predict its 68 knockdowns from H1's real effects for *other* genes plus the K562 genome-wide screen [Replogle2022], which covers the target itself. This is the 2026 situation: 13 of 300 targets are in H1, 272 are in K562. Cross-line splits follow [AhlmannEltze2025; Adduri2025]. Settings are tuned on one held-out line and scored on the other (nested selection [VarmaSimon2006]); every comparison is a paired bootstrap over knockdowns [Efron1979].

**What was tried (960 settings per line).**
- Per-knockdown amplitude (idea 3): C_k = C × s_k^γ, with s_k the K562 Anderson–Darling DE-gene count relative to its median, so knockdowns with little evidence of an effect are shrunk to the shared response (empirical-Bayes flavour [JamesStein1961; EfronMorris1975; Stephens2017; Zhu2019]).
- Expression gating (idea 6): no predicted change on genes detected in < 1% of the held-out line's controls [Hafemeister2019; Svensson2020].
- Low-rank similarity (idea 4.1): correlation in the top 50 PCs of the K562 screen [AhlmannEltze2025; CandesRecht2009].
- Stacking (idea 9): add the target's own K562 profile, converted to approximate log fold changes, to the borrowed H1 effect [Wolpert1992].
- Shared response from H1 or none (idea 5, simplest form), m = 10 or 20 neighbours, C from 1 to 16.

**Noise in this benchmark is high.** HepG2 has a median of 38 cells per knockdown and Jurkat 88. Split-half correlation of the real effects is median 0.17 (HepG2) and 0.23 (Jurkat); 21 and 15 knockdowns are below 0.1 (`results/xline_split_half_reliability.csv`). Even the same gene's *real* effect from the other line only reaches discrimination 0.54–0.57, and the oracle's overlap is 0.13–0.22. Differences are small in absolute terms; the bootstraps decide.

**Pseudobulk scores, tuned on the other line** (`results_xline/summary_by_objective.csv`; Δ vs ×8 with 95% CI)

| Held out | Method | Discrimination | Pearson | Pearson, specific part | MAE |
|---|---|---|---|---|---|
| HepG2 | Round-3 ×8 reference | 0.529 | 0.088 | 0.044 | 0.053 |
| HepG2 | Tuned for discrimination | **0.570** (+0.042 [+0.019, +0.067]) | 0.090 | 0.059 | 0.078 |
| HepG2 | Tuned for Pearson | 0.514 | **0.138** | **0.104** (+0.060 [+0.043, +0.077]) | **0.039** |
| HepG2 | Same gene, other line (not available in 2026) | 0.541 | 0.130 | 0.105 | 0.047 |
| Jurkat | Round-3 ×8 reference | 0.526 | 0.092 | 0.040 | 0.046 |
| Jurkat | Tuned for discrimination + Pearson | **0.549** (+0.023 [+0.007, +0.040]) | 0.108 | 0.072 | 0.057 |
| Jurkat | Tuned for Pearson | 0.518 | **0.162** | **0.122** (+0.082 [+0.063, +0.103]) | **0.032** |
| Jurkat | Same gene, other line | 0.550 | 0.134 | 0.110 | 0.044 |

**cell-eval on the held-out line** (`results_xline/celleval_*/`)

| Held out | Method | overlap_at_N ↑ | discrimination ↑ | mae ↓ |
|---|---|---|---|---|
| Jurkat | Control mean (repo baseline) | 0.019 | 0.510 | 0.085 |
| Jurkat | Round-3 ×8 reference | 0.000 | 0.529 | 0.054 |
| Jurkat | Tuned: amplitude + gate + K562 term (C 16) | **0.032** | **0.566** | 0.063 |
| Jurkat | Tuned: same, C 8 (disc − mae objective) | 0.030 | 0.546 | **0.048** |
| Jurkat | Same gene, other line | 0.003 | 0.568 | 0.050 |
| Jurkat | Oracle | 0.220 | 1.000 | 0.026 |
| HepG2 | Control mean (repo baseline) | **0.025** | 0.508 | 0.070 |
| HepG2 | Round-3 ×8 reference | 0.000 | 0.540 | 0.063 |
| HepG2 | Tuned for discrimination (no gate) | 0.000 | **0.575** | 0.087 |
| HepG2 | Tuned: C 4, amplitude + gate + K562 term | 0.011 | 0.540 | **0.053** |
| HepG2 | Same gene, other line | 0.000 | 0.557 | 0.057 |
| HepG2 | Oracle | 0.133 | 1.000 | 0.034 |

Paired bootstraps (68 knockdowns), Jurkat, tuned C 16 vs ×8: overlap +0.032 [+0.020, +0.047], discrimination +0.037 [+0.018, +0.057], mae +0.009. Tuned C 8 vs ×8: overlap +0.030, discrimination +0.017 [+0.005, +0.029], mae −0.006. Against the control-mean baseline the discrimination gains are not significant on Jurkat (+0.036 to +0.056, CIs cross 0); on HepG2 the discrimination-tuned setting is (+0.067 [+0.011, +0.120]), but the control mean keeps the best overlap there.

**What each improvement did**
- **Gating is required for overlap.** Every ungated method scores overlap 0.000 on both lines: changes on genes that are off in the new line produce infinite fold changes that fill the top of the DE ranking.
- **Per-knockdown amplitude helps.** Both lines independently chose the strength-scaled amplitude with γ = 2, and on Jurkat it adds discrimination over a fixed C (pseudobulk 0.549 vs 0.542; specific Pearson 0.072 vs 0.060).
- **Stacking the K562 profile helps most for Pearson and MAE**: the specific-part Pearson doubles or triples and reaches the "same gene, other line" level. Removing it raises discrimination on HepG2 only.
- **Low-rank similarity did not help**: PCA similarity was never chosen and scored at or below raw correlation.
- **The H1 shared response** was kept on Jurkat and dropped on HepG2 by the tuning, so it does not transfer reliably.

**Bottom line.** Across cell lines the signal is much weaker than within H1 (discrimination ~0.55 instead of 0.66), and the benchmark's own noise limits what can be shown. The improved settings beat the round-3 ×8 reference significantly on Jurkat and match or beat it on HepG2, and they reach the level of knowing the same gene's real effect in another line. They do not clearly beat the control-mean baseline across lines. Recommended default for 2026: C = 8 with strength-scaled amplitude (γ = 2), 1% expression gate, K562 term weight 1, H1 shared response.

**Not done yet:** ideas 5 (full shared/specific decomposition), 7 (realistic cell generation), 8 (State zero-shot), 4.2–4.4 (prior-knowledge and learned similarity).

## Caveats
- References in brackets are listed in `docs/IMPROVEMENTS_AUDIT.md`.
- C, m and shrinkage were picked on these same H1 knockdowns; held-out tuning (see Noise checks) gives the same gains.
- The 2025 support-set `competition_train.h5` stores log1p(counts); `celleval_screen_h1.py` converts back to counts. The earlier v2 run used the raw-count release, and its baselines reproduce here exactly.
- H1 is a stand-in: for 2026 the borrowed effects would be H1's (or the screens'), moved into contexts A/B/C.

## Files
See [README.md](README.md). Explainer page: https://claude.ai/artifact/6znsVo6PE54dxpW3fpWJLb
