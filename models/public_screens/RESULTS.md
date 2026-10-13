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

## Caveats
- C, m and shrinkage were picked on these same H1 knockdowns (few settings, gains stable across them).
- The 2025 support-set `competition_train.h5` stores log1p(counts); `celleval_screen_h1.py` converts back to counts. The earlier v2 run used the raw-count release, and its baselines reproduce here exactly.
- H1 is a stand-in: for 2026 the borrowed effects would be H1's (or the screens'), moved into contexts A/B/C.

## Files
See [README.md](README.md). Explainer page: https://claude.ai/artifact/6znsVo6PE54dxpW3fpWJLb
