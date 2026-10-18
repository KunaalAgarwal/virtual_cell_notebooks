"""Paired bootstrap and held-out tuning: does adding network similarity to screen similarity help?

For each candidate family (screens alone vs screens + network), the setting (m, shrink) is
tuned on a random half of the knockdowns and scored on the other half, 2,000 times; the
difference is also bootstrapped over knockdowns for one fixed best setting per family.

    python paired_check.py grn_h1_scores.csv net_name [net_name ...]
"""
import sys
import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1])
df = df[df.part == "B"]
P = df.pivot_table(index="pert", columns=["method", "top"], values="pearson_all")
rng = np.random.default_rng(0)


def family(prefix):
    return [c for c in P.columns if c[0] in (prefix, prefix + "_shrink0.5")]


fams = {"screens": family("screens")}
for nm in sys.argv[2:]:
    fams[f"screens+{nm}"] = family(f"screens+{nm}")
    fams[nm] = family(nm)
base = fams["screens"]
n = len(P)
for name, cols in fams.items():
    if name == "screens":
        continue
    wins, gains = 0, []
    for _ in range(2000):
        idx = rng.permutation(n); a, b = idx[: n // 2], idx[n // 2:]
        cb = max(base, key=lambda c: P[c].iloc[a].mean())
        cf = max(cols, key=lambda c: P[c].iloc[a].mean())
        g = P[cf].iloc[b].mean() - P[cb].iloc[b].mean()
        gains.append(g); wins += g > 0
    best_b = max(base, key=lambda c: P[c].mean()); best_f = max(cols, key=lambda c: P[c].mean())
    d = (P[best_f] - P[best_b]).to_numpy()
    bs = [d[rng.integers(0, n, n)].mean() for _ in range(2000)]
    print(f"{name:28s} best {best_f} {P[best_f].mean():.3f} vs screens {best_b} {P[best_b].mean():.3f}; "
          f"diff {d.mean():+.3f} CI [{np.percentile(bs, 2.5):+.3f}, {np.percentile(bs, 97.5):+.3f}]; "
          f"held-out half gain {np.mean(gains):+.3f}, better in {wins / 20:.0f}% of splits")
