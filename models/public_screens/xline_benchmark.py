"""Leave-one-cell-line-out benchmark for screen-similarity effect borrowing.

The 2026 challenge asks for knockdown effects in cell lines we have no knockdowns for, and
13 of its 300 targets were knocked down in H1, while 272 are in the Replogle K562 screen. This
benchmark mimics that (IMPROVEMENTS_AUDIT.md idea 1, after Ahlmann-Eltze et al. 2025, Nature
Methods, and Adduri et al. 2026, Cell): hold out a whole cell line from the 2025 VCC support set
(Nadig et al. 2025: HepG2 or Jurkat), give the methods only its control cells, and predict its
68 knockdowns from

  * H1's real knockdown effects for *other* genes (the same gene's H1 effect is never used), and
  * the K562 genome-wide screen (Replogle et al. 2022, Cell), which covers the target gene itself.

Settings are tuned on one held-out line and scored on the other (nested selection, Varma and
Simon 2006). Fast pseudobulk scores reproduce what cell-eval does on means: the predicted cells
are sampled control cells shifted by the predicted change and clipped at zero, and their mean
minus the control mean is the predicted effect. Discrimination is cell-eval's L1 rank,
target gene excluded. Methods are compared with paired bootstraps over knockdowns
(Efron 1979). cell-eval's overlap needs a DE test and is run separately (celleval_xline.py).

    python xline_benchmark.py --lines hepg2.h5 jurkat.h5 --h1 h1_pseudobulk.npz \
        --k562 K562_gwps_normalized_bulk_01.h5ad --out results_xline/
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

CTRL = "non-targeting"
H1_SELF_DROP = -0.633


# ----------------------------------------------------------------------------- data

def load_line(path, n_ctrl_sample=150, seed=0):
    """Target line: log-normalized cells (support-set .h5). Returns control matrix, real effects."""
    a = ad.read_h5ad(path)
    X = np.asarray(a.X.toarray() if hasattr(a.X, "toarray") else a.X, dtype=np.float32)
    t = a.obs["target_gene"].astype(str).to_numpy()
    ctrl = X[t == CTRL]
    perts = sorted(set(t) - {CTRL})
    cmean = ctrl.mean(0)
    real = np.stack([X[t == p].mean(0) - cmean for p in perts])
    n = np.array([(t == p).sum() for p in perts])
    det = (ctrl > 0).mean(0)
    rng = np.random.default_rng(seed)
    sample = ctrl[rng.choice(len(ctrl), min(n_ctrl_sample, len(ctrl)), replace=False)]
    return dict(name=Path(path).stem, genes=a.var_names.to_numpy(), perts=perts, real=real, n=n,
                ctrl_mean=cmean, ctrl_sample=sample, det=det)


def load_k562(path):
    """K562 GWPS z-scored pseudobulk -> per-gene profiles, approximate log changes, metadata."""
    a = ad.read_h5ad(path)
    Z = np.asarray(a.X, dtype=np.float64)
    Z[~np.isfinite(Z)] = 0
    np.clip(Z, -10, 10, out=Z)
    tgt = a.obs.index.str.split("_").str[1].to_numpy()
    df = pd.DataFrame(Z).groupby(tgt).mean()
    meta = a.obs.assign(tgt=tgt).groupby("tgt")[["energy_test_p_value", "anderson_darling_counts",
                                                  "fold_expr", "num_cells_filtered"]].mean()
    # z = change / control sd in normalized UMI units, so change / mean ~ log fold change.
    var = a.var
    lfc_scale = (var["std"] / var["mean"].clip(lower=1e-3)).to_numpy()
    return dict(Z=df, genes=var["gene_name"].to_numpy(), meta=meta, lfc_scale=lfc_scale)


def corr_rows(A, B):
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    A /= np.linalg.norm(A, axis=1, keepdims=True) + 1e-12
    B /= np.linalg.norm(B, axis=1, keepdims=True) + 1e-12
    return A @ B.T


# ----------------------------------------------------------------------------- features

def build_features(line, h1, k562, n_pcs=50):
    genes = line["genes"]; gi = {g: i for i, g in enumerate(genes)}
    T = line["perts"]
    h1_perts = [str(p) for p in h1["perts"]]
    D = h1["delta"].astype(np.float64)                         # H1 effects, same gene order
    Zdf = k562["Z"]
    in_k = [p in Zdf.index for p in T]
    src = [p for p in h1_perts if p in Zdf.index]               # H1 genes with a K562 profile
    Si = np.array([h1_perts.index(p) for p in src])
    readout_keep = ~np.isin(k562["genes"], list(set(T) | set(src)))
    Zall = Zdf.to_numpy()[:, readout_keep]
    zT = np.stack([Zdf.loc[p].to_numpy()[readout_keep] if ok else np.zeros(readout_keep.sum())
                   for p, ok in zip(T, in_k)])
    zS = Zdf.loc[src].to_numpy()[:, readout_keep]
    sim = {"corr": corr_rows(zT, zS)}
    # Low-rank similarity: top PCs of the whole genome-wide screen (audit idea 4.1;
    # Ahlmann-Eltze et al. 2025; Candes and Recht 2009).
    mu = Zall.mean(0)
    _, _, Vt = np.linalg.svd(Zall[np.random.default_rng(0).choice(len(Zall), 4000, replace=False)] - mu,
                             full_matrices=False)
    V = Vt[:n_pcs].T
    sim["pca"] = corr_rows((zT - mu) @ V, (zS - mu) @ V)
    for s in sim.values():
        s[~np.array(in_k)] = np.nan
    # same gene never borrowed from H1
    for i, p in enumerate(T):
        if p in src:
            for s in sim.values():
                s[i, src.index(p)] = np.nan
    # K562 direct: the target's own profile as an approximate log change on line genes
    kg = {g: c for c, g in enumerate(k562["genes"])}
    cols = np.array([gi[g] for g in k562["genes"] if g in gi]); kcols = np.array([kg[g] for g in k562["genes"] if g in gi])
    direct = np.zeros((len(T), len(genes)))
    for i, p in enumerate(T):
        if in_k[i]:
            direct[i, cols] = Zdf.loc[p].to_numpy()[kcols] * k562["lfc_scale"][kcols]
    np.clip(direct, -3, 3, out=direct)
    meta = k562["meta"].reindex(T)
    return dict(sim=sim, H=D[Si], h1_mean=D.mean(0), direct=direct, meta=meta, in_k=np.array(in_k),
                gi=gi, k562_mean=None)


# ----------------------------------------------------------------------------- predictor

DEFAULT = dict(sim="corr", m=10, alpha=0.5, C=8.0, shared="h1", amp="fixed", gamma=0.0,
               gate=0.0, w_direct=0.0, self_drop="h1")


def predict(line, F, p):
    """Predicted log change for every target knockdown under settings p."""
    S = F["sim"][p["sim"]]; H = F["H"]
    shared = {"h1": F["h1_mean"], "none": np.zeros_like(F["h1_mean"])}[p["shared"]]
    out = np.zeros((len(line["perts"]), H.shape[1]))
    # Per-knockdown amplitude (audit idea 3): C_k = C * s_k, s_k from K562 evidence of an effect.
    # Empirical-Bayes flavour (James and Stein 1961; Efron and Morris 1975; Stephens 2017):
    # knockdowns with little evidence are shrunk to the shared response.
    adc = np.log1p(F["meta"]["anderson_darling_counts"].fillna(0).to_numpy())
    strength = adc / np.nanmedian(adc[adc > 0]) if (adc > 0).any() else np.ones(len(adc))
    for i in range(len(out)):
        s = S[i]
        ok = np.flatnonzero(np.isfinite(s))
        spec = np.zeros(H.shape[1])
        if len(ok):
            order = ok[np.argsort(-s[ok])[: p["m"]]]
            w = np.clip(s[order], 0, None)
            if w.sum() > 0:
                borrowed = (w[:, None] * H[order]).sum(0) / w.sum()
                spec = p["alpha"] * (borrowed - F["h1_mean"])
        Ck = p["C"]
        if p["amp"] == "strength" and F["in_k"][i]:
            Ck = p["C"] * strength[i] ** p["gamma"]
        d = shared + Ck * spec + p["w_direct"] * F["direct"][i]
        out[i] = d
    if p["gate"] > 0:   # audit idea 6: no change on genes silent in this context (Hafemeister and Satija 2019)
        out[:, line["det"] < p["gate"]] = 0
    gi = F["gi"]
    fe = F["meta"]["fold_expr"].to_numpy()
    for i, g in enumerate(line["perts"]):
        if p["self_drop"] == "k562" and np.isfinite(fe[i]):
            # target's own drop from K562 knockdown strength, in this line's log units
            out[i, gi[g]] = np.log1p(np.expm1(line["ctrl_mean"][gi[g]]) * max(fe[i], 0.05)) - line["ctrl_mean"][gi[g]]
        else:
            out[i, gi[g]] = H1_SELF_DROP
    return out


# ----------------------------------------------------------------------------- scoring

def pred_effect(line, delta):
    """What cell-eval sees: mean of max(control + delta, 0) minus control mean."""
    C = line["ctrl_sample"]
    base = C.mean(0)
    return np.stack([np.clip(C + d.astype(np.float32), 0, None).mean(0) - base for d in delta])


def score(line, delta):
    P = pred_effect(line, delta); R = line["real"]; gi = {g: i for i, g in enumerate(line["genes"])}
    n = len(R)
    disc, pear, spec, mae = [], [], [], []
    Rm = R.mean(0); Pm = P.mean(0)
    for i, g in enumerate(line["perts"]):
        k = np.ones(R.shape[1], bool); k[gi[g]] = False
        d = np.abs(R[:, k] - P[i, k]).sum(1)
        disc.append(1 - (d < d[i]).sum() / (n - 1))
        pear.append(np.corrcoef(P[i, k], R[i, k])[0, 1] if P[i, k].std() > 0 else 0.0)
        a, b = P[i, k] - Pm[k], R[i, k] - Rm[k]
        spec.append(np.corrcoef(a, b)[0, 1] if a.std() > 0 else 0.0)
        mae.append(np.abs(P[i] - R[i]).mean())
    return pd.DataFrame({"pert": line["perts"], "discrimination": disc, "pearson": pear,
                         "pearson_specific": spec, "mae": mae})


def boot_diff(a, b, B=5000, seed=0):
    x = np.asarray(a) - np.asarray(b)
    bs = x[np.random.default_rng(seed).integers(0, len(x), (B, len(x)))].mean(1)
    return x.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5)


# ----------------------------------------------------------------------------- main

GRID = {
    # alpha is folded into C (only alpha * C matters); self_drop only moves the target gene,
    # which discrimination and Pearson exclude, so it is compared separately.
    "sim": ["corr", "pca"], "m": [10, 20], "alpha": [0.5], "C": [1, 2, 4, 8, 16],
    "shared": ["h1", "none"], "amp": ["fixed", "strength"], "gamma": [0.5, 1.0, 2.0],
    "gate": [0.0, 0.01], "w_direct": [0.0, 0.5, 1.0], "self_drop": ["k562"],
}


def configs():
    keys = list(GRID)
    for vals in itertools.product(*GRID.values()):
        p = dict(zip(keys, vals))
        if p["amp"] == "fixed" and p["gamma"] != GRID["gamma"][0]:
            continue
        yield p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines", nargs=2, required=True)
    ap.add_argument("--h1", required=True)
    ap.add_argument("--k562", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--objective", default="discrimination+pearson",
                    help="what nested tuning maximizes: discrimination, pearson, or discrimination+pearson")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    h1 = np.load(a.h1, allow_pickle=True)
    k562 = load_k562(a.k562)
    lines = [load_line(p) for p in a.lines]
    feats = [build_features(L, h1, k562) for L in lines]
    for L, F in zip(lines, feats):
        print(f"{L['name']}: {len(L['perts'])} knockdowns, {F['in_k'].sum()} in K562, "
              f"median {int(np.median(L['n']))} cells each")

    cfgs = list(configs())
    print(f"{len(cfgs)} configurations")
    allres = []
    for L, F in zip(lines, feats):
        for j, p in enumerate(cfgs):
            s = score(L, predict(L, F, p))
            allres.append({"line": L["name"], "cfg": j, **{k: s[k].mean() for k in
                           ["discrimination", "pearson", "pearson_specific", "mae"]}})
            if j % 200 == 0:
                print(L["name"], j, flush=True)
    grid = pd.DataFrame(allres)
    grid.to_csv(out / "grid.csv", index=False)
    json.dump(cfgs, open(out / "configs.json", "w"))

    def obj(df):
        if a.objective == "discrimination+pearson":
            return df["discrimination"] + df["pearson"]
        return df[a.objective]

    rows, per = [], []
    baselines = {
        "zero (resampled controls)": None,
        "H1 mean effect": dict(DEFAULT, C=0.0),
        "reference x8 (round 3)": DEFAULT,
        "reference x8, K562 self-drop": dict(DEFAULT, self_drop="k562"),
        "shuffled x8 (null)": "shuffle",
    }
    for ti, (L, F) in enumerate(zip(lines, feats)):
        other = lines[1 - ti]["name"]
        g = grid[grid.line == other]
        best = cfgs[int(g.loc[obj(g).idxmax(), "cfg"])]
        # ablations: best config, one improvement switched off at a time, all tuned on the other line
        methods = dict(baselines)
        methods["tuned on other line (all improvements)"] = best
        methods["tuned, fixed amplitude"] = _best_where(g, cfgs, obj, amp="fixed")
        methods["tuned, no K562 direct term"] = _best_where(g, cfgs, obj, w_direct=0.0)
        methods["tuned, corr similarity"] = _best_where(g, cfgs, obj, sim="corr")
        methods["tuned, no gating"] = _best_where(g, cfgs, obj, gate=0.0)
        # ceilings: the real effect, and the same gene's real effect in the other held-out line
        O = lines[1 - ti]
        methods["same gene, other line (not available in 2026)"] = "other"
        methods["oracle (true effect)"] = "oracle"
        scores = {}
        for name, p in methods.items():
            if p == "oracle":
                d = L["real"].copy()
            elif p == "other":
                d = np.stack([O["real"][O["perts"].index(q)] if q in O["perts"] else np.zeros(len(L["genes"]))
                              for q in L["perts"]])
            elif p is None:
                d = np.zeros((len(L["perts"]), len(L["genes"])))
            elif isinstance(p, str) and p == "shuffle":
                d = predict(L, F, DEFAULT)
                perm = np.random.default_rng(1).permutation(len(d))
                d = d[perm]
            else:
                d = predict(L, F, p)
            s = score(L, d); s["method"] = name; s["line"] = L["name"]
            scores[name] = s; per.append(s)
        ref = scores["reference x8 (round 3)"]
        for name, s in scores.items():
            r = {"line": L["name"], "method": name}
            for mtr in ["discrimination", "pearson", "pearson_specific", "mae"]:
                r[mtr] = s[mtr].mean()
                dm, lo, hi = boot_diff(s[mtr], ref[mtr])
                r[f"{mtr}_vs_ref"] = f"{dm:+.3f} [{lo:+.3f},{hi:+.3f}]"
            rows.append(r)
        json.dump({"tuned_on": other, "best": best, "ablations": {k: v for k, v in methods.items()
                   if isinstance(v, dict)}}, open(out / f"chosen_{L['name']}.json", "w"), indent=1)
    res = pd.DataFrame(rows)
    res.to_csv(out / "summary.csv", index=False)
    pd.concat(per).to_csv(out / "per_pert.csv", index=False)
    pd.set_option("display.width", 250)
    print(res.round(3).to_string())


def _best_where(g, cfgs, obj, **fixed):
    ok = [j for j, p in enumerate(cfgs) if all(p[k] == v for k, v in fixed.items())]
    gg = g[g.cfg.isin(ok)]
    return cfgs[int(gg.loc[obj(gg).idxmax(), "cfg"])]


if __name__ == "__main__":
    main()
