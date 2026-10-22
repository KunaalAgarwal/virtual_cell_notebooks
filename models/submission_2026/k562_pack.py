"""Ship every row of the Replogle 2022 K562 genome-wide screen in the repo, without figshare.

The final phase (contexts D, E, F) scores a perturbation list that is only published inside the
controls bundle, so the build cannot know in advance which K562 rows it needs. The full bulk file
(K562_gwps_normalized_bulk_01.h5ad, figshare article 20029387, file 35773217; 375 MB) is too big
for git, and figshare has failed to serve it from Colab (HTTP 500 pages saved as the .h5ad).

`pack` stores all 11,258 rows as int16 (z-score x 1000, after the same clip to +-10 and NaN -> 0
that xline_benchmark.load_k562 applies on load) in byte-shuffled, zlib-compressed .npz parts of
under 50 MB each. `unpack` rebuilds a full .h5ad with the original obs, var and row order, then
overwrites the rows that are also in k562_gwps_bulk_subset.h5ad (the validation targets, the H1
knockdowns and the controls) with their exact float32 values. So the validation-panel build is
unchanged, and any other target's profile differs from the original by at most 0.0005 z (the
median |z| in the screen is 0.06).

    python k562_pack.py pack --full K562_gwps_normalized_bulk_01.h5ad --out data/k562_gwps_all
    python k562_pack.py unpack --parts data/k562_gwps_all --subset data/k562_gwps_bulk_subset.h5ad \
        --out K562_gwps_rebuilt.h5ad
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SCALE = 1000.0
CLIP = 10.0


def _cols(df, prefix):
    return {f"{prefix}{c}": (df[c].to_numpy().astype(str) if df[c].dtype.name in ("category", "object")
                             else df[c].to_numpy()) for c in df.columns}


def pack(full, out, rows_per_part):
    a = ad.read_h5ad(full)
    X = np.asarray(a.X, dtype=np.float32)
    X = np.clip(np.where(np.isfinite(X), X, 0), -CLIP, CLIP)
    q = np.round(X * SCALE).astype("<i2")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    n = q.shape[0]
    starts = list(range(0, n, rows_per_part))
    for k, s in enumerate(starts):
        b = q[s:s + rows_per_part].view(np.uint8).reshape(-1, q.shape[1], 2)
        np.savez_compressed(out / f"part_{k}.npz", lo=b[..., 0], hi=b[..., 1], start=s)
    np.savez_compressed(out / "meta.npz", obs_index=a.obs.index.to_numpy().astype(str),
                        var_index=a.var.index.to_numpy().astype(str), n_parts=len(starts),
                        scale=SCALE, shape=np.array(q.shape), **_cols(a.obs, "obs__"), **_cols(a.var, "var__"))
    for p in sorted(out.glob("*.npz")):
        print(f"{p.name}: {p.stat().st_size / 1e6:.1f} MB")


def unpack(parts, subset, out):
    parts = Path(parts)
    m = np.load(parts / "meta.npz", allow_pickle=False)
    X = np.zeros(tuple(m["shape"]), dtype=np.float32)
    for k in range(int(m["n_parts"])):
        z = np.load(parts / f"part_{k}.npz", allow_pickle=False)
        b = np.stack([z["lo"], z["hi"]], axis=-1)
        s = int(z["start"])
        X[s:s + b.shape[0]] = np.ascontiguousarray(b).view("<i2")[..., 0].astype(np.float32) / float(m["scale"])
    obs = pd.DataFrame({k[5:]: m[k] for k in m.files if k.startswith("obs__")}, index=m["obs_index"])
    var = pd.DataFrame({k[5:]: m[k] for k in m.files if k.startswith("var__")}, index=m["var_index"])
    var["gene_name"] = var["gene_name"].astype("category")
    exact = 0
    if subset and Path(subset).exists():
        sub = ad.read_h5ad(subset)
        assert sub.var["gene_name"].astype(str).tolist() == var["gene_name"].astype(str).tolist()
        pos = pd.Series(np.arange(len(obs)), index=obs.index).reindex(sub.obs.index)
        assert pos.notna().all(), "subset rows missing from the packed screen"
        X[pos.to_numpy().astype(int)] = np.asarray(sub.X, dtype=np.float32)
        exact = sub.n_obs
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.uns["source"] = ("Replogle et al. 2022, K562_gwps_normalized_bulk_01.h5ad (figshare 20029387 / file "
                       f"35773217), rebuilt by k562_pack.py; {exact} rows exact, the rest to 1e-3 z")
    a.write_h5ad(out)
    print(f"wrote {out}: {a.n_obs} rows x {a.n_vars} genes ({exact} exact rows from the subset)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack")
    p.add_argument("--full", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--rows-per-part", type=int, default=3000)
    u = sub.add_parser("unpack")
    u.add_argument("--parts", required=True)
    u.add_argument("--subset", default=str(Path(__file__).with_name("data") / "k562_gwps_bulk_subset.h5ad"))
    u.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "pack":
        pack(a.full, a.out, a.rows_per_part)
    else:
        unpack(a.parts, a.subset, a.out)


if __name__ == "__main__":
    main()
