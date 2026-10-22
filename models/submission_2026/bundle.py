"""Read a VCC 2026 controls bundle without assuming which phase it is for.

The validation bundle (`vcc datasets download controls`) covers contexts A, B, C and the
final-phase bundle (`controls-final`) covers D, E, F with its own perturbation list. Both ship
gene_names.csv, pert_counts.csv, manifest.json and one control .h5ad per context. Everything the
build needs is read from the unpacked bundle folder:

  * contexts: manifest.json "contexts" (the same source `vcc prep` uses), else the `context`
    labels found in the control .h5ad files;
  * targets: pert_counts.csv "target_gene", per context when the file has a `context` column;
  * the control file of each context: context_<X>.h5ad, else the .h5ad whose obs `context` is X.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

CTRL = "non-targeting"


def read_manifest(data):
    p = Path(data) / "manifest.json"
    if not p.is_file():
        return None
    m = json.loads(p.read_text(encoding="utf-8"))
    return m if isinstance(m, dict) else None


def _h5ad_contexts(path):
    """The `context` labels in an h5ad's obs, read straight from HDF5 (no matrix load)."""
    with h5py.File(path, "r") as f:
        if "obs" not in f or "context" not in f["obs"]:
            return set()
        g = f["obs"]["context"]
        if isinstance(g, h5py.Group):  # categorical
            cats = g["categories"][()]
            used = np.unique(g["codes"][()])
            vals = cats[used[used >= 0]]
        else:
            vals = np.unique(g[()])
    return {v.decode() if isinstance(v, bytes) else str(v) for v in vals}


def control_files(data):
    """{context: path} for every control .h5ad in the bundle folder."""
    out = {}
    for p in sorted(Path(data).glob("*.h5ad")):
        labels = _h5ad_contexts(p)
        stem = p.stem
        if stem.startswith("context_") and (not labels or labels == {stem[len("context_"):]}):
            out[stem[len("context_"):]] = p
        elif len(labels) == 1:
            out.setdefault(next(iter(labels)), p)
    return out


def contexts(data, flag=None):
    """The contexts to predict: --contexts if given, else manifest.json, else the control files."""
    if flag:
        return [c.strip() for c in flag.split(",") if c.strip()]
    m = read_manifest(data)
    if m and isinstance(m.get("contexts"), list) and m["contexts"]:
        return [str(c) for c in m["contexts"]]
    found = sorted(control_files(data))
    if not found:
        raise SystemExit(f"no manifest.json and no control .h5ad files in {data}: is this the unzipped bundle?")
    return found


def control_file(data, ctx):
    p = Path(data) / f"context_{ctx}.h5ad"
    if p.is_file():
        return p
    files = control_files(data)
    if ctx not in files:
        raise SystemExit(f"no control .h5ad for context {ctx} in {data} (found contexts {sorted(files)})")
    return files[ctx]


def targets(data, ctxs):
    """{context: [target genes]} from pert_counts.csv (controls dropped, file order kept)."""
    df = pd.read_csv(Path(data) / "pert_counts.csv")
    df = df[df["target_gene"].notna()]
    df = df[df["target_gene"].astype(str).str.strip() != CTRL]
    if "context" in df.columns:
        out = {c: list(dict.fromkeys(df.loc[df["context"].astype(str) == c, "target_gene"].astype(str).str.strip()))
               for c in ctxs}
        missing = [c for c, t in out.items() if not t]
        if missing:
            raise SystemExit(f"pert_counts.csv lists no targets for context(s) {missing}")
        return out
    t = list(dict.fromkeys(df["target_gene"].astype(str).str.strip()))
    return {c: t for c in ctxs}


def describe(data, flag=None):
    """One line per fact, printed by the notebook so a wrong bundle is caught before the build."""
    m = read_manifest(data) or {}
    ctxs = contexts(data, flag)
    tg = targets(data, ctxs)
    lines = [f"panel {m.get('panel_id', '?')} (partition {m.get('partition', '?')}), contexts {', '.join(ctxs)}"]
    for c in ctxs:
        lines.append(f"  context {c}: {len(tg[c])} targets, controls {control_file(data, c).name}")
    return "\n".join(lines)


if __name__ == "__main__":
    import sys
    print(describe(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None))
