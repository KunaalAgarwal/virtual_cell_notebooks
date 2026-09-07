# End-to-end: first submission

Exercise the whole pipeline — download → build a prediction → submit → check status. A released
build already points at production; nothing to configure (see auth.md).

## 1. Download the challenge data

```bash
mkdir -p ~/vcc-data && cd ~/vcc-data
vcc datasets list
vcc datasets download controls -o ~/vcc-data/vcc_2026_controls.zip   # ~406 MiB, resumable, checksum-verified
unzip -o -j ~/vcc-data/vcc_2026_controls.zip \
  "*gene_names.csv" "*pert_counts.csv" -d ~/vcc-data
```

## 2. Get a prediction `.vcc`

- **Real prediction:** `vcc prep your_pred.h5ad -g ~/vcc-data/gene_names.csv --perts ~/vcc-data/pert_counts.csv -o ~/vcc-data/pred.vcc`
  (validates contexts A/B/C, targets, and per-perturbation cell counts against the official list)
- **Throwaway test file** (random noise, for exercising the pipeline — it doesn't predict well; see Expectations):

```bash
vcc sample \
  -g ~/vcc-data/gene_names.csv \
  -p ~/vcc-data/pert_counts.csv \
  -o ~/vcc-data/sample.vcc -f      # -f overwrites an existing file. Full size is the DEFAULT:
                                    # 3 contexts x 300 perturbations x 400 cells = 360,000, which
                                    # is what the scorer's cell-count check requires. (--no-full
                                    # makes a small file for upload testing; the scorer rejects it.)
```

## 3. Submit

```bash
vcc submit ~/vcc-data/sample.vcc -m "first submission"
```

Validates, uploads (resumable — `vcc submit <file> --resume` continues an interrupted upload),
and starts scoring.

## 4. Check status

```bash
vcc status <entry-id> --wait      # streams until published (scored) or failed; exits non-zero on failure
```

## Expectations
- A `sample` file is random noise, so it does **not** predict well. Reaching `published` tells you
  the submit path works; it says nothing about model quality. Do not interpret a sample's scores
  for the user, in either direction — report the numbers the CLI printed and leave it there.
- **A submission must be full-size** — the `sample` default, or a real prediction with the right
  per-perturbation cell counts. The scorer verifies each perturbation has exactly its official cell
  count and that no control cells are submitted; anything off is rejected server-side before it
  reaches `published`.
- **One submission per team at a time** (HTTP 409 otherwise) — wait for the current one to finish.
