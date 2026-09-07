# `vcc` command reference

All commands support `--json` for machine-readable output. Global options work **either before
the subcommand or on the subcommand itself** — both positions are accepted:

```bash
vcc --json submit pred.vcc   # before the subcommand
vcc submit pred.vcc --json   # on the subcommand — also fine
```

## Commands

| Command | What it does |
|---|---|
| `vcc version` | Print version (`--json` for details). |
| `vcc login` / `logout` / `whoami` | Authenticate (PAT), clear the credential, show who you are. See auth.md. |
| `vcc datasets list` | Show which reference datasets you can download. |
| `vcc datasets download <id>` | Download a dataset (resumable, checksum-verified). `controls` is the control-cells + gene-list + `pert_counts.csv` bundle (`pert_counts.csv` is a one-column list of the 300 perturbations). |
| `vcc prep <input.h5ad> -g genes.csv --perts pert_counts.csv -o out.vcc` | Validate (genes, contexts A/B/C, targets, per-perturbation cell counts, no controls, raw counts) → slim → package a real prediction `.h5ad` into a submittable `.vcc`. |
| `vcc sample -g genes.csv -p pert_counts.csv -o out.vcc` | Generate a random-but-valid **test** `.vcc` (scores poorly by design, but the scorer accepts it: every perturbation, 400 cells each). `-p` is required. `--no-full` makes a small file the scorer **rejects** — upload testing only. |
| `vcc submit <file.vcc> -m "msg"` | Upload (resumable, checksum-verified) and start scoring. `--resume` continues an interrupted upload. |
| `vcc status <entry-id>` | Show a submission's status/scores. `--wait` blocks until terminal and exits non-zero on failure. |

## Common flows

```bash
# Download data (control cells + gene_names.csv + pert_counts.csv)
vcc datasets download controls -o ~/vcc-data/vcc_2026_controls.zip

# Prep a real prediction
vcc prep pred.h5ad -g gene_names.csv --perts pert_counts.csv -o pred.vcc

# Or generate a throwaway test file (scorable by default: ~300k cells, ~240 MB)
vcc sample -g gene_names.csv -p pert_counts.csv -o sample.vcc -f

# Submit and watch
vcc submit pred.vcc -m "run 1"
vcc status <entry-id> --wait
```

## Rules the CLI enforces (explain these, don't fight them)

- **One in-flight submission per team.** A second submission while one is still
  uploading/scoring returns **HTTP 409** — wait for the first (`vcc status <id> --wait`) or let
  the server's staleness window clear a crashed one.
- **Resumable, verified transfers.** Interrupted `submit` uploads and `datasets download`s
  resume where they left off; the final object's checksum is verified, so a corrupted transfer
  is never marked ready.
- **`sample` is noise.** A `vcc sample` file exists to exercise the pipeline, not to score.

## Flags at a glance

`vcc <command> --help` is the authoritative, always-current list (it prints each flag's
description and default). This summary is for quick reference; when in doubt, run `--help`.

**Global (before the subcommand; `--profile`/`--endpoint` also work on the subcommand):** `--json` · `--profile <name>` · `--endpoint <url>` (env: `VCC_PROFILE`, `VCC_ENDPOINT`) · `--version` · `-h/--help`

**`login`** — `--token-stdin` (preferred; read token from stdin) · `--token <t>` (visible in history) · `--endpoint <url>` (remembered for the profile) · `--store-plaintext` (0600 file if no keychain) · `--no-store` (validate only) · `--json`

**`whoami`** — `--endpoint <url>` · `--profile <name>` · `--json`

**`logout`** — `--profile <name>` · `--json`

**`prep <input.h5ad>`** — `-i/--input` (alt to positional) · `-g/--genes` (required) · `--perts <pert_counts.csv>` (official list; **required by default** — targets are checked against it; `--no-verify-targets` to skip) · `-o/--output` · `-p/--pert-col` · `-c/--celltype-col` · `-n/--ntc-name` · `-P/--output-pert-col` · `-C/--output-celltype-col` · `--context-col` · `--contexts A,B,C` (`""` disables) · `--verify-targets/--no-verify-targets` · `--check-cell-counts/--no-check-cell-counts` · `--cells-per-pert <n>` (default 400 — the panel's per-perturbation count; an `n_cells` column in `--perts` wins over it, `-1` drops the expectation) · `--require-counts/--no-require-counts` (raw counts; `--no-` restores log-normalization) · `--reject-controls/--allow-controls` · `-e/--encoding {32|64}` (default 32) · `--allow-discrete` (legacy, with `--no-require-counts`) · `--expected-gene-dim <n>` (default 18533; `-1` disables) · `--max-cell-dim <n>` (default 400000; `-1` disables) · `--max-counts-per-cell <n>` (default 1000000 — no cell may total more counts across genes; `-1` disables) · `--dry-run` · `-f/--force` · `--json`

**`sample`** — `-g/--genes` (required) · `-p/--perts <csv>` (**required**) · `-o/--output` · `--full`/`--no-full` (default `--full`: the official 400 cells per perturbation; `--no-full` is rejected by the scorer) · `--cells-per-pert <n>` (default 400 with `--full`, 5 with `--no-full`; an `n_cells` column in `-p` wins) · `--ntc-cells <n>` · `--genes-per-cell <n>` · `--context-col` · `--contexts A,B,C` (`""` disables) · `--max-cell-dim <n>` (default 400000; `-1` disables) · `--seed <n>` · `--h5ad` (raw .h5ad instead of .vcc) · `-f/--force` · `--json`

**`submit <file>`** — `-m/--model-name` (required) · `-d/--description` (≤2000 chars) · `-g/--genes` (if FILE is a raw .h5ad) · `--perts <pert_counts.csv>` (if FILE is a raw .h5ad; required by default when prepping) · `--verify-targets/--no-verify-targets` · `--check-cell-counts/--no-check-cell-counts` · `--wait` · `--poll-interval <s>` (default 5) · `--wait-timeout <s>` · `--resume [entry-id]` · `--skip-limit-check` · `-f/--force` · `--endpoint` · `--profile` · `--json`

**`status <entry-id>`** — `--wait` · `--poll-interval <s>` · `--wait-timeout <s>` · `--endpoint` · `--profile` · `--json`

**`datasets list`** — `--endpoint` · `--profile` · `--json`

**`datasets download <id>`** — `-o/--output <file>` · `-d/--dir <dir>` (env: `VCC_DATA_DIR`) · `-f/--force` (re-download) · `--endpoint` · `--profile` · `--json`

**`skill install`** — `--agent {auto|claude|codex|gemini|all}` (default `auto` = every agent present on the machine) · `--dir <path>` · `--force` · `--json` · (`skill uninstall` same minus `--force`; `skill path --agent <a>`)

## Getting help

`vcc --help`, or `vcc <command> --help` for any command's flags.
