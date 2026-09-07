# Troubleshooting

## `zsh: command not found: vcc` (installed, but not found)

PATH issue, not a failed install — **do not reinstall**.

```bash
ls ~/.local/bin/vcc && echo "installed — PATH problem"
uv tool update-shell        # adds ~/.local/bin to your shell profile
rehash                      # zsh (bash: hash -r), or open a new terminal
vcc --version
```

If `~/.local/bin/vcc` is missing, the install landed in the wrong environment (e.g. `uv pip
install -e .` into a project `.venv`). Reinstall with `uv tool install …` (see install.md).

## `error: Failed to spawn: vcc`

The project isn't installed into the active environment. Use `uv tool install …` (global), not
`uv run --no-sync` (which skips installing).

## Git permission error during install

```
Permission denied (publickey). / Could not read from remote repository.
```
You're installing from a git URL you don't have access to. Install from PyPI instead —
`uv tool install vcc-cli` — which needs no repository access at all (see install.md).

## `vcc: Not logged in`

- Missing endpoint — `--endpoint`/`VCC_ENDPOINT` must match where you generated the key (dev vs
  prod). Re-run `whoami` with the right endpoint.
- Token not stored (no usable keychain) — set `VCC_TOKEN` as a fallback, or re-run `login`.

## `Your API token is invalid or has been revoked`

Generate a fresh key on the Credentials page and log in again. Generating a new key **revokes
the old one** — that's why another machine may suddenly stop working.

## macOS keychain prompt

First `login`/`whoami` may pop a keychain-access dialog. Click **Always Allow**.

## `403` / "not available until the challenge opens"

Submissions are gated until the Challenge opens, and registration must be approved. Check
`vcc whoami`, which names the specific blocker. This is a server-side gate — waiting or
completing registration is the only fix, and pointing the CLI somewhere else is not.

## `N cell(s) have a per-cell count total above the maximum of 1,000,000`

`prep` refused the file locally. No single cell may total more than **1,000,000 counts summed
across all genes**; scoring applies the same cap, so an over-cap cell fails either way. The
message names how many cells are over and the worst cell's row and total.

The cause is almost always the prediction's **library size**, not its cell count or gene count:
nothing about such a matrix looks wrong — values are integral, finite and in range — so a model
emitting a plausible-but-too-large per-cell total crosses it silently. Fix it by scaling the
prediction down, e.g. renormalizing each cell to a library size at or below the cap. Dropping
cells does not help and will break the per-perturbation cell-count check.

⚠️ **`--max-counts-per-cell -1` is almost never the right answer.** It disables only the *local*
check; the scoring job applies the cap regardless, so the submission still fails — just after a
multi-gigabyte upload and a scoring run instead of in seconds. Use it only to confirm that this
specific check is what's firing, never to get a submission through.

## Submission shows `failed`

`vcc status <entry-id>` prints the error + timestamp; `submit --wait`/`status` exit non-zero.

**Treat a scoring failure as a real failure and report it** with the entry id and the error text.
The scoring pipeline runs full-size submissions across all three contexts through to `published`,
so do not wave a failure off as a known backend issue. The most
common genuine causes are content problems the scorer catches that a `.vcc` container check cannot
(wrong per-perturbation cell counts, submitted control cells, log-normalized instead of raw
counts); `vcc prep --perts pert_counts.csv` catches all of those locally before an upload.

## HTTP 409 "already has a submission in progress"

Your team already has a non-terminal submission. Wait for it (`vcc status <id> --wait`) or let
the server's staleness window clear a crashed one. This guard is intentional.

## The skill itself seems stale or missing after installing the CLI

Run `vcc skill install` (re-run it after any `vcc` upgrade), then **restart the agent session**
(or `/reload`) so it picks up the refreshed `~/.claude/skills/vcc/` files.
