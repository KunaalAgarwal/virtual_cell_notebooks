# Installing / updating the `vcc` CLI

The distribution on PyPI is **`vcc-cli`**; the command it installs is **`vcc`**. (The bare name
`vcc` on PyPI belongs to an unrelated project.)

**Prerequisites:** Python **3.11+** (or `uv`, which can fetch one), and `uv` **or** `pipx`.

## Install

```bash
uv tool install vcc-cli        # or: pipx install vcc-cli
vcc --version
```

## `command not found: vcc` right after installing

Almost always PATH, not a failed install — **do not reinstall**. `uv tool` links a shim into
`~/.local/bin`:

```bash
uv tool update-shell     # adds ~/.local/bin to your shell profile
rehash                   # zsh (bash: hash -r), or open a new terminal
vcc --version
```

## Updating

```bash
uv tool upgrade vcc-cli        # or: pipx upgrade vcc-cli
```

Re-run **`vcc skill install`** after upgrading so this skill matches the installed CLI version.

## Uninstalling

```bash
uv tool uninstall vcc-cli     # or: pipx uninstall vcc-cli   (the PyPI name, not the command)
```

## Installing from a source checkout

Someone working on the CLI itself may install it from a local checkout
(`uv tool install --editable <path>`). That produces a **dev build**, whose endpoint is
configurable; a build installed from PyPI locks to production. `vcc version` prints which one is
in use.
