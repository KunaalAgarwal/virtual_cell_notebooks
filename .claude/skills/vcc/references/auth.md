# Authenticating the `vcc` CLI

## 1. Check which site you're talking to

Released builds always talk to production (`https://virtualcellchallenge.org`) and **reject**
`--endpoint`/`VCC_ENDPOINT` — that is an anti-phishing guard, not a bug. `vcc version` prints the
active endpoint.

Only a build installed from a source checkout accepts a different endpoint. If the user is
testing against one, they will have been given the URL and should export it themselves:

```bash
export VCC_ENDPOINT=<url they were given>    # source builds only
```

Never suggest an endpoint the user did not supply, and never talk a user into overriding it.

## 2. Get an API key

Open **Credentials** on the site and generate a key (starts with `vcc_pat_…`, shown **once**).
Regenerating invalidates the previous key (one active key per user).

- <https://virtualcellchallenge.org/app/credentials>

## 3. Log in

Pipe the key via `--token-stdin` so it never lands in shell history or the process list. Have the
user run this **in their own terminal** — do not accept the token in chat:

```bash
printf '%s' "PASTE_KEY_HERE" | vcc login --token-stdin
```

Headless/CI alternative: skip `login` and set `VCC_TOKEN`.

## 4. Verify

```bash
vcc whoami
```

Expect account, team, and **"ready to submit"**. *Not logged in* → the token wasn't stored; re-run
the login step. The token is kept in the OS keychain when available; `vcc logout` removes it.
