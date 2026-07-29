# Running This Repo in Fish Shell (Ubuntu)

Fish differs from bash in three ways that matter here: `source venv/bin/activate` needs the
`.fish` variant, environment variables are set with `set -x` instead of `export` or inline
`VAR=val cmd`, and there's no direct equivalent of `set -euo pipefail`. Commands below are
fish-native, not bash pasted in and hoped for.

## 0. Confirm fish is installed

```fish
fish --version
```

If that fails: `sudo apt update; and sudo apt install fish`

## 1. One virtual environment for the whole repo

```fish
cd tnbc-project
python3 -m venv .venv
source .venv/bin/activate.fish
pip install -r requirements.txt
```

`activate.fish` is created automatically alongside the regular `activate` script by
`python3 -m venv` — no separate install step needed. Confirm it took effect:

```fish
which python3   # should point inside .venv/bin/
```

Deactivate later with `deactivate` (same command name in fish and bash).

## 2. `src/` — CTS/HCOS/DepMap cohort-scale project

```fish
python3 src/survival/patient_level_survival_model.py   # runs the built-in smoke test
```

Once real data loaders exist (see each `src/*/README.md` for what's still a placeholder):

```fish
python3 scripts/run_track_a_survival.py
```

Local data setup — use the fish-native version, not the `.sh` one:

```fish
# first edit the four SOURCE_* paths near the top of the file
fish scripts/setup_local_data.fish
```

(`scripts/setup_local_data.sh` still works too if you'd rather run it via `bash
scripts/setup_local_data.sh` — fish can invoke a bash script fine, it just can't *source* one
into your current fish session, which is why this repo has a native `.fish` copy instead.)

## 3. `tnbc-genomics-agent/` — VCF pipeline

```fish
cd tnbc-genomics-agent
pip install -r requirements.txt
```

**Setting the Claude API key** (needed for Step 7's AI narrative; the rest of the pipeline
runs fine without it):

```fish
set -x ANTHROPIC_API_KEY sk-ant-...
```

`-x` exports it to child processes, matching what `export` does in bash. This only lasts for
the current fish session — add it to `~/.config/fish/config.fish` if you want it permanent:

```fish
echo 'set -x ANTHROPIC_API_KEY sk-ant-...' >> ~/.config/fish/config.fish
```

Run the pipeline:

```fish
python3 pipeline.py --vcf data/sample/patient_1.vcf --out reports/patient_1_report.json
```

**Config overrides** — bash lets you prefix a command with inline `VAR=val`; fish doesn't
support that syntax at all. Use `env`, or `set -x` for the duration of one command with `and`:

```fish
env TNBC_MIN_AF=0.10 TNBC_MIN_DEPTH=50 python3 pipeline.py --vcf my_patient.vcf
```

**Running the verified test suite:**

```fish
python3 _manual_test_harness.py   # verified 31/31 passing — see INTEGRATION_NOTES.md
```

If you have real network access and want to install pytest properly instead of relying on the
harness:

```fish
pip install pytest
python3 -m pytest tests/ -v
```

## 4. `mc-ore-prototype/` — MC-ORE/Hybrid-CORE/CL-MODE prototype

```fish
cd mc-ore-prototype
pip install pandas jupyter
jupyter nbconvert --to notebook --execute tnbc_combo_pipeline.ipynb
```

Same `set -x ANTHROPIC_API_KEY sk-ant-...` applies here if you want live rationale generation
instead of the offline template fallback (Phase 3 of the notebook).

## Common fish gotchas coming from bash

| bash | fish |
|---|---|
| `export VAR=value` | `set -x VAR value` |
| `VAR=value cmd` (inline, one-shot) | `env VAR=value cmd` |
| `source .venv/bin/activate` | `source .venv/bin/activate.fish` |
| `$()` command substitution | `()` — same idea, no `$` prefix |
| `&&` / `\|\|` | `and` / `or` |
| `set -euo pipefail` | no direct equivalent — handle errors per-command with `; or` |
| `~/.bashrc` | `~/.config/fish/config.fish` |
