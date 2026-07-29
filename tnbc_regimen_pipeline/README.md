# TNBC Regimen Pipeline

Kinase-alteration-driven combination therapy scoring for TNBC (TCGA-BRCA
cohort), extended with literature-mining discovery of novel drug
candidates for under-covered genes.

## Install

```bash
pip install -e ".[dev]"
```

## Structure

```
tnbc_regimen_pipeline/
├── tnbc_regimen_pipeline/
│   ├── discovery/           # PubMed search, INN-suffix extraction,
│   │                          mandatory DGIdb confirmation, caching,
│   │                          parallelization, provenance
│   ├── cohort/               # TNBC patient ID, MAF reading, cohort-wide
│   │                          MDCOE/HCOS scoring, frequency aggregation
│   ├── pipeline/              # orchestration, diffing, exports, visual
│   │                          summary, reproducibility stamping
│   ├── config/                # PipelineConfig + YAML/JSON loaders
│   ├── cli/                   # command-line entrypoint
│   └── utils/                 # logging setup, schema validation
├── tests/                      # pytest suite
├── pyproject.toml
├── setup.cfg
└── LICENSE                     # MIT -- placeholder; confirm this is the
                                   right choice for your institution
```

## Design discipline: this package doesn't own MDCOE/HCOS scoring

`resolve_tp53_fn`, `DrugGraph`, `SynergyNet`, `hcos_fn`, and `mdcoe_fn`
are NOT part of this package -- they're your existing, separately
validated scoring logic. Everything here (`cohort.run_cohort_wide_analysis`,
`pipeline.run_full_pipeline`, the CLI) takes them as parameters/plugin
callables rather than reimplementing them. If you extend this pipeline,
keep that boundary: add a new wrapping layer rather than editing
`discovery/` or `cohort/` directly, and re-run `pytest tests/` before
trusting anything built on top of a change to either.

## Running via Python

```python
from tnbc_regimen_pipeline.config import load_default_config
from tnbc_regimen_pipeline.pipeline import run_full_pipeline

config = load_default_config()
config.output_dir = "./outputs"

result = run_full_pipeline(
    genes=[...], curated_gene_drugs={...}, real_gene_drugs={...},
    patient_barcodes=[...], maf_glob_pattern="data/*.maf.gz", gene_panel=[...],
    resolve_tp53_fn=..., drug_graph_cls=..., synergy_net_cls=...,
    hcos_fn=..., mdcoe_fn=..., config=config,
)
```

## Running via CLI

Three ways to drive the pipeline from a terminal, in increasing order of
convenience:

**`discover`** — literature discovery only, no cohort scoring:

```bash
tnbc-pipeline discover \
  --genes EGFR,PTEN,TYK2 \
  --curated-json curated.json --real-json real.json \
  --output-dir outputs/
```

**`full-pipeline`** — full run, one flag per input:

```bash
tnbc-pipeline full-pipeline \
  --genes EGFR,PTEN,TYK2 \
  --curated-json curated.json --real-json real.json \
  --patients-file patients.txt --maf-glob "data/*.maf.gz" \
  --gene-panel EGFR,PTEN,TYK2 \
  --plugin my_mdcoe_plugin:get_pipeline_components \
  --output-dir outputs/
```

**`run`** — full run from a single YAML job-config file:

```bash
tnbc-pipeline run job_config.yaml
```

```yaml
# job_config.yaml
genes: [EGFR, PTEN, TYK2]
gene_panel: [EGFR, PTEN, TYK2]
curated_gene_drugs_json: curated.json   # or inline: curated_gene_drugs: {...}
real_gene_drugs_json: real.json
patients_file: patients.txt              # or inline: patient_barcodes: [...]
maf_glob_pattern: "data/*.maf.gz"
plugin: my_mdcoe_plugin:get_pipeline_components
output_dir: outputs/
```

All three need `--plugin`/`plugin:` rather than accepting `DrugGraph`,
`hcos_fn`, etc. directly, for the reason explained above: those are
Python objects, and neither a CLI flag nor a YAML file can hold one.
`cli/main.py`'s module docstring has the exact plugin contract and a
worked example.

## Testing

```bash
pytest tests/
```

**Honest note on how this suite was actually verified:** it was built and
run in a sandbox with no network access and no `pytest` installed (no
network to add it either). Every test file here uses standard pytest
syntax (`@pytest.fixture`, plain `assert`, `monkeypatch`, `tmp_path`) and
will run correctly under a real `pytest` install -- but in that sandbox,
verification was done with a small local harness
(`_manual_test_harness.py`, not part of the package) that implements just
enough of pytest's fixture-resolution mechanism to execute these exact
files. All 21 tests passed under that harness. **Run `pytest tests/` for
real in your environment before trusting this further** -- the harness is
a stand-in, not a substitute for the real tool.

What IS genuinely real in these tests, not mocked: MAF file parsing is
tested against an actual gzipped file written to `tmp_path`, and the CLI
was separately verified end-to-end (both `discover` and `full-pipeline`
subcommands, via `click.testing.CliRunner`) against a real gzipped MAF
file and a real plugin module. PubMed/DGIdb network calls are mocked
throughout, since there is no network access to test against the live
services.

## Known open items

- `resolve_tp53_fn`'s default assumption (every patient's TP53 alteration
  is `Nonsense_Mutation`) is only confirmed for one previously-inspected
  patient (TCGA-AO-A128) — verify each real patient's actual
  classification before trusting cohort-wide TP53-arm results.
- Schema validation (`utils/validation.py`) is hand-written, not
  `pydantic` — not installed/available when this was built. Fine for the
  current shallow `Dict[str, List[str]]` shape.
- The gene→drug→regimen flow diagram (`generate_gene_drug_regimen_flow_svg`)
  is a hand-rolled SVG, not a true Sankey — `plotly` wasn't available
  either. Edge widths aren't unified across both stages the way a real
  Sankey's would be.
- The `LICENSE` file is MIT as a placeholder — confirm this is actually
  the license you want before publishing, especially if this touches
  institutional or grant-funded IP policy.

## Data note

TCGA-BRCA clinical and MAF files are broadly open-access, but
patient-level derived outputs (regimen tables, discovery results) are
gitignored by default (`outputs/`) rather than judged case-by-case.
