# tnbc-genomics-agent

An AI-powered bioinformatics pipeline for analyzing RTK and non-receptor tyrosine kinase (nRTK) redundancy in Triple-Negative Breast Cancer (TNBC) genomic data. The pipeline parses patient VCF files, maps variants onto a curated kinase database, scores pathway redundancy, identifies bypass risk, and synthesizes findings via a Claude AI reasoning agent.

---

## Table of contents

- [Background](#background)
- [Project structure](#project-structure)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Pipeline steps](#pipeline-steps)
- [Modules](#modules)
- [Configuration](#configuration)
- [Running tests](#running-tests)
- [Adding your own VCF](#adding-your-own-vcf)
- [Extending the kinase database](#extending-the-kinase-database)
- [Dashboard](#dashboard)
- [License](#license)

---

## Background

Triple-Negative Breast Cancer (TNBC) lacks the three receptors (ER, PR, HER2) that standard therapies target, making it one of the hardest breast cancer subtypes to treat. A key challenge is **pathway redundancy**: when one RTK or nRTK is blocked by a drug, cancer cells often compensate by activating another kinase that feeds the same downstream pathway (typically MAPK/ERK or PI3K/AKT).

This pipeline automates the detection and scoring of that redundancy from a patient's VCF file, and uses Claude to synthesize the findings into a clinician-readable report.

---

## Project structure

```
tnbc-genomics-agent/
├── pipeline.py                  # Main orchestrator — run this
├── config.py                    # Centralised settings & thresholds
├── setup.py                     # Package installation
├── requirements.txt             # Python dependencies
│
├── modules/
│   ├── __init__.py
│   ├── kinase_db.py             # Curated RTK/nRTK database (10 kinases)
│   ├── vcf_parser.py            # VCF parsing & variant extraction
│   ├── redundancy_analyzer.py  # Pathway redundancy scoring
│   └── ai_agent.py             # Claude API reasoning layer
│
├── data/
│   ├── sample/
│   │   └── patient_1.vcf       # Demo VCF with 13 variants
│   └── real/                   # Place real patient VCFs here (gitignored)
│
├── reports/                     # JSON outputs land here (gitignored)
│   └── .gitkeep
│
├── tests/
│   ├── __init__.py
│   ├── test_kinase_db.py
│   ├── test_vcf_parser.py
│   └── test_redundancy_analyzer.py
│
└── notebooks/                   # Jupyter notebooks (optional)
```

---

## Installation

**Prerequisites:** Python 3.10+, pip

```bash
# 1. Clone the repo
git clone https://github.com/yourusername/tnbc-genomics-agent.git
cd tnbc-genomics-agent

# 2. Create a virtual environment
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. (Optional) Install as a package
pip install -e .

# 5. Set your Anthropic API key (required for AI analysis step)
export ANTHROPIC_API_KEY="sk-ant-..."
```

---

## Quick start

Run the full pipeline on the included sample VCF:

```bash
python pipeline.py --vcf data/sample/patient_1.vcf --out reports/patient_1_report.json
```

Expected output:

```
══════════════════════════════════════════════════════════════════════
  TNBC RTK/nRTK Redundancy Pipeline  |  2026-06-22 21:09:01
══════════════════════════════════════════════════════════════════════

Step 1 — VCF Parsing
  Total variants parsed   : 13
  RTK/nRTK variants found : 10
  Kinase genes detected   : ABL1, ALK, EGFR, ERBB2, FGFR1, FGFR2, FGFR3, MET, SRC

Step 3 — Pathway Redundancy Analysis
  MAPK/ERK      ████████████████████  0.865
  PI3K/AKT      ████████████████      0.848
  ...

Step 5 — Combination Therapy Suggestions
  ➤ EGFR + MET  [MAPK/ERK]  Evidence: STRONG
     Drugs: Erlotinib / Gefitinib  +  Crizotinib / Cabozantinib
  ...
```

---

## Pipeline steps

| Step | What it does |
|------|-------------|
| 1 | Parse VCF — extract all variants, filter to RTK/nRTK genes |
| 2 | Variant detail — annotate with kinase type, TNBC relevance, effect severity |
| 3 | Pathway redundancy — score how many independent genes feed each pathway |
| 4 | Bypass risk — per-gene estimate of monotherapy failure risk |
| 5 | Combination suggestions — rational drug pairs to block redundant pathways |
| 6 | High-severity variants — surface stop-gained, frameshift, missense hits |
| 7 | AI synthesis — Claude narrative report (requires `ANTHROPIC_API_KEY`) |

---

## Modules

### `modules/kinase_db.py`

A curated dictionary of 10 kinases (7 RTKs, 3 nRTKs) with:
- Full name, family, chromosomal location
- TNBC overexpression frequency
- Downstream pathway memberships
- FDA-approved and investigational inhibitors
- Known resistance mechanisms and redundant bypass receptors

To add a new gene, append an entry to `KINASE_DATABASE` following the existing schema.

### `modules/vcf_parser.py`

Parses standard VCF 4.x files. Requires `GENE` and `EFFECT` in the `INFO` field for full annotation; falls back to a chromosomal locus map for common kinase positions.

Key class: `VCFParser(vcf_path).parse()` → `.variants`, `.rtk_nrtk_variants`, `.get_summary()`

### `modules/redundancy_analyzer.py`

Computes:
- **Pathway redundancy score** (0–1) per pathway using a sigmoid scaling of gene count
- **Bypass risk per target** — how many co-affected genes could compensate if a given gene is inhibited
- **Combination therapy suggestions** — ranked pairs of gene targets with supporting drug options

Key class: `RedundancyAnalyzer(rtk_nrtk_variants).get_full_report()`

### `modules/ai_agent.py`

Builds the structured prompt and system context sent to `claude-sonnet-4-6`. The prompt includes the full VCF summary and redundancy report as JSON, and requests a six-section clinical narrative. Designed to be called both from Python and from the browser dashboard via the Anthropic API.

---

## Configuration

All thresholds are in `config.py` and can be overridden with environment variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `ANTHROPIC_API_KEY` | — | Required for AI analysis |
| `TNBC_AI_MODEL` | `claude-sonnet-4-6` | Claude model to use |
| `TNBC_MIN_AF` | `0.05` | Minimum allele frequency filter |
| `TNBC_MIN_DEPTH` | `20` | Minimum read depth filter |
| `TNBC_OUTPUT_DIR` | `reports` | Directory for JSON reports |

Example:

```bash
TNBC_MIN_AF=0.10 TNBC_MIN_DEPTH=50 python pipeline.py --vcf my_patient.vcf
```

---

## Running tests

With pytest installed:

```bash
pytest tests/ -v
```

Without pytest (built-in):

```bash
python -m unittest discover tests/ -v
```

The test suite covers:
- VCF parsing correctness (variant counts, gene detection, field validation)
- Redundancy scoring range and monotonicity
- Kinase database schema integrity
- Bypass risk ordering and self-reference exclusion

---

## Adding your own VCF

1. Place your VCF file in `data/real/` (this directory is gitignored — patient data never leaves your machine).
2. Ensure your VCF `INFO` column contains `GENE=<gene_symbol>` and `EFFECT=<vep_effect>` fields, or annotate with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) first.
3. Run:

```bash
python pipeline.py --vcf data/real/your_patient.vcf --out reports/your_patient_report.json
```

---

## Extending the kinase database

Open `modules/kinase_db.py` and add a new entry to `KINASE_DATABASE`:

```python
"NTRK1": {
    "full_name": "Neurotrophic Receptor Tyrosine Kinase 1",
    "type": "RTK",
    "family": "NTRK",
    "chromosome": "1q23.1",
    "tnbc_relevance": "moderate",
    "overexpression_freq_tnbc": 0.10,
    "pathways": ["MAPK/ERK", "PI3K/AKT"],
    "known_inhibitors": [
        {"drug": "Larotrectinib", "type": "small_molecule", "fda_approved": True, "trial_phase": "approved"},
        {"drug": "Entrectinib",   "type": "small_molecule", "fda_approved": True, "trial_phase": "approved"},
    ],
    "resistance_mechanisms": ["NTRK point mutations", "bypass via EGFR"],
    "redundant_pathways": ["EGFR", "MET"],
    "clinical_notes": "NTRK fusions rare in TNBC but highly actionable when present.",
},
```

The parser, redundancy analyzer, and AI agent will automatically incorporate it on the next run.

---

## Dashboard

An interactive browser dashboard is available as a React/HTML artifact (see the project conversation). It visualises:
- Variant table with allele frequencies
- Pathway redundancy bar chart
- Per-gene bypass risk scores
- Combination therapy suggestions
- Live AI analysis tab (calls the Anthropic API client-side)

---

## License

MIT License — see `LICENSE` for details.

> **Note on patient data:** This repository is designed for research use. Real patient VCF files are gitignored by default. Never commit identifiable genomic data to a public repository.
