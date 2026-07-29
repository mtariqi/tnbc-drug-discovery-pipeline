# src/predictive_models

GNN-based predictive components intended to eventually replace two heuristic
placeholders elsewhere in this pipeline:

1. **Kinase response / redundancy prediction** (`models/kinase_cell.py`) -- a
   three-encoder architecture (cell state, perturbation, kinase network) producing
   response, pathway-activity, and redundancy predictions. Redundancy/pathway
   training is currently disabled (loss weight 0 in `configs/kinase_cell.yaml`)
   because no real label exists yet for either -- see `docs/limitations.md`'s note
   on missing combination-synergy data. Do not enable those losses against a
   fabricated label; wait for real data (see item 2, which may supply one).

2. **Drug-pair synergy prediction** (`models/synergy_gnn.py`) -- a molecular-graph
   GNN over two drugs' structures plus cell-line context, predicting real Bliss/
   Loewe/HSA/ZIP synergy scores. This is meant to eventually populate `PairCTS`'s
   and `TripletCTS`'s synergy term (see `src/scoring/`), which currently has no
   real synergy signal at all. **Implemented and verified against synthetic data
   only -- not yet trained on real data.** See "Running on real data" below.

Both models reuse the same leakage-safe splitting discipline as the rest of this
project (`data/splits.py`, `data/synergy_splits.py`) -- verified, not assumed;
see the tests in `tests/predictive_models/`.

## Status, plainly

- Architecture: implemented, unit-tested against synthetic data with known planted
  signal (same discipline as CTS's own synthetic validation).
- Real-data training: **not yet run**. `data/real_kinase_loader.py` and
  `data/real_drugcomb_loader.py`/`data/assemble_synergy_dataset.py` are written
  against expected real schemas (DepMap's confirmed column layout for the kinase
  side; DrugComb's *unverified* layout for the synergy side) but were never
  executed against your actual files in the session that built them. Treat your
  first real run as the actual integration test.
- **Known gap:** the synergy predictor only handles small-molecule drugs with a
  valid SMILES. Biologics (e.g. trastuzumab, which appears in this project's own
  focal-patient regimen in the manuscript) have no SMILES and are silently
  excluded by `assemble_synergy_dataset.py` -- check its printed "UNRESOLVED" list
  before assuming full drug coverage.

## Running on synthetic data (smoke test, no real data needed)

```bash
# from the repo root
python -m src.predictive_models.data.kinase_synthetic --output data/processed/kinase_synthetic.npz
python -m src.predictive_models.train_kinase --config src/predictive_models/configs/kinase_cell.yaml

python -m src.predictive_models.data.synergy_synthetic --output data/processed/synergy_synthetic.npz
python -m src.predictive_models.train_synergy --config src/predictive_models/configs/synergy.yaml
```

## Running on real data

**Fastest path (recommended first): PyTDC's pre-resolved DrugSyn dataset.**
Confirmed, current info (checked live this session): PyTDC ships a DrugComb-derived
dataset with SMILES already resolved (queried from PubChem by TDC's own
maintainers), skipping the PubChem name-matching step entirely.

```bash
pip install PyTDC rdkit
python -m src.predictive_models.data.load_tdc_drugsyn --inspect-only   # confirm real column names first
python -m src.predictive_models.data.load_tdc_drugsyn --tdc-name DrugComb --output data/processed/synergy_tdc.npz
# edit src/predictive_models/configs/synergy.yaml: data.path -> data/processed/synergy_tdc.npz
python -m src.predictive_models.train_synergy --config src/predictive_models/configs/synergy.yaml
```

This gets you training on **real** data fastest, but it's a curated NCI-60 subset
(129 drugs, 59 cell lines, ~297K rows) with a single synergy label replicated
across all four Bliss/Loewe/HSA/ZIP slots (see the module's docstring) --
appropriate for a first real run and a sanity check, not your final result.

**Full-scale path: the raw DrugComb portal file.** Confirmed real download (checked
live this session): `summary_v_1_5.csv` (1.4 GB) at
https://zenodo.org/records/15235991, with real columns `block_id, drug_row,
drug_col, cell_line_name, synergy_zip, synergy_bliss, synergy_loewe, synergy_hsa,
css/css_ri, study_name` (independently confirmed against a third-party source
this session, not just assumed). This is the ~1.4M-record full dataset, but
requires resolving drug names to SMILES yourself (`fetch_smiles.py`'s PubChem
fallback) since the raw file doesn't include structures.

```bash
pip install rdkit requests
# 1. Download from https://zenodo.org/records/15235991/files/summary_v_1_5.csv?download=1
# 2. Confirm columns match the above before trusting anything downstream:
python -m src.predictive_models.data.real_drugcomb_loader --summary-csv summary_v_1_5.csv --inspect-only
# 3. Assemble (resolves SMILES via PubChem, featurizes via RDKit, joins everything):
python -m src.predictive_models.data.assemble_synergy_dataset --summary-csv summary_v_1_5.csv --output data/processed/synergy_real.npz
# 4. python -m src.predictive_models.train_synergy --config src/predictive_models/configs/synergy.yaml
```

**Kinase model:** point `configs/kinase_cell.yaml`'s `data.path` at the output of
`src/predictive_models/data/real_kinase_loader.py`, which wires directly into
`src/data_loaders/depmap_multiomic_loader.py` and
`src/data_loaders/depmap_supplemental_loader.py` -- copy or symlink those two
files (or add `src/data_loaders` to your Python path) so the import resolves.

## Layout

```
src/predictive_models/
  encoders/       cell_state.py, perturbation.py, kinase_network.py, molecule.py
  models/         kinase_cell.py, synergy_gnn.py
  data/           splits.py, synergy_splits.py (leakage-safe splitting),
                  kinase_synthetic.py, synergy_synthetic.py (smoke-test generators),
                  real_kinase_loader.py, real_drugcomb_loader.py,
                  fetch_smiles.py, assemble_synergy_dataset.py (real-data path)
  configs/        kinase_cell.yaml, synergy.yaml
  train_kinase.py, train_synergy.py
tests/predictive_models/
  test_kinase_synthetic.py, test_synergy_splits.py
```
