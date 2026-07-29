# mc-ore-prototype/

A **research proposal + runnable Phase-1 prototype** extending `tnbc-genomics-agent` with a
GNN synergy-prediction layer, a biomedical knowledge graph, and multi-agent LLM reasoning.
This is proposed/prototype work, not a validated pipeline — treat everything here as a
hypothesis-generation scaffold, same caution level as the rest of this repo's "computational
prioritization, not clinical validation" framing.

## A naming note, stated plainly rather than silently picked

The source documents use **three different names for what appears to be the same evolving
architecture**: the proposal document (`docs/breast_cancer_hybrid_core_report.pdf`) calls it
**MC-ORE** (Multi-Agent Combination Oncology Reasoning Engine) in its main architecture
section, then switches to **Hybrid-CORE** for the wet-lab validation section ("referred to
here under their working translational-program names"), while the separate graphical abstract
(`docs/CL-MODE_graphical_abstract.pdf`) calls it **CL-MODE** (Closed-Loop Multi-Modal Oncology
Design Engine) with a somewhat different layer breakdown (KG++ / MM-SynergyNet / Multi-Agent
Oncology Reasoning Collective / Adaptive Experimentation Loop, vs. MC-ORE's L1–L5). These
line up closely enough to be the same underlying idea at different drafting stages, but they
have not been reconciled into one canonical name or architecture diagram. Pick one before this
goes further, rather than let three names for the same thing propagate into code and docs.

## What's actually runnable vs. what's still a proposal

| | Status |
|---|---|
| `tnbc_combo_pipeline.ipynb` | **Runnable now.** A genuine Phase-1 scaffold — placeholder `SynergyModel` (rule-based stand-in for a real GNN), a real curated toxicity-screening table, offline-safe rationale generation (falls back to a template if no `ANTHROPIC_API_KEY` is set). Verified in this session: re-running it from scratch reproduces `data/tnbc_ranked_drug_combinations.csv` byte-for-byte. |
| GNN synergy model (SynerGNet/DeepDDS-style) | **Not built.** `SynergyModel.predict_pair()` is a documented placeholder interface — see the notebook's Phase 2 section for exactly what training data (DrugCombDB, AZ-DREAM) and validation protocol would be needed. |
| Knowledge graph (KG++) | **Not built.** The notebook's `GENE_PATHWAYS`/`GENE_DRUGS`/`DRUG_TOXICITY` are flat Python dicts, not a graph structure — a real KG++ would need the drugs/targets/pathways/trials node-and-edge structure the proposal describes. |
| Multi-agent LLM reasoning (Genomics/Pharmacology/Clinical-Evidence/Safety/Chair) | **Partially prototyped.** The notebook's `generate_rationale()` is a single-pass template/API call, not the proposal's adversarial multi-agent debate — the Critic-agent traces in the graphical abstract (EGFR+MET+SRC flagged as redundant vs. FGFR1+JAK2+mTOR passing clean) are illustrative reasoning examples, not executed code. |
| RL feedback layer (L5) | **Not built**, and explicitly scoped as future work in the proposal itself. |

## Real, resolved finding: a patient-identity error was caught and fixed (not a new problem) — now precisely confirmed

`docs/Grant_Proposal.pdf` (v1.1, July 2026) is the most mature, authoritative version of this
proposal — it names the exact patient-identity error for the first time: an initial focal
patient, **`TCGA-A8-A08B`**, was selected via an automated 20-gene mutation scan and presented
as an illustrative TNBC case. Cross-referencing against real PAM50 molecular subtype data later
revealed this patient to be **`BRCA_Her2` subtype, not TNBC**. This was corrected by adding an
explicit subtype-confirmation step (checking for `BRCA_Basal`, the PAM50 proxy for TNBC) before
treating any patient as a TNBC case — yielding the real, confirmed focal patient used
throughout this whole repo, **`TCGA-AO-A128`**, and the same `afatinib + alpelisib +
trastuzumab`, HCOS=0.450 result independently re-verified through real code execution multiple
times elsewhere (`src/regimen/`, `tnbc_regimen_pipeline/`).

The two accidental chat-fragment uploads (`total_88_kinases.txt`, `grounded_result.txt`) that
first surfaced this thread are now both fully resolved: one to a real, cited kinase-panel audit
with a confirmed open discrepancy (`LMTK2`/`LMTK3`, see `docs/limitations.md`), the other to
this precise, already-fixed patient-identity correction.

**One genuinely new finding from the fuller proposal, not previously documented**: the same
patient's TP53 nonsense mutation resolves, under this proposal's NMD-aware zygosity logic, to a
**second, independently-generated hypothesis** — `adavosertib + venetoclax` (WEE1i + BCL-2i) —
ranked lower under the current HCOS scoring not for biological reasons but because the heuristic
rewards shared-pathway overlap, and TP53's pathway doesn't overlap with the EGFR/PTEN/ERBB2
triplet's. The proposal states this as an explicit, self-identified scoring bias and a specific,
pre-registered wet-lab test (Aim 3): if the lower-HCOS regimen matches or beats the higher-HCOS
one in real synergy assays, that's direct evidence the pathway-overlap term is miscalibrated.

## Relationship to the rest of this repo

Three independent things now coexist under one repo, all addressing TNBC kinase-redundancy
and combination therapy from different angles — see each subproject's own docs for how they
diverge:

- **`src/` + `data/`** — CTS/HCOS/DepMap cohort-scale project (90-kinase panel, real
  STRING/DepMap/TCGA-BRCA/DGIdb data, patient-level survival ML)
- **`tnbc-genomics-agent/`** — single-patient VCF pipeline, 10-kinase curated DB, rule-based
  redundancy scoring, single-pass LLM narrative (its companion report is now in
  `tnbc-genomics-agent/docs/`)
- **`mc-ore-prototype/` (this folder)** — proposed extension of `tnbc-genomics-agent`
  specifically, not of the CTS/HCOS project — toward GNN-learned synergy + knowledge graph +
  multi-agent reasoning

None of the three have been cross-validated against each other. Before treating any one as
corroborating another, that needs to be checked explicitly.

## Running the prototype

```bash
cd mc-ore-prototype
pip install pandas
jupyter nbconvert --to notebook --execute tnbc_combo_pipeline.ipynb
```

Or without Jupyter: extract and run the code cells directly (this is exactly how the
byte-for-byte reproduction above was verified in this session, since no Jupyter tooling was
available here either).
