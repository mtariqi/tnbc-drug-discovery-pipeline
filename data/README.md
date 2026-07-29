# Data Provenance and Regeneration

Every result in this project derives from real, externally-sourced data — no synthetic or
placeholder data contributes to any reported finding (synthetic data is used only for smoke
tests, and never leaves `tests/`). This file exists so any of these datasets can be
re-downloaded from scratch if `data/raw/` is ever wiped or moved to a new machine.

**Real confirmed data root (from actually running `validate_pipeline_inputs.py` and
`run_discovery.py` against it, not guessed):** `~/rtk_nrtk_tnbc/data/`, with `raw/` and
`processed/` subfolders. **`data/raw/kinases/kinase_90_list.txt` is now the real file**
(90 confirmed unique kinase symbols, copied in and verified this session), not a placeholder —
e.g. `~/rtk_nrtk_tnbc/data/raw/kinases/kinase_90_list.txt` on the real machine.
`~/rtk_nrtk_tnbc/data/processed/dgidb/dgidb_interactions.tsv`,
`~/rtk_nrtk_tnbc/data/processed/depmap/depmap_tnbc_essentiality.tsv`,
`~/rtk_nrtk_tnbc/data/processed/tcga_brca/survival_stats.tsv`,
`~/rtk_nrtk_tnbc/data/processed/string/string_edges.tsv` +
`~/rtk_nrtk_tnbc/data/processed/string/community_map.tsv`,
`~/rtk_nrtk_tnbc/data/raw/drugs/drug_list.txt`, `~/rtk_nrtk_tnbc/data/processed/faers/` are
still not present in this repo (large/regenerable — see each source's own section below).

**Unresolved path discrepancy, flagged in the R script itself, not by me — now confirmed on the
real machine:** both `~/rtk_nrtk_tnbc/data/processed/tcga_brca/` (containing `clinical.csv` +
`expression_tpm.csv`, matching `assemble_tcga_brca_data.R`'s actual output convention exactly)
AND a separate `~/rtk_nrtk_tnbc/data/processed/tcga/` (containing a differently-named
`clinical_data.csv`) genuinely exist. `run_track_c_subtyping.py` uses the `tcga_brca/` files
specifically, since their naming matches the R script's real output — but the `tcga/` folder's
origin and purpose hasn't been investigated. Worth checking what wrote `tcga/clinical_data.csv`
and whether it's stale, a duplicate, or serves some other real purpose before assuming either
folder is safe to delete.

This repo's `data/raw/<source>/` layout below is this consolidated repo's own convention, not
identical to the real `~/rtk_nrtk_tnbc/` path — update `scripts/setup_local_data.sh`/`.fish` to
point there specifically rather than a generic placeholder path.

## data/raw/tcga_brca/

| File | Contents | Source |
|---|---|---|
| `clinical.tsv` (or `.csv`) | ~1095-1098 patients: demographics, stage, receptor status (ER/PR/HER2), vital status, follow-up time. **No treatment-response field (pCR/RECIST) — confirmed, not assumed; see `docs/limitations.md`.** | GDC Data Portal, TCGA-BRCA project |
| `expression_star_counts/` | STAR-Counts RNA-seq expression, all patients | GDC via `TCGAbiolinks` (R) |
| `maf/*.maf.gz` | Masked somatic mutation calls, one file per patient/sample (`Hugo_Symbol`, `Tumor_Sample_Barcode`, `Variant_Classification` columns confirmed real against the actual TCGA-AO-A128 file) | GDC Data Portal |

Run `inspect_clinical_columns()` (planned home: `src/data_loaders/tcga_brca_loader.py`) on
your real file before trusting any downstream column-name assumption — TCGA-BRCA clinical
exports vary in whether they carry a `BRCA_Subtype_PAM50` column directly or require deriving
TNBC status from the three receptor-status columns.

## data/raw/depmap/

Broad Institute DepMap Portal, restricted at load-time to the ~25 confirmed-TNBC cell lines
(full pan-cancer files are downloaded; TNBC filtering happens in code, not at download time —
this matters if you later train the pan-cancer model recommended in
`results/reports/Heuristic_vs_ML_Report.docx` §5.1):

- `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`
- `OmicsCNGeneWGS.csv` and `PortalOmicsCNGeneLog2.csv` (both — not yet compared against each other, open item)
- `CRISPRGeneDependency.csv`
- `OmicsSomaticMutationsMatrixDamaging.csv`, `...Hotspot.csv`
- `OmicsFusionFiltered.csv`
- `OmicsInferredMolecularSubtypes.csv` (downloaded, found not applicable to TNBC subtyping — pan-cancer hotspot-defined, not Lehmann-style)
- `CRISPRInferredCommonEssentials.csv`, `AchillesNonessentialControls.csv`
- `Gene.csv`, `CRISPRScreenMap.csv`, `CRISPRConfounders.csv`

## data/raw/string/

Live STRING REST API pull, cached locally: 464 edges among the 90-kinase RTK/NRTK panel, 9
Louvain communities, centrality metrics (betweenness, PageRank, degree).

## data/raw/dgidb/

`dgidb_interactions.tsv` — live DGIdb GraphQL v5 pull, 4655 raw records. Used for both the
CTS druggability term and independent literature-candidate confirmation (must be re-queried
live for confirmation, not read from this cache, since the whole point is independence from
whatever was true at initial-download time).

## data/raw/openfda/

Live openFDA FAERS pull, cached: adverse-event co-occurrence records for the HCOS toxicity
term.

## data/raw/literature_cache/

Cached PubMed ESearch/EFetch results from agentic discovery runs — a cache for repeat runs
against the same query, not a substitute for live confirmation of new candidates.

## data/processed/

Derived, code-generated intermediate files (merged feature matrices, cached CTS scores, etc.)
— nothing here should be hand-edited; regenerate from `data/raw/` via the scripts in
`scripts/` instead. Add `data/processed/*` to `.gitignore` if these get large.
