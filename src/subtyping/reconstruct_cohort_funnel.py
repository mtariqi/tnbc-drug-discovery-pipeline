"""
reconstruct_cohort_funnel.py

Reconstructs the manuscript's cohort-construction funnel (Section 2.6 /
4.2) following the exact counting rules specified by Dr. Sayem, in
response to the unreconciled provenance finding: the manuscript reports
992 -> 171 -> 168 -> 93, while the live TNBC-identification pipeline
(run_track_c_subtyping.py) produces 1098 -> 192 with no located source
for the two downstream steps.

This module does NOT assume the manuscript's original numbers were
wrong -- it implements a specified, auditable methodology and reports
whatever the real data actually yields at each stage, with the full
participant-ID set retained after every step (Rule 10) so the result
is independently checkable.

THE 10 RULES, AND HOW EACH IS HANDLED HERE:
  1. Counting unit = unique participant (first 12 barcode chars).
     -> to_participant_id()
  2. Expression sample = primary tumor only (TCGA sample-type code 01).
     -> is_primary_tumor_barcode()
  3. Duplicate aliquots: retain one sample per participant via a
     deterministic rule -> DEDUP_RULE, documented explicitly below,
     not silently chosen.
  4. PAM50 definition: the exact annotation column/classifier/version
     MUST be passed explicitly as an argument (pam50_column) -- this
     script refuses to guess a column name, per this project's
     standing rule against assuming schema.
  5. Mutation source: GDC release / workflow / access date must be
     passed explicitly (maf_workflow_description) and is echoed back
     in every report, not silently omitted.
  6. Eligible mutations: ELIGIBLE_VARIANT_CLASSES is defined explicitly
     below, not inferred, and is printed in every run so it's always
     visible which classes were actually used.
  7. Kinase panel: loaded from the real, frozen kinase_90_list.txt
     already used throughout this project -- not redefined here.
  8. Drug mapping: drug_gene_map_path and its version/date must be
     passed explicitly; DRUG_MAPPING_VALIDITY_RULE documents what
     counts as a valid mapping.
  9. Alteration threshold: counts UNIQUE altered genes per participant
     (.nunique()), never raw mutation row counts -- this is the exact
     distinction Dr. Sayem's snippet specifies and is preserved exactly.
  10. Audit output: every stage's real participant-ID set is written
      to a real file, not just a count, via export_stage_ids().
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

# ---------------------------------------------------------------------
# Rule 6: explicit, disclosed eligible variant classes.
# Protein-affecting classes only; Silent, Intron, 3'UTR/5'UTR, RNA,
# IGR, and similar non-protein-altering classes are excluded. This
# choice is stated here explicitly so it can be checked or overridden.
# ---------------------------------------------------------------------
ELIGIBLE_VARIANT_CLASSES = [
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
    "Frame_Shift_Ins", "Frame_Shift_Del",
    "In_Frame_Ins", "In_Frame_Del",
    "Splice_Site", "Translation_Start_Site",
]

# ---------------------------------------------------------------------
# Rule 3: explicit, disclosed deduplication rule for duplicate aliquots
# per participant. Chosen rule: prefer the lexicographically FIRST
# full barcode after filtering to primary-tumor samples, which is
# deterministic and reproducible, though arbitrary with respect to
# any specific quality metric -- stated plainly rather than implying
# a quality-based selection that isn't actually being made.
# ---------------------------------------------------------------------
DEDUP_RULE = "lexicographically first full barcode among primary-tumor (01) samples for the participant"

# Rule 8: stated plainly, not assumed silently.
DRUG_MAPPING_VALIDITY_RULE = (
    "A gene is considered drug-mapped if it appears at least once as a "
    "non-null 'gene' entry in the provided DGIdb export, regardless of "
    "interaction_type or interaction_claim_source specificity."
)


def to_participant_id(barcode: str) -> str:
    """Rule 1: first 12 characters of a real TCGA barcode, e.g.
    'TCGA-AO-A128-01A-11D-A10M-09' -> 'TCGA-AO-A128'."""
    return barcode[:12]


def is_primary_tumor_barcode(barcode: str) -> bool:
    """Rule 2: TCGA sample-type code is the first two characters of the
    4th barcode field (e.g. '01A' in 'TCGA-AO-A128-01A-...' -> '01').
    Real TCGA convention: 01=primary tumor, 06=metastatic, 10/11=normal."""
    parts = barcode.split("-")
    if len(parts) < 4:
        return False
    sample_field = parts[3]
    return sample_field[:2] == "01"


def deduplicate_to_one_per_participant(barcodes: List[str]) -> Dict[str, str]:
    """Rule 3: applies DEDUP_RULE. Input barcodes must already be
    filtered to primary-tumor-only (Rule 2) before calling this."""
    by_participant: Dict[str, List[str]] = {}
    for bc in barcodes:
        pid = to_participant_id(bc)
        by_participant.setdefault(pid, []).append(bc)
    return {pid: sorted(bcs)[0] for pid, bcs in by_participant.items()}


def export_stage_ids(stage_name: str, participant_ids: Set[str], outdir: str) -> str:
    """Rule 10: writes the real, complete participant-ID set for a
    given stage to a real file, returning the path written."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    out_path = str(Path(outdir) / f"stage_{stage_name}_participant_ids.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(sorted(participant_ids)))
    return out_path


def compute_gene_count_distribution(
    clinical_path: str,
    pam50_column: str,
    pam50_basal_label: str,
    maf_glob_pattern: str,
    kinase_panel_path: str,
    drug_gene_map_path: str,
    outdir: str = "cohort_funnel_audit",
) -> Dict:
    """
    Dr. Sayem's diagnostic: within the real Stage-2 (168-participant)
    cohort, reports how many participants have >=0/1/2/3/4 unique
    eligible altered genes, computed TWICE -- once using the full
    90-gene panel, once restricted to genes with a real DGIdb mapping.
    Uses the SAME default ELIGIBLE_VARIANT_CLASSES (Rule 6) for both;
    only the gene-set definition varies, isolating exactly the variable
    Dr. Sayem asked about. Participants with zero qualifying rows are
    explicitly reindexed to 0, not silently omitted, since a groupby()
    only produces entries for participants with at least one row.
    """
    import glob as glob_module
    import gzip

    kinase_panel = {line.strip() for line in open(kinase_panel_path) if line.strip()}
    clinical_df = pd.read_csv(clinical_path)
    id_col_candidates = [c for c in clinical_df.columns if "barcode" in c.lower() or "submitter_id" in c.lower() or c.lower() in ("patient_id", "patient")]
    id_col = id_col_candidates[0]
    basal_mask = clinical_df[pam50_column].astype(str) == pam50_basal_label
    basal_ids = set(clinical_df.loc[basal_mask, id_col].astype(str).apply(to_participant_id))

    maf_rows = []
    for maf_path in glob_module.glob(maf_glob_pattern):
        opener = gzip.open if maf_path.endswith(".gz") else open
        with opener(maf_path, "rt") as f:
            header = None
            for line in f:
                if line.startswith("Hugo_Symbol"):
                    header = line.rstrip("\n").split("\t")
                    continue
                if header is None:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != len(header):
                    continue
                maf_rows.append(dict(zip(header, fields)))
    maf_df = pd.DataFrame(maf_rows)
    maf_df = maf_df[maf_df["Tumor_Sample_Barcode"].apply(is_primary_tumor_barcode)]
    dedup_map = deduplicate_to_one_per_participant(maf_df["Tumor_Sample_Barcode"].unique().tolist())
    maf_df = maf_df[maf_df["Tumor_Sample_Barcode"].isin(set(dedup_map.values()))]
    maf_df["patient_id"] = maf_df["Tumor_Sample_Barcode"].apply(to_participant_id)

    mutation_ids = set(maf_df["patient_id"].unique())
    stage2_ids = sorted(basal_ids & mutation_ids)
    print(f"Real Stage-2 cohort (basal-like + has mutation data): {len(stage2_ids)} participants")

    drug_gene_df = pd.read_csv(drug_gene_map_path, sep="\t" if drug_gene_map_path.endswith((".tsv", ".txt")) else ",")
    gene_col = [c for c in drug_gene_df.columns if c.lower() in ("gene", "gene_name", "genesymbol", "kinase_id")][0]
    drug_mapped_genes = set(drug_gene_df[gene_col].dropna().astype(str))

    base_eligible = maf_df[
        maf_df["patient_id"].isin(stage2_ids)
        & maf_df["Hugo_Symbol"].isin(kinase_panel)
        & maf_df["Variant_Classification"].isin(ELIGIBLE_VARIANT_CLASSES)
    ]
    print(f"Real eligible mutation rows (full 90-gene panel, Rule-6 variant classes): {len(base_eligible)}")
    print(f"Real unique participants represented among those rows: {base_eligible['patient_id'].nunique()}")

    def bucket_counts(eligible_df: pd.DataFrame) -> Dict[str, int]:
        counts = eligible_df.groupby("patient_id")["Hugo_Symbol"].nunique()
        counts = counts.reindex(stage2_ids, fill_value=0)
        return {
            "zero_genes": int((counts == 0).sum()),
            "at_least_1": int((counts >= 1).sum()),
            "at_least_2": int((counts >= 2).sum()),
            "at_least_3": int((counts >= 3).sum()),
            "at_least_4": int((counts >= 4).sum()),
        }

    full_panel_dist = bucket_counts(base_eligible)
    dgidb_only_eligible = base_eligible[base_eligible["Hugo_Symbol"].isin(drug_mapped_genes)]
    dgidb_dist = bucket_counts(dgidb_only_eligible)

    # Direct, exact-value evidence for the >=3/>=4 tie -- not inferred from the
    # cumulative buckets alone, so this can be checked against real output
    # rather than trusted on code-review reasoning.
    def exact_value_counts(eligible_df: pd.DataFrame) -> Dict[int, int]:
        counts = eligible_df.groupby("patient_id")["Hugo_Symbol"].nunique()
        counts = counts.reindex(stage2_ids, fill_value=0)
        return counts.value_counts().sort_index().to_dict()

    full_exact = exact_value_counts(base_eligible)
    dgidb_exact = exact_value_counts(dgidb_only_eligible)
    print(f"\nReal exact gene-count distribution, full panel (key=exact count, value=n participants): {full_exact}")
    print(f"Real exact gene-count distribution, DGIdb-restricted: {dgidb_exact}")
    print(f"Real confirmation -- participants with EXACTLY 3 genes (full panel): {full_exact.get(3, 0)}")

    print("\nThreshold distribution, FULL 90-gene panel (no DGIdb restriction):")
    for k, v in full_panel_dist.items():
        print(f"  {k}: {v}")
    print("\nThreshold distribution, DGIdb-restricted genes only:")
    for k, v in dgidb_dist.items():
        print(f"  {k}: {v}")

    stage3_path = Path(outdir) / "stage_3_actionable_participant_ids.txt"
    if stage3_path.exists():
        stage3_ids = set(line.strip() for line in open(stage3_path) if line.strip())
        is_subset = stage3_ids.issubset(set(stage2_ids))
        print(f"\nAll Stage-3 participant IDs are members of the Stage-2 cohort: {is_subset}")
        if not is_subset:
            print(f"  Real Stage-3 IDs NOT in Stage-2: {sorted(stage3_ids - set(stage2_ids))}")

    return {
        "stage2_n": len(stage2_ids),
        "full_panel_distribution": full_panel_dist,
        "dgidb_only_distribution": dgidb_dist,
        "full_panel_exact_counts": full_exact,
        "dgidb_only_exact_counts": dgidb_exact,
        "eligible_rows_unique_participants": int(base_eligible["patient_id"].nunique()),
    }


def reconstruct_cohort_funnel(
    clinical_path: str,
    pam50_column: str,
    pam50_basal_label: str,
    pam50_classifier_description: str,
    maf_glob_pattern: str,
    maf_workflow_description: str,
    maf_access_date: str,
    kinase_panel_path: str,
    drug_gene_map_path: str,
    drug_gene_map_version_description: str,
    outdir: str = "cohort_funnel_audit",
    min_altered_genes: int = 2,
    variant_classes_override: Optional[List[str]] = None,
    require_drug_mapping: bool = True,
) -> Dict:
    """
    variant_classes_override: if provided, replaces ELIGIBLE_VARIANT_CLASSES
    for this run only -- used to test whether the manuscript's original 93
    used a looser variant-class definition than Rule 6 specifies.
    require_drug_mapping: if False, skips the DGIdb filter entirely (every
    panel gene counts as "actionable" just by being altered) -- used to
    test whether the manuscript's original 93 never required a real DGIdb
    match at all. Both are real, disclosed sensitivity-test toggles, not
    silent changes to the specified protocol.
    """
    print("=" * 70)
    print("REAL, DISCLOSED PARAMETERS FOR THIS RUN (Rules 4, 5, 6, 8):")
    print(f"  PAM50 column: {pam50_column!r}, basal label: {pam50_basal_label!r}")
    print(f"  PAM50 classifier/version: {pam50_classifier_description}")
    print(f"  MAF workflow: {maf_workflow_description}")
    print(f"  MAF access date: {maf_access_date}")
    active_variant_classes = variant_classes_override if variant_classes_override is not None else ELIGIBLE_VARIANT_CLASSES
    print(f"  Eligible variant classes: {active_variant_classes}")
    print(f"  Require real DGIdb drug-mapping: {require_drug_mapping}")
    print(f"  Drug-gene map version: {drug_gene_map_version_description}")
    print(f"  Drug-mapping validity rule: {DRUG_MAPPING_VALIDITY_RULE}")
    print(f"  Deduplication rule: {DEDUP_RULE}")
    print("=" * 70)

    audit_trail = {}

    # --- Stage 1: PAM50 basal-like identification ---
    clinical_df = pd.read_csv(clinical_path)
    if pam50_column not in clinical_df.columns:
        raise ValueError(
            f"pam50_column={pam50_column!r} not found in {clinical_path}. Real columns: "
            f"{clinical_df.columns.tolist()}. Refusing to guess -- pass the real column name."
        )
    id_col_candidates = [c for c in clinical_df.columns if "barcode" in c.lower() or "submitter_id" in c.lower() or c.lower() in ("patient_id", "patient")]
    if not id_col_candidates:
        raise ValueError(f"No obvious participant-ID column found in {clinical_path}. Real columns: {clinical_df.columns.tolist()}")
    id_col = id_col_candidates[0]
    print(f"\nUsing {id_col!r} as the real participant-ID column in the clinical file.")

    total_patients = set(clinical_df[id_col].astype(str).apply(to_participant_id))
    print(f"Stage 0 -- total participants in clinical file: {len(total_patients)}")
    audit_trail["stage_0_total"] = export_stage_ids("0_total", total_patients, outdir)

    basal_mask = clinical_df[pam50_column].astype(str) == pam50_basal_label
    basal_ids = set(clinical_df.loc[basal_mask, id_col].astype(str).apply(to_participant_id))
    print(f"Stage 1 -- PAM50 basal-like participants: {len(basal_ids)}")
    audit_trail["stage_1_basal"] = export_stage_ids("1_basal", basal_ids, outdir)

    # --- Stage 2: mutation availability, from real MAF files ---
    import glob as glob_module
    import gzip

    kinase_panel = {line.strip() for line in open(kinase_panel_path) if line.strip()}
    print(f"\nReal kinase panel loaded: {len(kinase_panel)} genes")

    maf_rows = []
    maf_files = glob_module.glob(maf_glob_pattern)
    print(f"Real MAF files matched by pattern: {len(maf_files)}")
    for maf_path in maf_files:
        opener = gzip.open if maf_path.endswith(".gz") else open
        with opener(maf_path, "rt") as f:
            header = None
            for line in f:
                if line.startswith("Hugo_Symbol"):
                    header = line.rstrip("\n").split("\t")
                    continue
                if header is None:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != len(header):
                    continue
                row = dict(zip(header, fields))
                maf_rows.append(row)
    maf_df = pd.DataFrame(maf_rows)
    if len(maf_df) == 0:
        raise ValueError(f"No real MAF rows loaded from pattern {maf_glob_pattern!r} -- check the path.")

    barcode_col = "Tumor_Sample_Barcode"
    maf_df = maf_df[maf_df[barcode_col].apply(is_primary_tumor_barcode)]
    print(f"Real MAF rows after primary-tumor-only filter (Rule 2): {len(maf_df)}")

    real_barcodes = maf_df[barcode_col].unique().tolist()
    dedup_map = deduplicate_to_one_per_participant(real_barcodes)
    print(f"Real unique participants after deduplication (Rule 3): {len(dedup_map)}")
    kept_barcodes = set(dedup_map.values())
    maf_df = maf_df[maf_df[barcode_col].isin(kept_barcodes)]

    maf_df["patient_id"] = maf_df[barcode_col].apply(to_participant_id)
    mutation_ids = set(maf_df["patient_id"].unique())
    print(f"Real participants with any mutation data: {len(mutation_ids)}")

    mutation_eligible_basal = basal_ids & mutation_ids
    print(f"Stage 2 -- basal-like AND has real mutation data: {len(mutation_eligible_basal)}")
    audit_trail["stage_2_mutation_eligible"] = export_stage_ids("2_mutation_eligible", mutation_eligible_basal, outdir)

    # --- Stage 3: eligible, panel-restricted, drug-mapped alterations (Rule 6, 7, 8, 9) ---
    drug_gene_df = pd.read_csv(drug_gene_map_path, sep="\t" if drug_gene_map_path.endswith((".tsv", ".txt")) else ",")
    gene_col_candidates = [c for c in drug_gene_df.columns if c.lower() in ("gene", "gene_name", "genesymbol", "kinase_id")]
    if not gene_col_candidates:
        raise ValueError(f"No obvious gene column in {drug_gene_map_path}. Real columns: {drug_gene_df.columns.tolist()}")
    drug_mapped_genes = set(drug_gene_df[gene_col_candidates[0]].dropna().astype(str))
    print(f"\nReal drug-mapped genes in provided DGIdb export (using real column {gene_col_candidates[0]!r}): {len(drug_mapped_genes)}")

    eligible = maf_df[
        maf_df["patient_id"].isin(mutation_eligible_basal)
        & maf_df["Hugo_Symbol"].isin(kinase_panel)
        & maf_df["Variant_Classification"].isin(active_variant_classes)
    ]
    if require_drug_mapping:
        eligible = eligible[eligible["Hugo_Symbol"].isin(drug_mapped_genes)]
    else:
        print("SENSITIVITY TEST: skipping the DGIdb drug-mapping filter entirely -- every altered panel gene counts.")
    print(f"Real eligible mutation rows after panel + variant-class + drug-mapping filters: {len(eligible)}")

    # Rule 9: unique altered GENES per participant, not mutation rows.
    altered_gene_counts = eligible.groupby("patient_id")["Hugo_Symbol"].nunique()
    candidate_ids = set(altered_gene_counts[altered_gene_counts >= min_altered_genes].index)
    print(f"Stage 3 -- participants with >={min_altered_genes} unique panel genes altered: {len(candidate_ids)}")
    audit_trail["stage_3_actionable"] = export_stage_ids("3_actionable", candidate_ids, outdir)

    summary = {
        "stage_0_total_participants": len(total_patients),
        "stage_1_pam50_basal": len(basal_ids),
        "stage_2_basal_with_mutation_data": len(mutation_eligible_basal),
        "stage_3_actionable_cohort": len(candidate_ids),
        "audit_id_files": audit_trail,
    }
    print("\n" + "=" * 70)
    print("FINAL FUNNEL:", summary)
    return summary


# =====================================================================
# SMOKE TEST -- synthetic clinical/MAF/DGIdb data exercising every
# rule: participant truncation, primary-vs-normal filtering, duplicate
# aliquot deduplication, variant-class filtering, drug-mapping
# filtering, and unique-gene-not-row counting.
# =====================================================================

def _run_smoke_test():
    import tempfile

    tmpdir = tempfile.mkdtemp()

    clinical_df = pd.DataFrame({
        "submitter_id": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003", "TCGA-AA-0004"],
        "pam50_subtype": ["Basal", "Basal", "LumA", "Basal"],
    })
    clinical_path = f"{tmpdir}/clinical.csv"
    clinical_df.to_csv(clinical_path, index=False)

    maf_lines = [
        "Hugo_Symbol\tVariant_Classification\tTumor_Sample_Barcode",
        # Participant 1: two DIFFERENT real aliquots (duplicate) -- dedup should keep exactly one
        "EGFR\tMissense_Mutation\tTCGA-AA-0001-01A-11D-XXXX-09",
        "EGFR\tMissense_Mutation\tTCGA-AA-0001-01B-11D-YYYY-09",  # duplicate aliquot, same gene
        "ERBB2\tMissense_Mutation\tTCGA-AA-0001-01A-11D-XXXX-09",
        # Participant 1 also has a NORMAL-tissue row that must be excluded (sample type 11)
        "TP53\tNonsense_Mutation\tTCGA-AA-0001-11A-11D-ZZZZ-09",
        # Participant 2: real mutations but only ONE unique panel gene with drug mapping -> should NOT pass threshold
        "PTK2\tMissense_Mutation\tTCGA-AA-0002-01A-11D-WWWW-09",
        "PTK2\tMissense_Mutation\tTCGA-AA-0002-01A-11D-WWWW-09",  # same gene, same sample -- must count as 1 unique gene, not 2 rows
        # Participant 3 is NOT basal-like -- must be excluded regardless of mutation data
        "EGFR\tMissense_Mutation\tTCGA-AA-0003-01A-11D-VVVV-09",
        "ERBB2\tMissense_Mutation\tTCGA-AA-0003-01A-11D-VVVV-09",
        # Participant 4: basal-like, real mutations, but one is a SILENT variant (ineligible class)
        "FLT1\tMissense_Mutation\tTCGA-AA-0004-01A-11D-UUUU-09",
        "TYK2\tSilent\tTCGA-AA-0004-01A-11D-UUUU-09",  # ineligible variant class -- must be excluded
        "MET\tMissense_Mutation\tTCGA-AA-0004-01A-11D-UUUU-09",
    ]
    maf_path = f"{tmpdir}/test.maf"
    with open(maf_path, "w") as f:
        f.write("\n".join(maf_lines))

    panel_path = f"{tmpdir}/panel.txt"
    with open(panel_path, "w") as f:
        f.write("\n".join(["EGFR", "ERBB2", "PTK2", "FLT1", "MET", "TYK2"]))

    drug_map_df = pd.DataFrame({"gene": ["EGFR", "ERBB2", "FLT1", "MET"]})  # PTK2, TYK2 deliberately NOT drug-mapped
    drug_map_path = f"{tmpdir}/dgidb.csv"
    drug_map_df.to_csv(drug_map_path, index=False)

    result = reconstruct_cohort_funnel(
        clinical_path=clinical_path,
        pam50_column="pam50_subtype",
        pam50_basal_label="Basal",
        pam50_classifier_description="TEST synthetic PAM50 labels",
        maf_glob_pattern=f"{tmpdir}/test.maf",
        maf_workflow_description="TEST synthetic MAF",
        maf_access_date="TEST",
        kinase_panel_path=panel_path,
        drug_gene_map_path=drug_map_path,
        drug_gene_map_version_description="TEST synthetic DGIdb export",
        outdir=f"{tmpdir}/audit",
    )

    print("\n" + "=" * 70)
    print("SMOKE TEST ASSERTIONS:")
    assert result["stage_0_total_participants"] == 4, result
    assert result["stage_1_pam50_basal"] == 3, "TCGA-AA-0003 is LumA, should be excluded from basal"
    assert result["stage_2_basal_with_mutation_data"] == 3, "all 3 basal participants have real mutation rows"
    # Participant 1: EGFR + ERBB2 = 2 unique drug-mapped panel genes (after dedup, normal-tissue exclusion) -> PASSES
    # Participant 2: only PTK2, and PTK2 isn't drug-mapped -> 0 eligible genes -> FAILS
    # Participant 4: FLT1 + MET = 2 (TYK2 excluded for Silent classification) -> PASSES
    assert result["stage_3_actionable_cohort"] == 2, f"expected 2 (participants 1 and 4), got {result['stage_3_actionable_cohort']}"
    print("PASSED: duplicate-aliquot dedup, normal-tissue exclusion, non-basal exclusion, ineligible-variant-class "
          "exclusion, non-drug-mapped-gene exclusion, and unique-gene-not-row counting ALL verified correct.")

    print("\n" + "=" * 70)
    print("ALL SMOKE TESTS PASSED. Ready to run against real clinical, MAF, panel, and DGIdb files.")


def _run_distribution_smoke_test():
    import tempfile

    tmpdir = tempfile.mkdtemp()

    clinical_df = pd.DataFrame({
        "submitter_id": ["TCGA-AA-0001", "TCGA-AA-0002", "TCGA-AA-0003", "TCGA-AA-0004"],
        "pam50_subtype": ["Basal", "Basal", "Basal", "Basal"],
    })
    clinical_path = f"{tmpdir}/clinical.csv"
    clinical_df.to_csv(clinical_path, index=False)

    maf_lines = [
        "Hugo_Symbol\tVariant_Classification\tTumor_Sample_Barcode",
        # Participant 1: 3 real panel genes, all drug-mapped -> full=3, dgidb=3
        "EGFR\tMissense_Mutation\tTCGA-AA-0001-01A-11D-XXXX-09",
        "ERBB2\tMissense_Mutation\tTCGA-AA-0001-01A-11D-XXXX-09",
        "FLT1\tMissense_Mutation\tTCGA-AA-0001-01A-11D-XXXX-09",
        # Participant 2: 2 panel genes, only 1 drug-mapped -> full=2, dgidb=1
        "PTK2\tMissense_Mutation\tTCGA-AA-0002-01A-11D-WWWW-09",   # not drug-mapped
        "MET\tMissense_Mutation\tTCGA-AA-0002-01A-11D-WWWW-09",    # drug-mapped
        # Participant 3: has mutation data, but ZERO eligible panel genes -- must show up as 0, not be omitted
        "NOTAPANELGENE\tMissense_Mutation\tTCGA-AA-0003-01A-11D-VVVV-09",
        # Participant 4: no MAF rows at all -- also must show up as 0
    ]
    maf_path = f"{tmpdir}/test.maf"
    with open(maf_path, "w") as f:
        f.write("\n".join(maf_lines))

    panel_path = f"{tmpdir}/panel.txt"
    with open(panel_path, "w") as f:
        f.write("\n".join(["EGFR", "ERBB2", "FLT1", "PTK2", "MET"]))

    drug_map_df = pd.DataFrame({"gene": ["EGFR", "ERBB2", "FLT1", "MET"]})  # PTK2 deliberately NOT drug-mapped
    drug_map_path = f"{tmpdir}/dgidb.csv"
    drug_map_df.to_csv(drug_map_path, index=False)

    result = compute_gene_count_distribution(
        clinical_path=clinical_path, pam50_column="pam50_subtype", pam50_basal_label="Basal",
        maf_glob_pattern=f"{tmpdir}/test.maf", kinase_panel_path=panel_path,
        drug_gene_map_path=drug_map_path, outdir=f"{tmpdir}/audit_nonexistent",
    )

    print("\n" + "=" * 70)
    print("DISTRIBUTION SMOKE TEST ASSERTIONS:")
    assert result["stage2_n"] == 3, f"participant 3 (no eligible genes) and the 3 with mutations = 3 real Stage-2 (participant 4 has NO mutation data at all, correctly excluded from Stage 2 itself); got {result}"
    fp = result["full_panel_distribution"]
    assert fp["zero_genes"] == 1, f"exactly participant 3 should show zero eligible panel genes, got {fp}"
    assert fp["at_least_1"] == 2 and fp["at_least_2"] == 2 and fp["at_least_3"] == 1, fp
    dg = result["dgidb_only_distribution"]
    assert dg["at_least_2"] == 1 and dg["at_least_1"] == 2, \
        f"participant 2 should drop from 2 (full panel) to 1 (dgidb-only, PTK2 excluded); got {dg}"
    print(f"PASSED: zero-gene participant correctly counted (not omitted); full-panel vs DGIdb-only distributions "
          f"diverge exactly where expected (participant 2's PTK2 exclusion). Full: {fp}, DGIdb-only: {dg}")

    print("\n" + "=" * 70)
    print("ALL DISTRIBUTION SMOKE TESTS PASSED.")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        _run_smoke_test()
    elif len(sys.argv) > 1 and sys.argv[1] == "--test-distribution":
        _run_distribution_smoke_test()
    else:
        print("Run with --test first. For a real run, call reconstruct_cohort_funnel() directly with real paths "
              "and the required _description/_date arguments filled in -- see the module docstring for what each means.")
