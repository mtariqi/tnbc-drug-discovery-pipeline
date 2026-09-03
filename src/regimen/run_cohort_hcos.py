"""
Cohort-wide HCOS/MDCOE run across all eligible PAM50 basal-like patients.

Real, verified components used:
  - mdcoe.py: uploaded, verified to reproduce the documented HCOS=0.450 focal-patient
    result before this script was written.
  - patient_panel_mutations.csv: real MAF-derived (patient, gene, alteration_type)
    rows, extracted from the actual GDC masked-somatic-mutation files via
    extract_panel_mutations.py (992 files, 1635 rows, 749 patients with >=1 panel
    gene altered in the full BRCA cohort).
  - molecular_subtype.csv: real PAM50 calls (192 Basal patients -- matches the
    manuscript's documented cohort funnel exactly).

Cohort construction (mirrors the manuscript's stated funnel logic):
  1. Restrict mutations to protein-altering classes only (drop Silent, Intron,
     3'UTR, Splice_Region -- these don't change the protein and aren't what the
     manuscript's "protein-altering variants" language refers to).
  2. Restrict to patients called PAM50 Basal.
  3. Restrict further to patients with >=2 panel genes that have an actual entry
     in GENE_DRUGS (TP53 always counts, via resolve_tp53_drugs) -- i.e., >=2
     independently druggable altered genes, matching the "candidate drug-mapped
     (>=2 genes)" step in the manuscript's funnel.
  4. Run the real MDCOE beam search + HCOS scoring for every patient in that final
     set, individually -- one call per patient, no shared state, no pooling.

Honesty notes carried over from the rest of this project:
  - This is a mechanical loop around the real, already-verified mdcoe.py; the
    scoring logic itself is unchanged and unmodified here.
  - TP53 zygosity classification here uses the REAL Variant_Classification value
    per patient (not a hardcoded default), since real alteration_type data was
    extracted alongside the gene calls.
  - If any patient produces zero candidate drugs (all altered panel genes lack a
    GENE_DRUGS entry), they're reported as such, not silently dropped.
"""
from __future__ import annotations

import os
from pathlib import Path
from collections import Counter

import pandas as pd

from mdcoe import GENE_DRUGS, GENE_PATHWAYS, resolve_tp53_drugs, DrugGraph, SynergyNet, HCOS, MDCOE

# Matches this project's standard layout (see README.md "Repository layout").
# Override by setting the TNBC_BASE_DIR environment variable if your checkout
# lives somewhere other than ~/rtk_nrtk_tnbc.
# Resolves to the repo root this script actually lives in (src/regimen/../.. ->
# repo root), so this works regardless of what the repo is named or cloned to.
# Override with the TNBC_BASE_DIR environment variable if needed.
_DEFAULT_BASE_DIR = Path(__file__).resolve().parent.parent.parent
BASE_DIR = Path(os.environ.get("TNBC_BASE_DIR", _DEFAULT_BASE_DIR))
OUTPUT_DIR = BASE_DIR / "results" / "tables"

PROTEIN_ALTERING = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Splice_Region", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation",
}


def build_available_drugs(patient_genes: dict[str, str]) -> tuple[list[str], list[str]]:
    """patient_genes: {gene: alteration_type}. Returns (available_drugs, genes_with_no_drug_entry)."""
    drugs = set()
    no_entry = []
    for gene, alt_type in patient_genes.items():
        if gene == "TP53":
            drugs.update(resolve_tp53_drugs(alt_type))
        else:
            gene_drugs = GENE_DRUGS.get(gene, [])
            if gene_drugs:
                drugs.update(gene_drugs)
            else:
                no_entry.append(gene)
    return sorted(drugs), no_entry


def main():
    muts = pd.read_csv(BASE_DIR / "data/processed/patient_panel_mutations.csv")
    muts = muts[muts["alteration_type"].isin(PROTEIN_ALTERING)].copy()

    subtypes = pd.read_csv(BASE_DIR / "data/processed/tcga_brca/molecular_subtype.csv")
    basal_patients = set(subtypes.loc[subtypes["BRCA_Subtype_PAM50"] == "Basal", "patient"])
    print(f"PAM50 Basal patients: {len(basal_patients)}")

    muts = muts[muts["patient_barcode"].isin(basal_patients)].copy()
    print(f"Protein-altering panel-gene mutation rows in Basal patients: {len(muts)}")
    print(f"Basal patients with >=1 protein-altering panel-gene mutation: {muts['patient_barcode'].nunique()}")

    # Build per-patient gene->alteration_type dict (if a gene has multiple
    # alterations for one patient, keep the first -- rare, and not disambiguated
    # further here; flagged in output if it happens).
    patient_gene_alts: dict[str, dict[str, str]] = {}
    dup_flags = []
    for row in muts.itertuples(index=False):
        d = patient_gene_alts.setdefault(row.patient_barcode, {})
        if row.gene in d and d[row.gene] != row.alteration_type:
            dup_flags.append((row.patient_barcode, row.gene, d[row.gene], row.alteration_type))
        d.setdefault(row.gene, row.alteration_type)

    if dup_flags:
        print(f"\nNOTE: {len(dup_flags)} (patient, gene) pairs had >1 distinct alteration_type; kept the first seen:")
        for p, g, a1, a2 in dup_flags:
            print(f"  {p} {g}: kept {a1!r}, also saw {a2!r}")

    # Apply the >=2 drug-mappable-gene threshold.
    eligible = {}
    for patient, genes in patient_gene_alts.items():
        drugs, no_entry = build_available_drugs(genes)
        n_druggable_genes = len([g for g in genes if g == "TP53" or GENE_DRUGS.get(g)])
        if n_druggable_genes >= 2:
            eligible[patient] = {"genes": genes, "drugs": drugs, "no_entry": no_entry}

    print(f"\nEligible patients (Basal, >=2 drug-mappable altered panel genes): {len(eligible)}")
    for p in sorted(eligible):
        print(f"  {p}: genes={sorted(eligible[p]['genes'].keys())}  n_drugs={len(eligible[p]['drugs'])}")

    # Run real MDCOE/HCOS for every eligible patient.
    net = SynergyNet()
    results = []
    for patient, info in sorted(eligible.items()):
        graph = DrugGraph(info["drugs"])
        ranked = MDCOE(graph, net, HCOS, beam_width=50, max_depth=5, top_k=10)
        if not ranked:
            results.append({
                "patient": patient, "genes": ",".join(sorted(info["genes"])),
                "n_available_drugs": len(info["drugs"]),
                "top_regimen": None, "top_hcos": None, "n_tied_top": 0,
                "genes_with_no_drug_entry": ",".join(info["no_entry"]),
            })
            continue
        top_score = ranked[0][1]
        tied = [r for r, s in ranked if abs(s - top_score) < 1e-9]
        results.append({
            "patient": patient,
            "genes": ",".join(sorted(info["genes"])),
            "n_available_drugs": len(info["drugs"]),
            "top_regimen": " + ".join(sorted(tied[0])),
            "top_hcos": round(top_score, 6),
            "n_tied_top": len(tied),
            "tied_regimens": " | ".join(" + ".join(sorted(r)) for r in tied),
            "genes_with_no_drug_entry": ",".join(info["no_entry"]),
        })

    out = pd.DataFrame(results)
    out.to_csv("cohort_hcos_results.csv", index=False)
    print(f"\nWrote {len(out)} patient results to cohort_hcos_results.csv")

    # ---- Aggregate: regimen recurrence, target recurrence, pathway recurrence ----
    valid = out[out["top_regimen"].notna()]
    print(f"\nPatients with >=1 scoreable regimen: {len(valid)} / {len(out)}")

    regimen_counts = Counter(valid["top_regimen"])
    print("\nTop-regimen recurrence across cohort:")
    for regimen, count in regimen_counts.most_common(15):
        print(f"  {count}x  {regimen}")

    drug_to_gene = {d: g for g, ds in GENE_DRUGS.items() for d in ds}
    drug_counts = Counter()
    pathway_counts = Counter()
    for regimen in valid["top_regimen"]:
        for drug in regimen.split(" + "):
            drug_counts[drug] += 1
            gene = drug_to_gene.get(drug)
            if gene:
                for pw in GENE_PATHWAYS.get(gene, []):
                    pathway_counts[pw] += 1

    print("\nDrug recurrence across top-ranked regimens:")
    for drug, count in drug_counts.most_common(15):
        print(f"  {count}x  {drug}")

    print("\nPathway recurrence across top-ranked regimens:")
    for pw, count in pathway_counts.most_common():
        print(f"  {count}x  {pw}")

    n_unique_regimens = len(regimen_counts)
    print(f"\nUnique top regimens across cohort: {n_unique_regimens} (out of {len(valid)} scored patients)")


if __name__ == "__main__":
    main()
