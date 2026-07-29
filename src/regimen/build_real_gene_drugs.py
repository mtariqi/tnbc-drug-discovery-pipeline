"""
Real DGIdb-Sourced GENE_DRUGS for MDCOE
==========================================

mdcoe.py's GENE_DRUGS is a small, hand-curated dictionary (~20 genes,
1-3 drugs each). This builds a REAL, DGIdb-sourced alternative from your
actual dgidb_interactions.tsv -- but does NOT simply dump raw DGIdb
output in, for two reasons already proven necessary this session:

    1. Raw DGIdb output is heavily contaminated with research/lab compound
       codes and even occasional gene-symbol-as-drug-name noise (see the
       TNBC-drug-regimen-discovery repo's build_drug_list*.py iterations --
       it took 3 rounds of filtering to get a clean list). This reuses
       that exact, validated filter rather than re-inventing it.

    2. Uncapped, DGIdb often returns dozens of drugs per gene -- this
       would make MDCOE's beam search space explode. Capped to a
       reasonable number per gene, similar to the drug_list.txt curation
       approach.

IMPORTANT: this produces an ADDITIONAL, separate dictionary
(REAL_GENE_DRUGS), not a silent replacement of mdcoe.py's GENE_DRUGS --
your already-validated, real MDCOE results (the afatinib+alpelisib+
trastuzumab finding) depend on the hand-curated dictionary and must
remain reproducible exactly as before. Use REAL_GENE_DRUGS as an
explicit, opt-in alternative or supplement.
"""

from __future__ import annotations

import re
from typing import Dict, List

import pandas as pd

# Reused exactly from TNBC-drug-regimen-discovery's build_drug_list4.py,
# validated against real DGIdb output across 3 rounds of testing this session.
GENERIC_STOPLIST = {
    "human", "agent", "factor", "growth", "cytokine", "inhibitor", "recombinant",
    "antineoplastic", "autologous", "therapeutic", "isomer", "mab", "pegylated",
    "protein-bound", "tosylate", "esylate", "s-malate", "car-specific", "t-lymphocytes",
    "ii", "glulisine", "insulin", "regular", "arginase", "2", "i", "anhydrous", "fibroblast",
    "hydrochloride", "monohydrate", "sodium", "citrate", "tartrate", "succinate",
}

INN_SUFFIXES = (
    "nib", "rafenib", "tinib", "parib", "ciclib", "mab", "zumab", "ximab", "tuzumab",
    "umab", "stat", "degib", "lisib", "setron", "sertib", "clax", "dostat", "metinib",
    "afenib", "bulin", "gene", "kinra", "cept", "mustine", "platin", "rubicin",
)


def _clean(name: str) -> str:
    return name.lower().strip().rstrip(",").strip()


def _looks_like_real_drug(raw_name: str, exclude_gene_names: set = frozenset()) -> bool:
    n = _clean(raw_name)
    if not n or n in GENERIC_STOPLIST:
        return False
    if n in exclude_gene_names:  # e.g. 'syk', 'ror1' showing up as a "drug" -- noise in DGIdb's own source aggregation
        return False
    if n.endswith(INN_SUFFIXES):
        return True
    if re.search(r"\[pmid", n) or re.search(r"\bet al\.", n) or re.search(r"^compound\s", n):
        return False
    if re.search(r"chembl[:\s]?\d", n):
        return False
    if re.match(r"^\d+-[a-z]+$", n):  # legitimate locant prefix: 5-fluorouracil
        return True
    stripped = re.sub(r"[-\s]", "", n)
    if re.fullmatch(r"[a-z0-9]+", stripped) and re.search(r"[a-z]", stripped) and re.search(r"\d", stripped):
        return False  # lab code: nvp-bhg712, tg100-801, xmd8-92
    if re.search(r"-in-\d+$", n) or re.search(r"inhibitor \d+[a-z]?$", n):
        return False
    return True


def build_real_gene_drugs(
    dgidb_interactions_path: str,
    genes: List[str],
    max_drugs_per_gene: int = 3,
    known_validated: Dict[str, List[str]] = None,
    all_kinase_panel: List[str] = None,
) -> Dict[str, List[str]]:
    """
    Builds a REAL_GENE_DRUGS dict: gene -> [drug, drug, ...], filtered to
    real drugs only (via the validated stoplist/INN-suffix filter) and
    capped to `max_drugs_per_gene` by source count, same methodology as
    the drug_list.txt curation.

    all_kinase_panel: your FULL kinase panel (e.g. all 90), used as the
    gene-name-exclusion reference -- IMPORTANT: this must be the full
    panel, not just `genes`. DGIdb's noise (a kinase gene symbol
    mislabeled as a "drug") can involve ANY kinase in your panel, not
    only the ones you happen to be building entries for in this specific
    call. Defaults to `genes` if not given, but pass the full panel
    explicitly whenever you're only scoring a subset.

    known_validated: optional {gene: [drug, ...]} to union in explicitly,
    so already-confirmed-real drugs (e.g. from mdcoe.py's own curated
    GENE_DRUGS) can never be silently dropped by the automated filter,
    same discipline used for drug_list.txt's known_validated set.
    """
    exclude_gene_names = set(g.lower() for g in (all_kinase_panel or genes))

    dgidb = pd.read_csv(dgidb_interactions_path, sep="\t")
    dgidb["drug"] = dgidb["drug"].astype(str).str.lower()
    dgidb = dgidb[dgidb["kinase_id"].isin(genes)]

    dgidb_clean = dgidb[dgidb["drug"].apply(lambda d: _looks_like_real_drug(d, exclude_gene_names=exclude_gene_names))].copy()
    n_dropped = dgidb["drug"].nunique() - dgidb_clean["drug"].nunique()
    print(f"Dropped {n_dropped} junk/code-named entries "
          f"({dgidb['drug'].nunique()} -> {dgidb_clean['drug'].nunique()} unique drugs) "
          f"across {dgidb_clean['kinase_id'].nunique()} genes")

    dgidb_clean["n_sources"] = dgidb_clean["sources"].fillna("").apply(
        lambda s: len([x for x in str(s).split(",") if x.strip()])
    )
    top_per_gene = (
        dgidb_clean.sort_values("n_sources", ascending=False)
        .groupby("kinase_id")
        .head(max_drugs_per_gene)
    )

    real_gene_drugs: Dict[str, List[str]] = {}
    for gene, group in top_per_gene.groupby("kinase_id"):
        real_gene_drugs[gene] = sorted(group["drug"].unique())

    if known_validated:
        for gene, drugs in known_validated.items():
            existing = set(real_gene_drugs.get(gene, []))
            real_gene_drugs[gene] = sorted(existing | set(d.lower() for d in drugs))

    genes_with_no_real_drug = [g for g in genes if g not in real_gene_drugs]
    if genes_with_no_real_drug:
        print(f"Note: {len(genes_with_no_real_drug)} genes have no confirmed real drug after filtering "
              f"(consistent with mdcoe.py's own design -- genes with no real inhibitor are deliberately "
              f"left empty rather than forcing a fabricated mapping): {genes_with_no_real_drug[:10]}"
              f"{'...' if len(genes_with_no_real_drug) > 10 else ''}")

    return real_gene_drugs


# =====================================================================
# SMOKE TEST -- against the exact contamination patterns confirmed real
# =====================================================================

def _run_smoke_test():
    rows = [
        ("EGFR", "afatinib", "DrugBank,ChEMBL,FDA"),
        ("EGFR", "erlotinib", "DrugBank,ChEMBL"),
        ("EGFR", "gefitinib", "DrugBank"),
        ("EGFR", "cp-459632", "ChEMBL"),                    # lab code -- must be dropped
        ("EGFR", "compound 19a [pmid: 21855335]", "ChEMBL"),  # literature citation -- must be dropped
        ("ROR1", "syk", "SomeNoisySource"),                  # gene name masquerading as drug -- must be dropped
        ("OBSCURE_KINASE", "azd1480", "ChEMBL"),             # lab code, only candidate for this gene -- dropped, gene ends up empty
    ]
    df = pd.DataFrame(rows, columns=["kinase_id", "drug", "sources"])
    df.to_csv("/tmp/test_dgidb_for_mdcoe.tsv", sep="\t", index=False)

    genes = ["EGFR", "ROR1", "OBSCURE_KINASE"]
    full_panel = genes + ["SYK", "JAK1", "ABL1"]  # realistic: the full 90-kinase panel is always the exclusion reference, even when scoring a subset
    result = build_real_gene_drugs(
        "/tmp/test_dgidb_for_mdcoe.tsv", genes, max_drugs_per_gene=3,
        known_validated={"EGFR": ["afatinib"]},  # already-known-real, must survive regardless
        all_kinase_panel=full_panel,
    )
    print(result)
    print()

    assert result["EGFR"] == ["afatinib", "erlotinib", "gefitinib"], f"got {result['EGFR']}"
    assert "cp-459632" not in result["EGFR"]
    assert "compound 19a [pmid: 21855335]" not in result["EGFR"]
    assert "ROR1" not in result or "syk" not in result.get("ROR1", [])
    assert "OBSCURE_KINASE" not in result, "gene with only a junk-filtered candidate should end up with no entry, not a fabricated one"

    print("PASSED: real drugs correctly kept and capped, lab codes and literature")
    print("citations correctly dropped, gene-name-as-drug noise correctly dropped,")
    print("and a gene left with zero real candidates correctly has no entry")
    print("(matching mdcoe.py's own design philosophy) rather than a fabricated one.")


if __name__ == "__main__":
    _run_smoke_test()
