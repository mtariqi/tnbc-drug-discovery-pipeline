"""
Provenance-tracked merge of literature-discovered candidates.

Same merge rule as agentic_regimen_discovery.merge_into_candidate_pool
(ONLY dgidb_confirmed candidates are merged) but additionally returns a
provenance dict: provenance[gene][drug] -> list of PMIDs that supported
it. An empty list means the drug came from the curated dict directly (no
literature source); a non-empty list means it was literature-discovered
and independently DGIdb-confirmed, with the specific PMIDs that produced it.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Tuple

import pandas as pd

from ..utils.validation import validate_gene_drug_schema

logger = logging.getLogger(__name__)


def build_merged_pool_with_provenance(
    discovery_result: pd.DataFrame,
    curated_gene_drugs: Dict[str, List[str]],
) -> Tuple[Dict[str, List[str]], Dict[str, Dict[str, List[str]]]]:
    validate_gene_drug_schema(curated_gene_drugs, "curated_gene_drugs")

    merged = {gene: list(drugs) for gene, drugs in curated_gene_drugs.items()}
    provenance: Dict[str, Dict[str, List[str]]] = {
        gene: {drug: [] for drug in drugs} for gene, drugs in curated_gene_drugs.items()
    }

    if discovery_result is None or discovery_result.empty:
        return merged, provenance

    confirmed = discovery_result[discovery_result["dgidb_confirmed"]]
    for _, row in confirmed.iterrows():
        gene, drug, pmid = row["gene"], row["candidate_drug"], row["pmid"]
        merged.setdefault(gene, [])
        provenance.setdefault(gene, {})
        if drug not in merged[gene]:
            merged[gene].append(drug)
        provenance[gene].setdefault(drug, [])
        if pmid not in provenance[gene][drug]:
            provenance[gene][drug].append(pmid)
        logger.info(f"merged literature candidate gene={gene} drug={drug} pmid={pmid}")

    return merged, provenance
