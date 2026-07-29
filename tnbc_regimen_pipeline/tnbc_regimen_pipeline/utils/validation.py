"""
Manual schema validation for the Dict[str, List[str]] gene->drug shape
used throughout this package (curated_gene_drugs, real_gene_drugs,
merged_gene_drugs).

Hand-written rather than pydantic: pydantic was not available in the
sandbox this package was originally built in (no network to install it),
and the structure being validated is a shallow dict of lists of strings
-- proportionate to check by hand. If the schema grows more nested,
revisit this decision.
"""

from __future__ import annotations

from typing import Dict, List


def validate_gene_drug_schema(d: Dict[str, List[str]], name: str = "gene_drug_dict") -> None:
    """Raises ValueError with a specific, actionable message on the first
    structural problem found. Deliberately does not attempt partial
    recovery -- a malformed input should fail before reaching
    merging/scoring, not produce a silently wrong regimen downstream."""
    if not isinstance(d, dict):
        raise ValueError(f"{name} must be a dict, got {type(d).__name__}")
    for gene, drugs in d.items():
        if not isinstance(gene, str) or not gene.strip():
            raise ValueError(f"{name}: gene keys must be non-empty strings, got {gene!r}")
        if not isinstance(drugs, list):
            raise ValueError(f"{name}[{gene!r}] must be a list of drug names, got {type(drugs).__name__}")
        for drug in drugs:
            if not isinstance(drug, str) or not drug.strip():
                raise ValueError(f"{name}[{gene!r}] contains a non-string or empty drug entry: {drug!r}")
