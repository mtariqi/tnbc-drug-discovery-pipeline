"""
Before/after regimen diffing.

Runs cohort-wide scoring twice -- once with an original gene->drug dict,
once with a modified (e.g. discovery-merged) one -- and reports only the
patients whose actual top-ranked regimen or HCOS score changed, plus the
magnitude of that change and a specific reason. This is the honest way to
report a change's real impact: not "N candidates were confirmed" (a fact
about a discovery table) but "M patients' actual optimum changed, by how
much, because of what" (a fact about the thing that matters).
"""

from __future__ import annotations

import logging
from typing import Callable, Dict, List, Optional

import pandas as pd

from ..cohort.cohort_wide_regimen_analysis import run_cohort_wide_analysis

logger = logging.getLogger(__name__)


def _regimen_to_set(regimen: Optional[str]) -> set:
    return set(regimen.split(" + ")) if regimen else set()


def _reason_for_change(row) -> str:
    before_set = _regimen_to_set(row["top_regimen_before"])
    after_set = _regimen_to_set(row["top_regimen_after"])
    added = sorted(after_set - before_set)
    removed = sorted(before_set - after_set)
    parts = []
    if added:
        parts.append(f"added: {', '.join(added)}")
    if removed:
        parts.append(f"removed: {', '.join(removed)}")
    if not parts:
        if row["hcos_before"] != row["hcos_after"]:
            return "same drugs, HCOS score changed"
        return "no change"
    return "; ".join(parts)


def diff_regimens_with_and_without_discovery(
    patient_barcodes: List[str],
    maf_glob_pattern: str,
    gene_panel: List[str],
    curated_gene_drugs: Dict[str, List[str]],
    merged_gene_drugs: Dict[str, List[str]],
    resolve_tp53_fn: Callable,
    drug_graph_cls,
    synergy_net_cls,
    hcos_fn: Callable,
    mdcoe_fn: Callable,
    max_patients: Optional[int] = None,
) -> pd.DataFrame:
    """Returns a DataFrame with one row per patient whose top_regimen or
    hcos changed, including hcos_delta (after - before) and
    reason_for_change (which specific drugs were added/removed, or 'same
    drugs, HCOS score changed' if only the score moved)."""
    logger.info(f"running before/after diff across {len(patient_barcodes)} patients")

    before = run_cohort_wide_analysis(
        patient_barcodes, maf_glob_pattern, gene_panel, curated_gene_drugs,
        resolve_tp53_fn, drug_graph_cls, synergy_net_cls, hcos_fn, mdcoe_fn,
        max_patients=max_patients,
    ).set_index("patient")
    after = run_cohort_wide_analysis(
        patient_barcodes, maf_glob_pattern, gene_panel, merged_gene_drugs,
        resolve_tp53_fn, drug_graph_cls, synergy_net_cls, hcos_fn, mdcoe_fn,
        max_patients=max_patients,
    ).set_index("patient")

    joined = before.join(after, lsuffix="_before", rsuffix="_after")
    changed = joined[
        (joined["top_regimen_before"] != joined["top_regimen_after"])
        | (joined["hcos_before"] != joined["hcos_after"])
    ].copy()

    if changed.empty:
        logger.info("no patients' top regimen or HCOS score changed")
        changed["hcos_delta"] = pd.Series(dtype=float)
        changed["reason_for_change"] = pd.Series(dtype=str)
        return changed.reset_index()

    changed["hcos_delta"] = changed["hcos_after"] - changed["hcos_before"]
    changed["reason_for_change"] = changed.apply(_reason_for_change, axis=1)
    logger.info(f"{len(changed)} patient(s) changed")
    return changed.reset_index()
