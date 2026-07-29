"""
Parallel literature discovery across genes.

Runs coverage-gap identification once (fast, local), then fires one
run_discovery_loop() call PER GAP GENE concurrently via
ThreadPoolExecutor -- network-bound PubMed/DGIdb calls parallelize well
since threads release the GIL during I/O wait.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List

import pandas as pd

from . import agentic_regimen_discovery as ard

logger = logging.getLogger(__name__)

_EMPTY_COLUMNS = ["gene", "candidate_drug", "pmid", "dgidb_confirmed", "abstract_snippet"]


def parallel_discovery(
    genes: List[str],
    curated_gene_drugs: Dict[str, List[str]],
    real_gene_drugs: Dict[str, List[str]],
    max_workers: int = 8,
    **kwargs,
) -> pd.DataFrame:
    """
    Drop-in replacement for calling run_discovery_loop() directly: same
    return schema, but each gap gene's search+confirm loop runs in its
    own thread. A failure for one gene is logged and excluded rather than
    aborting the whole run.
    """
    gaps = ard.identify_coverage_gaps(genes, curated_gene_drugs, real_gene_drugs)
    logger.info(f"coverage gaps identified: {len(gaps)}/{len(genes)} genes -> {gaps}")

    if not gaps:
        return pd.DataFrame(columns=_EMPTY_COLUMNS)

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            gene: executor.submit(ard.run_discovery_loop, [gene], curated_gene_drugs, real_gene_drugs, **kwargs)
            for gene in gaps
        }
        for gene, future in futures.items():
            try:
                results.append(future.result())
            except Exception:
                logger.exception(f"discovery failed for gene={gene}; excluding it from this run's results")

    if not results:
        return pd.DataFrame(columns=_EMPTY_COLUMNS)
    return pd.concat(results, ignore_index=True)
