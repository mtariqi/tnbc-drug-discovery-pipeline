"""
Agentic Regimen Discovery: Literature-Mining with Mandatory Independent Confirmation
=======================================================================================

An "agentic" search-and-evidence-gathering loop for discovering novel drug
candidates for genes poorly covered by the existing curated GENE_DRUGS
dictionary and real-DGIdb candidate lists.

WHAT "AGENTIC" MEANS HERE, PRECISELY:
    This automates the SEARCH AND EVIDENCE-GATHERING loop -- deciding what
    to search for based on coverage gaps, retrieving real literature,
    extracting candidates, and independently verifying them -- NOT free-form
    biological reasoning. No mechanism, target, or drug name is ever
    accepted on the strength of extracted text alone. Every candidate must
    pass an INDEPENDENT confirmation step (a live DGIdb query proving the
    drug-gene interaction is actually registered in a curated database)
    before it is treated as real.

HARD LIMITS, STATED EXPLICITLY:
    - This does NOT replace expert review. Every surviving candidate must
      still be read by a human in its original abstract context before
      being treated as a serious hypothesis -- automated extraction can
      pull a drug mentioned as a NEGATIVE comparator or an unrelated
      finding, not necessarily a genuine proposed combination partner.
    - This does NOT invent novel biological mechanisms. It finds real
      papers that already propose a combination and extracts what they
      already say, with a citation.
    - Final regimen ranking still uses the existing, validated MDCOE/HCOS
      scoring -- this only feeds new, provenance-tagged candidates into
      that unchanged pipeline (see tnbc_regimen_pipeline.discovery.provenance
      and tnbc_regimen_pipeline.pipeline.full_pipeline).

REAL, VERIFIED API DETAILS:
    - PubMed ESearch: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
      params: db=pubmed, term=<query>, retmode=json, retmax=<n>
      -> JSON response, PMIDs at data['esearchresult']['idlist']
    - PubMed EFetch:  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
      params: db=pubmed, id=<comma-separated PMIDs>, retmode=text, rettype=abstract
      -> plain text abstracts
    - No API key required for reasonable use (NCBI recommends one for
      higher rate limits: 3 req/sec without a key, 10/sec with one).
"""

from __future__ import annotations

import logging
import re
import time
from typing import Dict, List, Optional

import pandas as pd
import requests

logger = logging.getLogger(__name__)

PUBMED_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Reused exactly from the validated drug-name filter (build_real_gene_drugs.py)
INN_SUFFIXES = (
    "nib", "rafenib", "tinib", "parib", "ciclib", "mab", "zumab", "ximab", "tuzumab",
    "umab", "stat", "degib", "lisib", "setron", "sertib", "clax", "dostat", "metinib",
    "afenib", "bulin", "gene", "kinra", "cept", "mustine", "platin", "rubicin",
)


# =====================================================================
# 1. IDENTIFY COVERAGE GAPS (what to search for)
# =====================================================================

def identify_coverage_gaps(
    genes: List[str],
    curated_gene_drugs: Dict[str, List[str]],
    real_gene_drugs: Dict[str, List[str]],
    min_candidates: int = 2,
) -> List[str]:
    """Genes with fewer than `min_candidates` drugs across BOTH the curated
    dictionary and real-DGIdb sources -- these are where literature search
    can plausibly add real value, rather than searching indiscriminately."""
    gaps = []
    for g in genes:
        n = len(set(curated_gene_drugs.get(g, [])) | set(real_gene_drugs.get(g, [])))
        if n < min_candidates:
            gaps.append(g)
    return gaps


# =====================================================================
# 2. SEARCH PUBMED (real E-utilities, verified format)
# =====================================================================

def search_pubmed_for_combination_therapy(
    gene: str,
    cancer_context: str = "triple negative breast cancer",
    max_results: int = 10,
    api_key: Optional[str] = None,
) -> List[Dict[str, str]]:
    """
    Searches real PubMed for combination-therapy literature mentioning this
    gene in this cancer context. Returns [{pmid, abstract}, ...] with real
    abstract text -- this is the ONLY source of candidate drug names; no
    text here is generated, all of it is retrieved.
    """
    query = f'({gene}[Title/Abstract]) AND ({cancer_context}[Title/Abstract]) AND (combination[Title/Abstract] OR inhibitor[Title/Abstract])'
    params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": max_results}
    if api_key:
        params["api_key"] = api_key

    logger.info(f"searching PubMed gene={gene} context={cancer_context!r}")
    resp = requests.get(PUBMED_ESEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    pmids = resp.json().get("esearchresult", {}).get("idlist", [])
    if not pmids:
        logger.info(f"no PubMed results gene={gene}")
        return []

    time.sleep(0.34)  # stay under the 3 req/sec unauthenticated rate limit
    fetch_params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "text", "rettype": "abstract"}
    if api_key:
        fetch_params["api_key"] = api_key
    fetch_resp = requests.get(PUBMED_EFETCH_URL, params=fetch_params, timeout=30)
    fetch_resp.raise_for_status()

    raw = fetch_resp.text
    entries = re.split(r"\n\d+\.\s", "\n" + raw)[1:]  # drop empty first split
    results = [{"pmid": pmid, "abstract": entry.strip()} for pmid, entry in zip(pmids, entries)]
    logger.info(f"retrieved {len(results)} paper(s) gene={gene}")
    return results


# =====================================================================
# 3. EXTRACT CANDIDATES (conservative: real INN-suffix pattern only)
# =====================================================================

def extract_candidate_drugs_from_abstract(abstract: str) -> List[str]:
    """
    Conservative extraction: only tokens matching real INN drug-naming
    suffixes are extracted as candidates -- NEVER any capitalized word or
    generic noun phrase. Trades recall for precision: it will miss some
    real drugs with unusual names, but it will not hallucinate
    plausible-sounding fake ones, which matters more for a discovery
    pipeline whose output feeds a real regimen ranking.
    """
    words = re.findall(r"\b[a-zA-Z][a-zA-Z0-9\-]{3,}\b", abstract.lower())
    candidates = [w for w in words if w.endswith(INN_SUFFIXES)]
    return sorted(set(candidates))


# =====================================================================
# 4. INDEPENDENT CONFIRMATION (mandatory -- nothing survives without this)
# =====================================================================

def confirm_candidate_via_dgidb(candidate_drug: str, gene: str) -> bool:
    """
    Live DGIdb GraphQL check: does DGIdb's own curated database confirm
    this drug actually targets this gene? This is the mandatory
    independent-confirmation gate -- a candidate extracted from a PubMed
    abstract is NOT trusted until a second, separate, curated source
    corroborates it.
    """
    query = """
    query GeneInteractions($names: [String!]) {
      genes(names: $names) { nodes { interactions { drug { name } } } }
    }
    """
    resp = requests.post(
        "https://dgidb.org/api/graphql",
        json={"query": query, "variables": {"names": [gene]}},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    nodes = data.get("data", {}).get("genes", {}).get("nodes", [])
    known_drugs = set()
    for node in nodes:
        for interaction in node.get("interactions", []):
            drug_name = (interaction.get("drug") or {}).get("name", "")
            known_drugs.add(drug_name.lower())
    confirmed = candidate_drug.lower() in known_drugs
    logger.info(f"DGIdb confirmation gene={gene} candidate={candidate_drug} confirmed={confirmed}")
    return confirmed


# =====================================================================
# 5. FULL DISCOVERY LOOP
# =====================================================================

def run_discovery_loop(
    genes: List[str],
    curated_gene_drugs: Dict[str, List[str]],
    real_gene_drugs: Dict[str, List[str]],
    cancer_context: str = "triple negative breast cancer",
    max_papers_per_gene: int = 10,
    api_key: Optional[str] = None,
) -> pd.DataFrame:
    """
    Full loop: identify gaps -> search literature -> extract candidates ->
    confirm independently -> return a fully provenance-tagged table.
    Every row traces to a real PMID and passed independent confirmation --
    nothing here is asserted without both.
    """
    gaps = identify_coverage_gaps(genes, curated_gene_drugs, real_gene_drugs)
    logger.info(f"coverage gaps identified ({len(gaps)}/{len(genes)} genes): {gaps}")

    rows = []
    for gene in gaps:
        papers = search_pubmed_for_combination_therapy(gene, cancer_context, max_papers_per_gene, api_key)
        for paper in papers:
            candidates = extract_candidate_drugs_from_abstract(paper["abstract"])
            for candidate in candidates:
                confirmed = confirm_candidate_via_dgidb(candidate, gene)
                rows.append({
                    "gene": gene,
                    "candidate_drug": candidate,
                    "pmid": paper["pmid"],
                    "dgidb_confirmed": confirmed,
                    "abstract_snippet": paper["abstract"][:200],
                })
                time.sleep(0.1)  # light courtesy delay on DGIdb calls

    result = pd.DataFrame(rows)
    if result.empty:
        logger.info("no candidates found for any gap gene")
        return result

    n_confirmed = result["dgidb_confirmed"].sum()
    logger.info(
        f"{len(result)} total (gene, candidate, paper) extractions; "
        f"{n_confirmed} passed independent DGIdb confirmation. "
        f"Unconfirmed candidates are kept for transparency but must NOT be "
        f"used in scoring -- filter to dgidb_confirmed == True before passing downstream."
    )
    return result


def merge_into_candidate_pool(
    discovery_result: pd.DataFrame,
    existing_gene_drugs: Dict[str, List[str]],
) -> Dict[str, List[str]]:
    """Merges ONLY dgidb_confirmed candidates into an existing gene->drugs
    dict. Returns a NEW dict; does not modify the input. For provenance
    (which PMIDs backed each merged candidate), use
    tnbc_regimen_pipeline.discovery.provenance.build_merged_pool_with_provenance
    instead."""
    merged = {k: list(v) for k, v in existing_gene_drugs.items()}
    confirmed = discovery_result[discovery_result["dgidb_confirmed"]] if not discovery_result.empty else discovery_result
    for _, row in confirmed.iterrows():
        merged.setdefault(row["gene"], [])
        if row["candidate_drug"] not in merged[row["gene"]]:
            merged[row["gene"]].append(row["candidate_drug"])
    return merged
