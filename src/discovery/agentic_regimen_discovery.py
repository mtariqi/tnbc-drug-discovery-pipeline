"""
Agentic Regimen Discovery: Literature-Mining with Mandatory Independent Confirmation
=======================================================================================

An "agentic" search-and-evidence-gathering loop for discovering novel drug
candidates for genes poorly covered by the existing curated GENE_DRUGS
dictionary and real-DGIdb candidate lists (see build_real_gene_drugs.py,
compare_curated_vs_real_gene_drugs.py).

WHAT "AGENTIC" MEANS HERE, PRECISELY:
    This automates the SEARCH AND EVIDENCE-GATHERING loop -- deciding what
    to search for based on coverage gaps, retrieving real literature,
    extracting candidates, and independently verifying them -- NOT free-form
    biological reasoning. No mechanism, target, or drug name is ever
    accepted on the strength of extracted text alone. Every candidate must
    pass an INDEPENDENT confirmation step (a live DGIdb query proving the
    drug-gene interaction is actually registered in a curated database)
    before it is treated as real. This mirrors the same discipline already
    applied throughout this project to DGIdb's own noise (lab codes,
    literature citations, gene-names-as-drugs).

HARD LIMITS, STATED EXPLICITLY:
    - This does NOT replace expert review. Every surviving candidate must
      still be read by a human (ideally your supervisor) in its original
      abstract context before being treated as a serious hypothesis --
      automated extraction can pull a drug mentioned as a NEGATIVE
      comparator or an unrelated finding, not necessarily a genuine
      proposed combination partner.
    - This does NOT invent novel biological mechanisms. It finds real
      papers that already propose a combination and extracts what they
      already say, with a citation. It cannot and does not reason beyond
      that.
    - Final regimen ranking still uses the existing, validated MDCOE/HCOS
      scoring -- this only feeds new, provenance-tagged candidates into
      that unchanged pipeline.

REAL, VERIFIED API DETAILS (checked directly, not from memory):
    - PubMed ESearch: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
      params: db=pubmed, term=<query>, retmode=json, retmax=<n>
      -> JSON response, PMIDs at data['esearchresult']['idlist']
    - PubMed EFetch:  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
      params: db=pubmed, id=<comma-separated PMIDs>, retmode=text, rettype=abstract
      -> plain text abstracts
    - No API key required for reasonable use (NCBI recommends one for
      higher rate limits: 3 req/sec without a key, 10/sec with one --
      see https://www.ncbi.nlm.nih.gov/books/NBK25497/ if running this at
      volume).
"""

from __future__ import annotations

import re
import time
from typing import Dict, List, Optional

import pandas as pd
import requests

PUBMED_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Reused exactly from the validated drug-name filter (build_real_gene_drugs.py)
INN_SUFFIXES = (
    "nib", "rafenib", "tinib", "parib", "ciclib", "mab", "zumab", "ximab", "tuzumab",
    "umab", "stat", "degib", "lisib", "setron", "sertib", "clax", "dostat", "metinib",
    "afenib", "bulin", "gene", "kinra", "cept", "mustine", "platin", "rubicin",
    "apopt",  # p53 reactivators, e.g. eprenetapopt -- added after a real run missed it
    # NOTE: p53-reactivator/degrader compounds are still often code-named (PC14586,
    # COTI-2) rather than following an established INN suffix -- suffix matching
    # cannot catch these regardless of how many endings are added. This is a real,
    # acknowledged ceiling on this extraction method, not something a longer suffix
    # list fully solves.
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

    resp = requests.get(PUBMED_ESEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    pmids = resp.json().get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return []

    time.sleep(0.34)  # stay under the 3 req/sec unauthenticated rate limit
    fetch_params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "text", "rettype": "abstract"}
    if api_key:
        fetch_params["api_key"] = api_key
    fetch_resp = requests.get(PUBMED_EFETCH_URL, params=fetch_params, timeout=30)
    fetch_resp.raise_for_status()

    # EFetch returns all abstracts concatenated as plain text, separated by
    # blank lines and numbered entries -- split on the numbered-entry pattern.
    raw = fetch_resp.text
    entries = re.split(r"\n\d+\.\s", "\n" + raw)[1:]  # drop empty first split
    results = []
    for pmid, entry in zip(pmids, entries):
        results.append({"pmid": pmid, "abstract": entry.strip()})
    return results


# =====================================================================
# 3. EXTRACT CANDIDATES (conservative: real INN-suffix pattern only)
# =====================================================================

EXTRACTION_STOPLIST = {
    "gene", "proto-oncogene", "oncogene",  # bare/compound generic biology terms ending in a real gene-therapy suffix
    "celgene", "genmab", "alphamab", "numab",  # company/product-line names, not drugs, that happen to end in a real suffix
}


def extract_candidate_drugs_from_abstract(abstract: str) -> List[str]:
    """
    Conservative extraction: only tokens matching real INN drug-naming
    suffixes are extracted as candidates -- NEVER any capitalized word or
    generic noun phrase. This deliberately trades recall for precision:
    it will miss some real drugs with unusual names, but it will not
    hallucinate plausible-sounding fake ones, which is the more important
    property for a discovery pipeline whose output feeds a real regimen
    ranking.

    REAL BUG FOUND AND FIXED (see chat): a live run extracted "gene",
    "proto-oncogene", and several company names ("celgene", "genmab") as
    candidates, because the suffix check was `word.endswith(suffix)` --
    and the bare word "gene" trivially satisfies endswith("gene") since it
    IS that suffix, with nothing in front of it. Every one of these was
    correctly rejected downstream by the mandatory DGIdb confirmation
    step (none had dgidb_confirmed=True), so nothing false ever reached
    scoring -- but they added pure noise to the transparency table. Now
    requires at least 2 real characters before the suffix, plus a small
    explicit stoplist for the exact terms this real run surfaced.
    """
    words = re.findall(r"\b[a-zA-Z][a-zA-Z0-9\-]{3,}\b", abstract.lower())
    candidates = [
        w for w in words
        if w.endswith(INN_SUFFIXES)
        and w not in EXTRACTION_STOPLIST
        and any(len(w) >= len(suffix) + 2 for suffix in INN_SUFFIXES if w.endswith(suffix))
    ]
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
    corroborates it. Reuses the same DGIdb v5 GraphQL endpoint already
    validated in kinase_data_fetchers.py.
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
    return candidate_drug.lower() in known_drugs


# =====================================================================
# 4b. MENTION-CONTEXT ASSESSMENT (real gap found this session, see chat)
# =====================================================================

# Real, concrete finding this addresses: a first live run extracted
# "olaparib" from an abstract whose actual point was overcoming olaparib
# RESISTANCE, and "pembrolizumab" from an abstract explicitly describing
# its combination trial as showing "limited efficacy" -- both drugs are
# real, both were independently DGIdb-confirmed, but both were being
# recommended backwards relative to what the source paper actually says.
# This is a lightweight, pattern-based flag -- NOT free-form reasoning --
# to catch that specific, now-demonstrated failure mode. It does not
# replace the mandatory human read of the full abstract; it prioritizes
# what to look at with suspicion first.

CAUTIONARY_PATTERNS = (
    r"resistan\w*\s+to\s+", r"\bresistant\b", r"cross-resistan\w*",
    r"limited efficacy", r"failed to improve", r"did not improve",
    r"lack of\s+\w*\s*efficacy", r"no improvement", r"acquired resistance",
)
SUPPORTIVE_PATTERNS = (
    r"synerg\w*", r"combin\w*(?:\s+\w+){0,3}\s+with\b", r"combination of",
    r"enhanc\w+\s+cytotoxicity", r"potentiat\w*", r"promising",
)


def assess_mention_context(abstract: str, candidate_drug: str) -> str:
    """
    Classifies each SENTENCE that mentions candidate_drug as cautionary
    (resistance/limited-efficacy framing) or supportive (synergy/
    combination framing), then compares total counts across all such
    sentences. Returns 'cautionary', 'supportive', or 'unclear' -- a
    triage label, not a verdict; 'unclear' does not mean safe, it means no
    strong signal was found and a human read is still required regardless
    of label.

    THREE real bugs found and fixed while testing against actual PubMed
    abstracts, not synthetic text -- kept here as a record of why the
    design looks like this:
    (1) An earlier fixed-character-window version let a single cautionary
        phrase anywhere near ANY mention override clear, later synergy
        language for the same drug in the same abstract.
    (2) A drug mentioned several times close together produced heavily
        OVERLAPPING character windows, multiply-counting the same
        surrounding phrase and skewing the ratio.
    (3) Even after fixing (1) and (2) by merging windows, a fixed character
        radius in dense, multi-drug scientific text still pulled in
        sentences about a DIFFERENT drug entirely (two drugs discussed in
        adjacent sentences bled into each other's context). Switching to
        whole-sentence attribution -- the natural unit of a single claim
        -- fixes this: a sentence not mentioning the drug is never counted
        for it, regardless of proximity.

    This remains a pattern-based heuristic, not language understanding --
    a sentence discussing two drugs together can still misattribute a
    claim between them. It is a triage aid, not a substitute for reading
    the source.
    """
    drug_lower = candidate_drug.lower()
    sentences = re.split(r"(?<=[.!?])\s+", abstract)

    cautionary_count = 0
    supportive_count = 0
    found = False
    for sent in sentences:
        sent_lower = sent.lower()
        if drug_lower not in sent_lower:
            continue
        found = True
        cautionary_count += sum(1 for p in CAUTIONARY_PATTERNS if re.search(p, sent_lower))
        supportive_count += sum(1 for p in SUPPORTIVE_PATTERNS if re.search(p, sent_lower))

    if not found or (cautionary_count == 0 and supportive_count == 0):
        return "unclear"
    if supportive_count > cautionary_count:
        return "supportive"
    return "cautionary"  # strictly more cautionary hits, or a tie -- safer default either way


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
    min_candidates: int = 2,
) -> pd.DataFrame:
    """
    Full loop: identify gaps -> search literature -> extract candidates ->
    confirm independently -> return a fully provenance-tagged table.
    Every row traces to a real PMID and passed independent confirmation --
    nothing here is asserted without both.

    REAL BUG FOUND AND FIXED (see chat): this function used to call
    identify_coverage_gaps() with no min_candidates argument at all,
    silently using its default of 2 regardless of what a caller computed
    or passed elsewhere -- a manual identify_coverage_gaps(..., min_candidates=4)
    call showed 5 gap genes, but this function then recomputed gaps on its
    own with the hardcoded default and found 0. min_candidates is now a
    real parameter here, passed straight through.
    """
    gaps = identify_coverage_gaps(genes, curated_gene_drugs, real_gene_drugs, min_candidates=min_candidates)
    print(f"Coverage gaps identified ({len(gaps)}/{len(genes)} genes, min_candidates={min_candidates}): {gaps}")

    rows = []
    for gene in gaps:
        print(f"Searching PubMed for {gene} + {cancer_context}...")
        papers = search_pubmed_for_combination_therapy(gene, cancer_context, max_papers_per_gene, api_key)
        print(f"  {len(papers)} paper(s) found")

        for paper in papers:
            candidates = extract_candidate_drugs_from_abstract(paper["abstract"])
            for candidate in candidates:
                confirmed = confirm_candidate_via_dgidb(candidate, gene)
                context = assess_mention_context(paper["abstract"], candidate)
                rows.append({
                    "gene": gene,
                    "candidate_drug": candidate,
                    "pmid": paper["pmid"],
                    "dgidb_confirmed": confirmed,
                    "mention_context": context,
                    "abstract_snippet": paper["abstract"][:200],
                })
                time.sleep(0.1)  # light courtesy delay on DGIdb calls

    result = pd.DataFrame(rows)
    if result.empty:
        print("No candidates found for any gap gene.")
        return result

    n_confirmed = result["dgidb_confirmed"].sum()
    n_cautionary = (result["mention_context"] == "cautionary").sum()
    print(f"\n{len(result)} total (gene, candidate, paper) extractions; "
          f"{n_confirmed} passed independent DGIdb confirmation; "
          f"{n_cautionary} flagged 'cautionary' (mentioned near resistance/limited-efficacy language).")
    print("UNCONFIRMED candidates are NOT deleted from this table (kept for "
          "transparency) but must NOT be used in scoring -- filter to "
          "dgidb_confirmed == True before passing anything downstream. "
          "'cautionary' candidates need a full manual read before use, even "
          "if DGIdb-confirmed -- see mention_context column.")
    return result


def merge_into_candidate_pool(
    discovery_result: pd.DataFrame,
    existing_gene_drugs: Dict[str, List[str]],
    require_supportive: bool = True,
) -> Dict[str, List[str]]:
    """
    Merges ONLY dgidb_confirmed candidates into an existing gene->drugs
    dict (e.g. to feed into mdcoe.py's MDCOE search).

    REAL BUG FOUND AND FIXED (see chat): the original version only
    excluded candidates explicitly labeled 'cautionary', which silently
    let 'unclear' candidates through unchanged -- 'unclear' means the
    heuristic found NO signal either way, not that the candidate is safe.
    A real run demonstrated this exactly: an olaparib mention describing
    baseline drug mechanism (no resistance or synergy language at all) was
    correctly labeled 'unclear' and then merged in anyway, undetected,
    because 'unclear' was treated as equivalent to 'supportive'.

    By default (require_supportive=True), only candidates explicitly
    labeled 'supportive' are merged -- 'cautionary' AND 'unclear' both now
    require the same manual override (require_supportive=False) after you
    have personally read the source abstract and confirmed the candidate
    is genuinely being recommended.
    """
    merged = {k: list(v) for k, v in existing_gene_drugs.items()}
    if discovery_result.empty:
        return merged
    survivors = discovery_result[discovery_result["dgidb_confirmed"]]
    if require_supportive:
        n_before = len(survivors)
        survivors = survivors[survivors["mention_context"] == "supportive"]
        n_excluded = n_before - len(survivors)
        if n_excluded:
            print(f"Excluded {n_excluded} DGIdb-confirmed candidate(s) not labeled 'supportive' "
                  f"(cautionary OR unclear) from the merge -- review them manually and pass "
                  f"require_supportive=False if you've confirmed a specific one is genuinely "
                  f"recommended, not being worked around or merely mentioned in passing.")
    for _, row in survivors.iterrows():
        merged.setdefault(row["gene"], [])
        if row["candidate_drug"] not in merged[row["gene"]]:
            merged[row["gene"]].append(row["candidate_drug"])
    return merged


# =====================================================================
# SMOKE TESTS -- mocked PubMed/DGIdb responses, no live network needed
# =====================================================================

def _run_smoke_tests():
    from unittest.mock import patch, MagicMock

    print("=== Testing extract_candidate_drugs_from_abstract() ===")
    fake_abstract = (
        "In this study, we evaluated the combination of alpelisib and "
        "trastuzumab in PIK3CA-mutant, HER2-low breast cancer models. "
        "Patients were also treated with standard chemotherapy and pain "
        "medication. The combination showed synergistic reduction in tumor growth."
    )
    candidates = extract_candidate_drugs_from_abstract(fake_abstract)
    print(f"  Extracted: {candidates}")
    assert "alpelisib" in candidates
    assert "trastuzumab" in candidates
    assert "chemotherapy" not in candidates, "generic terms must never be extracted as drug candidates"
    print("  PASSED: real INN-suffix drugs extracted, generic terms correctly excluded.\n")

    print("=== Testing search_pubmed_for_combination_therapy() with mocked E-utilities ===")
    fake_esearch_response = {"esearchresult": {"idlist": ["12345678", "87654321"]}}
    fake_efetch_response = (
        "\n1. Journal of Oncology. 2023.\n"
        "Title: Alpelisib plus trastuzumab in TNBC.\n"
        "Alpelisib combined with trastuzumab showed activity in PIK3CA-altered TNBC models.\n\n"
        "2. Cancer Research. 2022.\n"
        "Title: PTEN-loss combination strategies.\n"
        "Capivasertib demonstrated efficacy when combined with standard therapy in PTEN-null models.\n"
    )
    with patch("requests.get") as mock_get:
        def side_effect(url, params=None, timeout=None):
            resp = MagicMock()
            resp.raise_for_status.return_value = None
            if "esearch" in url:
                resp.json.return_value = fake_esearch_response
            else:
                resp.text = fake_efetch_response
            return resp
        mock_get.side_effect = side_effect
        papers = search_pubmed_for_combination_therapy("PTEN", max_results=2)
    print(f"  {len(papers)} papers retrieved")
    assert len(papers) == 2
    assert papers[0]["pmid"] == "12345678"
    assert "alpelisib" in papers[0]["abstract"].lower()
    print("  PASSED: real E-utilities two-step flow (esearch -> efetch) correctly parsed.\n")

    print("=== Testing confirm_candidate_via_dgidb() with mocked GraphQL response ===")
    fake_dgidb_response = {
        "data": {"genes": {"nodes": [
            {"interactions": [{"drug": {"name": "ALPELISIB"}}, {"drug": {"name": "CAPIVASERTIB"}}]}
        ]}}
    }
    with patch("requests.post") as mock_post:
        resp = MagicMock()
        resp.raise_for_status.return_value = None
        resp.json.return_value = fake_dgidb_response
        mock_post.return_value = resp
        confirmed_true = confirm_candidate_via_dgidb("alpelisib", "PTEN")
        confirmed_false = confirm_candidate_via_dgidb("totally_made_up_drug", "PTEN")
    assert confirmed_true is True
    assert confirmed_false is False
    print("  PASSED: real DGIdb match correctly confirmed; unconfirmed candidate correctly rejected.\n")

    print("=== Testing identify_coverage_gaps() ===")
    curated = {"EGFR": ["afatinib", "erlotinib"], "PTEN": ["alpelisib"]}
    real = {"EGFR": ["afatinib"], "TP53": []}
    gaps = identify_coverage_gaps(["EGFR", "PTEN", "TP53"], curated, real, min_candidates=2)
    print(f"  Gaps found: {gaps}")
    assert "PTEN" in gaps  # only 1 candidate total
    assert "TP53" in gaps  # 0 candidates
    assert "EGFR" not in gaps  # 2 candidates, meets threshold
    print("  PASSED: correctly identifies under-covered genes as search targets.\n")

    print("=== Testing merge_into_candidate_pool() only merges CONFIRMED + explicitly SUPPORTIVE candidates ===")
    fake_discovery = pd.DataFrame([
        {"gene": "PTEN", "candidate_drug": "capivasertib", "pmid": "111", "dgidb_confirmed": True, "mention_context": "supportive", "abstract_snippet": "..."},
        {"gene": "PTEN", "candidate_drug": "fake_unconfirmed_drug", "pmid": "222", "dgidb_confirmed": False, "mention_context": "supportive", "abstract_snippet": "..."},
        {"gene": "PTEN", "candidate_drug": "resisted_drug", "pmid": "333", "dgidb_confirmed": True, "mention_context": "cautionary", "abstract_snippet": "..."},
        {"gene": "PTEN", "candidate_drug": "mentioned_in_passing_drug", "pmid": "444", "dgidb_confirmed": True, "mention_context": "unclear", "abstract_snippet": "..."},
    ])
    merged = merge_into_candidate_pool(fake_discovery, {"PTEN": ["alpelisib"]})
    print(f"  Merged PTEN candidates: {merged['PTEN']}")
    assert "capivasertib" in merged["PTEN"]
    assert "fake_unconfirmed_drug" not in merged["PTEN"], "unconfirmed candidates must never reach the scoring pipeline"
    assert "resisted_drug" not in merged["PTEN"], "cautionary-flagged candidates must be excluded by default"
    assert "mentioned_in_passing_drug" not in merged["PTEN"], "unclear-flagged candidates must ALSO be excluded by default -- this is the real bug just found and fixed"
    print("  PASSED: only explicitly DGIdb-confirmed AND supportive candidates merged; both cautionary AND unclear correctly excluded.\n")

    print("All offline smoke tests passed.")


if __name__ == "__main__":
    _run_smoke_tests()
