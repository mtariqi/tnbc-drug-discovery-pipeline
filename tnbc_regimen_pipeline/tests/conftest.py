"""
Shared fixtures for tnbc_regimen_pipeline's test suite.

All network-touching functions (PubMed E-utilities, DGIdb GraphQL) are
mocked here -- this sandbox has no network access, and even with network
access, hitting live PubMed/DGIdb on every test run would make tests
slow, flaky, and dependent on literature that changes over time. What IS
real: MAF file parsing uses an actual gzipped file on disk (tmp_path),
not a mock, since that logic is worth testing against real I/O.
"""

from __future__ import annotations

import gzip
import itertools
from unittest.mock import MagicMock

import pytest


# ---------------------------------------------------------------------
# Fake MDCOE/HCOS components (stand-ins for the real, separately-
# validated scoring pipeline that lives outside this package)
# ---------------------------------------------------------------------

class FakeDrugGraph:
    def __init__(self, drugs):
        self.drugs = drugs


class FakeSynergyNet:
    pass


@pytest.fixture
def fake_hcos():
    def _hcos(regimen, net):
        return len(regimen) * 0.1
    return _hcos


@pytest.fixture
def fake_mdcoe():
    def _mdcoe(graph, net, hcos_fn, beam_width, max_depth, top_k):
        best = None
        for size in range(2, min(len(graph.drugs), max_depth) + 1):
            for combo in itertools.combinations(graph.drugs, size):
                score = hcos_fn(combo, net)
                if best is None or score > best[1]:
                    best = (combo, score)
        return [best] if best else []
    return _mdcoe


@pytest.fixture
def fake_resolve_tp53():
    def _resolve(alt_type):
        return ["adavosertib", "venetoclax"]
    return _resolve


@pytest.fixture
def fake_drug_graph_cls():
    return FakeDrugGraph


@pytest.fixture
def fake_synergy_net_cls():
    return FakeSynergyNet


# ---------------------------------------------------------------------
# Gene/drug dictionaries and patient data
# ---------------------------------------------------------------------

@pytest.fixture
def curated_gene_drugs():
    return {
        "EGFR": ["afatinib", "erlotinib"],
        "PTEN": ["alpelisib"],
        "TYK2": ["deucravacitinib", "baricitinib"],
    }


@pytest.fixture
def real_gene_drugs():
    return {"EGFR": ["afatinib"], "PTEN": [], "TYK2": ["deucravacitinib"]}


@pytest.fixture
def fake_patient_genes():
    return {"PATIENT-001": ["EGFR", "PTEN"], "PATIENT-002": ["EGFR"]}


# ---------------------------------------------------------------------
# Mocked PubMed/DGIdb responses
# ---------------------------------------------------------------------

@pytest.fixture
def fake_esearch_response():
    return {"esearchresult": {"idlist": ["999"]}}


@pytest.fixture
def fake_efetch_response():
    return (
        "\n1. Fake Journal. 2024.\n"
        "Capivasertib combined with alpelisib showed activity in PTEN-null TNBC models.\n"
    )


@pytest.fixture
def fake_dgidb_response():
    return {"data": {"genes": {"nodes": [{"interactions": [{"drug": {"name": "CAPIVASERTIB"}}]}]}}}


@pytest.fixture
def mock_requests_get(fake_esearch_response, fake_efetch_response):
    """A requests.get side_effect distinguishing PubMed ESearch vs EFetch by URL."""
    def _side_effect(url, params=None, timeout=None):
        resp = MagicMock()
        resp.raise_for_status.return_value = None
        if "esearch" in url:
            resp.json.return_value = fake_esearch_response
        else:
            resp.text = fake_efetch_response
        return resp
    return _side_effect


@pytest.fixture
def mock_requests_post(fake_dgidb_response):
    resp = MagicMock()
    resp.raise_for_status.return_value = None
    resp.json.return_value = fake_dgidb_response
    return MagicMock(return_value=resp)


# ---------------------------------------------------------------------
# Real (not mocked) gzipped MAF file on disk
# ---------------------------------------------------------------------

@pytest.fixture
def real_maf_file(tmp_path, fake_patient_genes):
    """Writes an actual gzipped MAF file matching fake_patient_genes, so
    MAF-parsing tests exercise real gzip + tab-delimited I/O rather than
    a mock of it."""
    maf_path = tmp_path / "cohort.maf.gz"
    header = "Hugo_Symbol\tVariant_Classification\tTumor_Sample_Barcode\n"
    rows = []
    for patient, genes in fake_patient_genes.items():
        for gene in genes:
            rows.append(f"{gene}\tMissense_Mutation\t{patient}-TUMOR\n")
    with gzip.open(maf_path, "wt") as f:
        f.write(header)
        for row in rows:
            f.write(row)
    return str(maf_path)
