"""
MAF (Mutation Annotation Format) reading.

Extracts a single patient's altered genes (restricted to a gene panel)
from real, gzipped GDC MAF files. Uses the CONFIRMED real column names
(Hugo_Symbol, Tumor_Sample_Barcode) verified against an actual downloaded
TCGA MAF file -- not assumed from documentation alone.
"""

from __future__ import annotations

import gzip
import glob
import logging
from typing import List

logger = logging.getLogger(__name__)


def extract_patient_altered_genes(patient_barcode: str, maf_glob_pattern: str, gene_panel: List[str]) -> List[str]:
    """
    Searches real, gzipped GDC MAF files matching maf_glob_pattern for
    this patient's barcode and returns which genes in gene_panel show a
    variant for them.
    """
    altered = set()
    matched_files = glob.glob(maf_glob_pattern)
    if not matched_files:
        logger.warning(f"no MAF files matched glob pattern {maf_glob_pattern!r}")

    for maf_path in matched_files:
        opener = gzip.open if maf_path.endswith(".gz") else open
        with opener(maf_path, "rt") as f:
            header = None
            for line in f:
                if line.startswith("Hugo_Symbol"):
                    header = line.rstrip("\n").split("\t")
                    continue
                if header is None:
                    continue
                if patient_barcode not in line:
                    continue
                fields = line.rstrip("\n").split("\t")
                row = dict(zip(header, fields))
                gene = row.get("Hugo_Symbol", "")
                barcode_field = row.get("Tumor_Sample_Barcode", "")
                if patient_barcode in barcode_field and gene in gene_panel:
                    altered.add(gene)

    return sorted(altered)
