"""
extract_panel_mutations.py

Run this LOCALLY (not in Claude's sandbox -- it needs your real MAF files on disk).
Parses every GDC masked-somatic-mutation MAF.gz file under the given glob pattern,
keeps only mutations in the real 90/102-gene RTK/NRTK panel (imported directly from
your own mdcoe.py -- not re-typed, so it can't drift out of sync with the real panel),
and writes one clean CSV: one row per (patient, altered gene), with the real MAF
Variant_Classification preserved (needed for mdcoe.py's resolve_tp53_drugs zygosity
logic).

Usage:
    python3 extract_panel_mutations.py \
        --maf-glob "/home/mtariq/rtk_nrtk_tnbc/data/raw/tcga_brca/TCGA-BRCA/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/*/*.maf.gz" \
        --mdcoe-path /home/mtariq/rtk_nrtk_tnbc/src/mdcoe.py \
        --out patient_panel_mutations.csv

Then upload the resulting patient_panel_mutations.csv back into the chat.
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import sys
from pathlib import Path

import pandas as pd


def load_gene_panel(mdcoe_path: str) -> set[str]:
    """Imports GENE_PATHWAYS directly from your real mdcoe.py file (given by path,
    not by package name, so this works regardless of your repo's import setup)."""
    spec = importlib.util.spec_from_file_location("mdcoe_real", mdcoe_path)
    mdcoe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mdcoe)
    return set(mdcoe.GENE_PATHWAYS.keys())


def parse_one_maf(path: str, gene_panel: set[str]) -> pd.DataFrame:
    """Reads one gzipped MAF file, skipping the '#'-prefixed header comment lines
    GDC MAFs always include before the real column header, and returns only rows
    whose Hugo_Symbol is in gene_panel."""
    df = pd.read_csv(
        path, sep="\t", comment="#", low_memory=False,
        usecols=lambda c: c in {"Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"},
    )
    if df.empty:
        return df
    df = df[df["Hugo_Symbol"].isin(gene_panel)].copy()
    if df.empty:
        return df
    # Tumor_Sample_Barcode is the full aliquot barcode (e.g.
    # TCGA-XX-XXXX-01A-11D-XXXX-XX); truncate to the 12-char patient barcode
    # (TCGA-XX-XXXX) so this joins cleanly against clinical/subtype tables.
    df["patient_barcode"] = df["Tumor_Sample_Barcode"].str[:12]
    return df[["patient_barcode", "Hugo_Symbol", "Variant_Classification"]].rename(
        columns={"Hugo_Symbol": "gene", "Variant_Classification": "alteration_type"}
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maf-glob", required=True, help="Glob pattern matching all .maf.gz files")
    ap.add_argument("--mdcoe-path", required=True, help="Path to your real mdcoe.py")
    ap.add_argument("--out", default="patient_panel_mutations.csv",
                     help="Output path. Point this at your repo's data/processed/ "
                          "directory explicitly (e.g. --out data/processed/patient_panel_mutations.csv) "
                          "if you want run_cohort_hcos.py to pick it up automatically.")
    args = ap.parse_args()

    gene_panel = load_gene_panel(args.mdcoe_path)
    print(f"Loaded {len(gene_panel)}-gene panel from {args.mdcoe_path}")

    files = sorted(glob.glob(args.maf_glob))
    print(f"Found {len(files)} MAF files matching pattern")
    if not files:
        print("No files matched -- check --maf-glob", file=sys.stderr)
        sys.exit(1)

    frames = []
    for i, f in enumerate(files, 1):
        try:
            frames.append(parse_one_maf(f, gene_panel))
        except Exception as e:
            print(f"  [skip] {f}: {e}", file=sys.stderr)
        if i % 100 == 0:
            print(f"  processed {i}/{len(files)} files...")

    result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["patient_barcode", "gene", "alteration_type"]
    )
    result = result.drop_duplicates()

    result.to_csv(args.out, index=False)
    print(f"\nWrote {len(result)} (patient, gene) alteration rows to {args.out}")
    print(f"Unique patients with at least one panel-gene alteration: {result['patient_barcode'].nunique()}")
    print("\nPer-gene alteration counts (top 15):")
    print(result["gene"].value_counts().head(15))


if __name__ == "__main__":
    main()
