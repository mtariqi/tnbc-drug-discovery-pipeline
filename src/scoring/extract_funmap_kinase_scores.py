"""
extract_funmap_kinase_scores.py

Streams through the large (1.76 GB compressed, ~55 million row) FunMap
all-pairs prediction-score file from Zenodo (record 10080764,
all_pair_prediction_score.tsv.gz) WITHOUT loading it fully into memory,
and extracts only the rows where both genes are in the 90-kinase panel.

Real, confirmed file format (verified against actual output from the
collaborator's real download, not assumed):
    index                                       ALL_RNA    ALL_PRO   ALL_RNA_PRO
    ('ENSG00000000003', 'ENSG00000000005')      0.599...   0.251...  0.485...

The first column is a Python-tuple-formatted string containing two Ensembl
gene IDs (not gene symbols, and not two separate columns). There are three
score columns; ALL_RNA_PRO (the combined RNA+protein model) is used by
default, matching the model most likely underlying the published funmap.tsv
edge list, but this is a real design choice -- not verified against the
paper's own methods section -- and should be treated as such until confirmed.

Usage:
    python3 extract_funmap_kinase_scores.py \
        <path_to_all_pair_prediction_score.tsv.gz> \
        <path_to_ensembl_to_symbol_mapping.tsv> \
        <path_to_kinase_90_list.txt> \
        <output_path> \
        [--score-column ALL_RNA_PRO]

The Ensembl-to-symbol mapping file must be a two-column TSV (no header
required, but if present it will be skipped): ensembl_id, gene_symbol.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

TUPLE_PATTERN = re.compile(r"^\('([^']+)',\s*'([^']+)'\)$")


def load_kinase_ensembl_ids(mapping_path: str, kinase_panel_path: str):
    """Builds a symbol->Ensembl and Ensembl->symbol lookup restricted to the
    90-kinase panel, from a real, local Ensembl-to-symbol mapping file."""
    kinase_symbols = set()
    with open(kinase_panel_path) as f:
        for line in f:
            line = line.strip()
            if line:
                kinase_symbols.add(line.upper())

    ensembl_to_symbol = {}
    with open(mapping_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            ensembl_id, symbol = parts[0].strip(), parts[1].strip()
            # strip Ensembl version suffix if present (e.g., ENSG00000141736.13 -> ENSG00000141736)
            ensembl_id_base = ensembl_id.split(".")[0]
            if symbol.upper() in kinase_symbols:
                ensembl_to_symbol[ensembl_id_base] = symbol.upper()

    missing = kinase_symbols - set(ensembl_to_symbol.values())
    if missing:
        print(f"WARNING: {len(missing)} of {len(kinase_symbols)} kinase panel genes "
              f"were not found in the mapping file: {sorted(missing)}", file=sys.stderr)

    print(f"Resolved {len(ensembl_to_symbol)} of {len(kinase_symbols)} kinase panel genes "
          f"to Ensembl IDs.", file=sys.stderr)
    return ensembl_to_symbol


def stream_filter(gz_path: str, ensembl_to_symbol: dict, score_column: str, output_path: str):
    ensembl_ids = set(ensembl_to_symbol.keys())
    rows_written = 0
    rows_scanned = 0
    header_cols = None

    with gzip.open(gz_path, "rt") as fin, open(output_path, "w") as fout:
        fout.write("gene1\tgene2\tscore\n")
        for line in fin:
            rows_scanned += 1
            if header_cols is None:
                header_cols = line.rstrip("\n").split("\t")
                try:
                    score_idx = header_cols.index(score_column)
                except ValueError:
                    raise SystemExit(
                        f"Score column '{score_column}' not found in header: {header_cols}"
                    )
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) <= score_idx:
                continue

            m = TUPLE_PATTERN.match(parts[0])
            if not m:
                continue
            g1, g2 = m.group(1), m.group(2)

            if g1 in ensembl_ids and g2 in ensembl_ids:
                sym1 = ensembl_to_symbol[g1]
                sym2 = ensembl_to_symbol[g2]
                score = parts[score_idx]
                fout.write(f"{sym1}\t{sym2}\t{score}\n")
                rows_written += 1

            if rows_scanned % 5_000_000 == 0:
                print(f"  ...scanned {rows_scanned:,} rows, matched {rows_written} so far", file=sys.stderr)

    print(f"Done. Scanned {rows_scanned:,} total rows, wrote {rows_written} kinase-panel pairs "
          f"to {output_path}.", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gz_path")
    ap.add_argument("mapping_path")
    ap.add_argument("kinase_panel_path")
    ap.add_argument("output_path")
    ap.add_argument("--score-column", default="ALL_RNA_PRO")
    args = ap.parse_args()

    ensembl_to_symbol = load_kinase_ensembl_ids(args.mapping_path, args.kinase_panel_path)
    if len(ensembl_to_symbol) < 2:
        raise SystemExit("Fewer than 2 kinase genes resolved to Ensembl IDs -- check the mapping file.")

    stream_filter(args.gz_path, ensembl_to_symbol, args.score_column, args.output_path)


if __name__ == "__main__":
    main()
