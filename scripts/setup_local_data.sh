#!/usr/bin/env bash
# Run this ON YOUR OWN MACHINE, from the repo root, after editing the SOURCE_* paths below
# to point at wherever you already downloaded TCGA-BRCA / DepMap / STRING / DGIdb data.
# It copies (not moves) your existing files into data/raw/ in the layout data/README.md
# documents, so nothing you already have gets deleted or overwritten by accident.

set -euo pipefail

# --- EDIT THESE FOUR PATHS -------------------------------------------------
SOURCE_TCGA_DIR="/path/to/your/tcga_brca_downloads"
SOURCE_DEPMAP_DIR="/path/to/your/depmap_downloads"
SOURCE_STRING_FILE="/path/to/your/string_network_cache.tsv"
SOURCE_DGIDB_FILE="/path/to/your/dgidb_interactions.tsv"
# ---------------------------------------------------------------------------

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "Copying TCGA-BRCA files into data/raw/tcga_brca/ ..."
mkdir -p "$REPO_ROOT/data/raw/tcga_brca"
cp -n "$SOURCE_TCGA_DIR"/*.tsv "$REPO_ROOT/data/raw/tcga_brca/" 2>/dev/null || true
cp -n "$SOURCE_TCGA_DIR"/*.csv "$REPO_ROOT/data/raw/tcga_brca/" 2>/dev/null || true
cp -rn "$SOURCE_TCGA_DIR"/maf "$REPO_ROOT/data/raw/tcga_brca/" 2>/dev/null || true

echo "Copying DepMap files into data/raw/depmap/ ..."
mkdir -p "$REPO_ROOT/data/raw/depmap"
cp -n "$SOURCE_DEPMAP_DIR"/*.csv "$REPO_ROOT/data/raw/depmap/" 2>/dev/null || true

echo "Copying STRING cache into data/raw/string/ ..."
mkdir -p "$REPO_ROOT/data/raw/string"
cp -n "$SOURCE_STRING_FILE" "$REPO_ROOT/data/raw/string/" 2>/dev/null || true

echo "Copying DGIdb cache into data/raw/dgidb/ ..."
mkdir -p "$REPO_ROOT/data/raw/dgidb"
cp -n "$SOURCE_DGIDB_FILE" "$REPO_ROOT/data/raw/dgidb/" 2>/dev/null || true

echo ""
echo "Done. Nothing was overwritten (cp -n skips files that already exist at the destination)."
echo "Check data/README.md for the exact expected filenames if anything above found nothing to copy."
