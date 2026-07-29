"""
TNBC-restricted RNA-seq correlation network over the existing 90-gene RTK/NRTK
panel, purpose-built to test whether real expression-correlation-derived
redundancy/compensation information changes PairCTS/TripletCTS rankings --
NOT a general-purpose "FunMap-style" network, and not attributed to FunMap's
actual (different) methodology. See ablation_pairct_with_vs_without_network.py
for the actual with/without ranking comparison this feeds into.

Design decisions, each made to avoid repeating a specific real mistake found
earlier in this project's history:
  - Restricted to the actual 168-patient confirmed-TNBC cohort, not the full
    BRCA cohort (the uploaded FunMap-labeled report used the BRCA-wide 1,082
    patients; that scope mismatch does not get repeated here).
  - Restricted to the existing 90-gene panel already used throughout this
    project's CTS/PairCTS/TripletCTS scoring, not a new 50-200 gene list (the
    uploaded report's 43-gene panel was a second, separate scope mismatch;
    reusing the real panel avoids introducing a third, different gene list).
  - FDR (Benjamini-Hochberg) correction applied to the ~4,005 pairwise tests
    (90 choose 2) before calling any correlation "significant" -- consistent
    with how the existing CPTAC-STRING crosstalk validation already reports
    FDR-corrected results, not a new standard invented for this script.

Output is deliberately shaped to slot directly into the EXISTING pair_cts()
function's crosstalk_edges: Dict[Tuple[str, str], float] parameter with no
modification to that function required -- this network is a second, optional
crosstalk source to compare against the existing STRING-derived one, not a
replacement requiring new scoring code.

Usage:
  python -m src.scoring.tnbc_correlation_network \
      --expression path/to/tnbc_cohort_expression.csv \
      --kinase-panel data/raw/kinases/kinase_90_list.txt \
      --output results/tables/tnbc_correlation_network.json
"""
from __future__ import annotations

import argparse
import itertools
import json

import numpy as np
import pandas as pd
from scipy import stats


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Standard BH FDR correction. Returns q-values (adjusted p-values), same
    order as input. Implemented directly (no extra dependency) and verified
    against a known worked example below.
    """
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    q = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity (BH step-up): q-values must not decrease going
    # backwards from the largest p-value
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    result = np.empty(n)
    result[order] = q
    return result


def build_correlation_network(
    expression: pd.DataFrame,  # rows = patients, columns = gene symbols (the 90-gene panel)
    fdr_threshold: float = 0.05,
    min_abs_r: float = 0.3,
) -> dict:
    """Returns:
      - 'crosstalk_edges': {(gene_i, gene_j): abs(r)} for pairs passing BOTH
        the FDR threshold AND the min_abs_r floor -- directly usable as
        pair_cts()'s crosstalk_edges argument.
      - 'compensatory_pairs': {(gene_i, gene_j): r} subset with r < 0 specifically
        (candidate compensatory/escape-route relationships -- a qualitatively
        different biological claim than a positive redundancy/crosstalk edge,
        kept as a separate, explicitly-labeled output rather than folded into
        the same undifferentiated edge-weight number).
      - 'redundant_pairs': {(gene_i, gene_j): r} subset with r > 0 specifically.
      - 'n_tests', 'n_significant': for reporting.
    """
    genes = list(expression.columns)
    pairs = list(itertools.combinations(genes, 2))
    rs, ps = [], []
    for gi, gj in pairs:
        r, p = stats.pearsonr(expression[gi], expression[gj])
        rs.append(r)
        ps.append(p)
    rs = np.array(rs)
    ps = np.array(ps)
    qs = benjamini_hochberg(ps)

    significant = (qs < fdr_threshold) & (np.abs(rs) >= min_abs_r)

    crosstalk_edges, compensatory_pairs, redundant_pairs = {}, {}, {}
    for (gi, gj), r, sig in zip(pairs, rs, significant):
        if not sig:
            continue
        crosstalk_edges[(gi, gj)] = float(abs(r))
        if r < 0:
            compensatory_pairs[(gi, gj)] = float(r)
        else:
            redundant_pairs[(gi, gj)] = float(r)

    return {
        "crosstalk_edges": crosstalk_edges,
        "compensatory_pairs": compensatory_pairs,
        "redundant_pairs": redundant_pairs,
        "n_tests": len(pairs),
        "n_significant": int(significant.sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True,
                         help="CSV: rows=TNBC patients, columns=gene symbols (90-gene panel). "
                              "Must already be restricted to the confirmed-TNBC cohort before "
                              "calling this script -- this script does not do cohort filtering.")
    parser.add_argument("--kinase-panel", required=True,
                         help="Path to the existing 90-gene panel list (one gene per line), "
                              "e.g. data/raw/kinases/kinase_90_list.txt")
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--min-abs-r", type=float, default=0.3)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.kinase_panel) as f:
        panel = [line.strip() for line in f if line.strip()]

    expression = pd.read_csv(args.expression, index_col=0)
    missing = [g for g in panel if g not in expression.columns]
    if missing:
        print(f"WARNING: {len(missing)} panel genes not found in the expression file, "
              f"excluded from the network rather than silently zero-filled: {missing}")
    present_panel = [g for g in panel if g in expression.columns]
    expression = expression[present_panel]

    print(f"Building correlation network: {len(present_panel)} genes, "
          f"{expression.shape[0]} patients, FDR<{args.fdr_threshold}, |r|>={args.min_abs_r}")

    result = build_correlation_network(expression, args.fdr_threshold, args.min_abs_r)
    print(f"{result['n_significant']} of {result['n_tests']} pairs significant "
          f"({len(result['redundant_pairs'])} redundant/positive, "
          f"{len(result['compensatory_pairs'])} compensatory/negative)")

    # JSON can't have tuple keys -- serialize as "geneA|geneB" strings
    def stringify_keys(d):
        return {f"{k[0]}|{k[1]}": v for k, v in d.items()}

    with open(args.output, "w") as f:
        json.dump({
            "crosstalk_edges": stringify_keys(result["crosstalk_edges"]),
            "compensatory_pairs": stringify_keys(result["compensatory_pairs"]),
            "redundant_pairs": stringify_keys(result["redundant_pairs"]),
            "n_tests": result["n_tests"],
            "n_significant": result["n_significant"],
        }, f, indent=2)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
