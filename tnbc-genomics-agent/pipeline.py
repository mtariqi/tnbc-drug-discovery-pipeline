"""
pipeline.py
End-to-end TNBC RTK/nRTK redundancy analysis pipeline orchestrator.
Run: python pipeline.py --vcf data/sample/patient_1.vcf
"""

import argparse
import json
import os
import sys
from datetime import datetime

# Allow running from project root
sys.path.insert(0, os.path.dirname(__file__))

from modules.vcf_parser import VCFParser
from modules.redundancy_analyzer import RedundancyAnalyzer
from modules.kinase_db import KINASE_DATABASE


def print_section(title: str, char: str = "═", width: int = 70):
    print(f"\n{char * width}")
    print(f"  {title}")
    print(f"{char * width}")


def run_pipeline(vcf_path: str, output_json: str = None, verbose: bool = True):
    """
    Full pipeline:
      1. Parse VCF
      2. Extract RTK/nRTK variants
      3. Run redundancy analysis
      4. Print structured report
      5. Optionally dump JSON for downstream use (e.g., dashboard)
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if verbose:
        print_section(f"TNBC RTK/nRTK Redundancy Pipeline  |  {timestamp}")
        print(f"  Input VCF : {vcf_path}")

    # ── Step 1: Parse VCF ────────────────────────────────────────────────────
    if verbose:
        print_section("Step 1 — VCF Parsing", char="─")

    parser = VCFParser(vcf_path).parse()
    summary = parser.get_summary()

    if verbose:
        print(f"  Total variants parsed   : {summary['total_variants']}")
        print(f"  RTK/nRTK variants found : {summary['rtk_nrtk_variants']}")
        print(f"  Kinase genes detected   : {', '.join(sorted(summary['genes_detected']))}")

    # ── Step 2: Variant Detail ────────────────────────────────────────────────
    if verbose:
        print_section("Step 2 — RTK/nRTK Variant Details", char="─")
        for v in parser.rtk_nrtk_variants:
            info = KINASE_DATABASE.get(v.gene, {})
            print(
                f"  [{info.get('type','?'):4s}] {v.gene:<8}  "
                f"chr{v.chrom}:{v.pos}  "
                f"AF={v.allele_freq:.2f}  "
                f"Effect={v.effect:<25}  "
                f"Severity={v.severity_score:.2f}  "
                f"TNBC={info.get('tnbc_relevance','?')}"
            )

    # ── Step 3: Redundancy Analysis ──────────────────────────────────────────
    if verbose:
        print_section("Step 3 — Pathway Redundancy Analysis", char="─")

    analyzer = RedundancyAnalyzer(parser.rtk_nrtk_variants)
    report = analyzer.get_full_report()

    if verbose:
        print("\n  Pathway Redundancy Scores (higher = harder to block):")
        for pathway, score in sorted(
            report["pathway_redundancy_scores"].items(), key=lambda x: x[1], reverse=True
        ):
            bar = "█" * int(score * 20)
            print(f"  {pathway:<25}  {bar:<20}  {score:.3f}")

    # ── Step 4: Bypass Risk per Target ───────────────────────────────────────
    if verbose:
        print_section("Step 4 — Bypass Risk per Target", char="─")
        for item in report["bypass_risk_per_target"]:
            risk = item["bypass_risk_score"]
            flag = "🔴" if risk >= 0.6 else "🟡" if risk >= 0.3 else "🟢"
            print(f"\n  {flag}  {item['gene']} ({item['type']}) — Bypass Risk: {risk:.3f}")
            print(f"     {item['recommendation']}")
            if item["bypass_genes"]:
                bg = ", ".join(f"{b['gene']}({'+'.join(b['shared_pathways'])})" for b in item["bypass_genes"])
                print(f"     Potential compensators: {bg}")

    # ── Step 5: Combination Therapy Suggestions ──────────────────────────────
    if verbose:
        print_section("Step 5 — Combination Therapy Suggestions", char="─")
        for sug in report["combination_therapy_suggestions"][:5]:
            print(
                f"\n  ➤ {sug['target_gene_1']} + {sug['target_gene_2']}  "
                f"[{sug['shared_pathway']}]  "
                f"Evidence: {sug['evidence_strength'].upper()}"
            )
            print(f"     {sug['rationale']}")
            d1 = " / ".join(sug["suggested_drugs_gene1"])
            d2 = " / ".join(sug["suggested_drugs_gene2"])
            print(f"     Drugs: {d1}  +  {d2}")

    # ── Step 6: High-Severity Variants ───────────────────────────────────────
    if verbose and summary["high_severity_variants"]:
        print_section("Step 6 — High-Severity Variants (score ≥ 0.7)", char="─")
        for v in summary["high_severity_variants"]:
            print(
                f"  {v['gene']:<8}  {v['effect']:<30}  "
                f"AF={v['allele_freq']:.2f}  Severity={v['severity_score']:.2f}"
            )

    # ── Final JSON output ────────────────────────────────────────────────────
    full_output = {
        "pipeline_version": "1.0.0",
        "timestamp": timestamp,
        "input_vcf": vcf_path,
        "vcf_summary": summary,
        "redundancy_report": report,
    }

    if output_json:
        with open(output_json, "w") as fh:
            json.dump(full_output, fh, indent=2)
        if verbose:
            print(f"\n  ✅  JSON report saved → {output_json}")

    if verbose:
        print_section("Pipeline Complete", char="═")

    return full_output


# ── CLI entry point ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="TNBC RTK/nRTK Redundancy Pipeline")
    ap.add_argument("--vcf",  required=True,  help="Path to input VCF file")
    ap.add_argument("--out",  default=None,   help="Optional path to save JSON output")
    ap.add_argument("--quiet", action="store_true", help="Suppress console output")
    args = ap.parse_args()

    run_pipeline(
        vcf_path=args.vcf,
        output_json=args.out,
        verbose=not args.quiet,
    )
