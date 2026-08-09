"""
run_real_essentiality_comparison.py

Orchestrates the real, already-tested modules found in this project's
codebase (depmap_multiomic_loader.py, depmap_supplemental_loader.py,
dependency_ml_model.py) against REAL DepMap files, using the REAL,
confirmed 34-ID TNBC cohort (tnbc_model_ids.txt) and the real 90-gene
panel -- to produce a genuinely reproducible version of the gene-identity
vs. per-sample multi-omic model comparison, since an exhaustive real
search across this project's entire codebase found NO locatable script
that actually computed the previously-cited 0.337/-0.029 numbers from
real data. This does not assume those numbers are wrong; it computes a
fresh, traceable result instead of continuing to cite an unverifiable one.

REAL DESIGN DECISION, stated explicitly rather than silently assumed:
CRISPRGeneDependency.csv (a 0-1 dependency PROBABILITY) is used as the
real training LABEL, matching dependency_ml_model.py's own docstring
("dependency_prob") and parse_depmap_dependency()'s explicit mapping --
NOT CRISPRGeneEffect.csv (Chronos score), even though Chronos is what
CTS's own essentiality component uses elsewhere in this project. This
means the R^2 computed here answers "how well is dependency PROBABILITY
predicted", a different (related) question than "how well is Chronos
SCORE predicted" -- stated here so this distinction is never silently
lost. If a Chronos-target version is wanted instead for closer
consistency with CTS itself, swap parse_depmap_dependency() for
parse_depmap_gene_effect() (Chronos loader, kinase_data_fetchers.py) with
value_name="dependency_prob" (the label column name train_and_evaluate()
expects) and rerun -- both are legitimate; this module defaults to the
one actually matching the pre-existing code's stated design.

REAL COHORT NOTE:
tnbc_model_ids.txt (34 real IDs) is used as the starting cohort. This
module reports how many of the 34 survive the real per-file
completeness/QC filtering built into the loaders (e.g., the "default
entry per model" filter on the mutation matrices, and dropna() on the
final dependency-label column in build_feature_matrix()) -- rather than
assuming any specific final count. If the real final n comes out to 25,
that would be the first real, traceable confirmation of that number this
whole project has produced; if it comes out to something else, that
number is what should be cited going forward instead.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd


def load_real_kinase_panel(path: str) -> List[str]:
    return [g.strip() for g in open(path) if g.strip()]


def load_real_tnbc_cohort(path: str) -> List[str]:
    return [line.strip() for line in open(path) if line.strip()]


def run_real_comparison(
    depmap_dir: str,
    kinase_panel_path: str,
    tnbc_cohort_path: str,
    loader_module_dir: str,
    label_source: str = "dependency_probability",
    common_essentials_path: Optional[str] = None,
) -> dict:
    """
    label_source: "dependency_probability" (default, uses CRISPRGeneDependency.csv,
    matching the pre-existing code's own stated design) or "chronos" (uses
    CRISPRGeneEffect.csv instead, for closer consistency with CTS's own
    essentiality component -- see module docstring for the real tradeoff).
    """
    sys.path.insert(0, loader_module_dir)
    from depmap_multiomic_loader import (
        parse_depmap_expression, parse_depmap_cnv, parse_depmap_dependency,
        parse_depmap_damaging_mutations, parse_depmap_wide_matrix,
    )
    from dependency_ml_model import build_feature_matrix, train_and_evaluate, compute_gene_identity_baseline

    kinase_panel = load_real_kinase_panel(kinase_panel_path)
    tnbc_ids = load_real_tnbc_cohort(tnbc_cohort_path)
    print(f"Starting real cohort: {len(tnbc_ids)} TNBC model IDs, {len(kinase_panel)} kinase-panel genes.")

    d = Path(depmap_dir)
    print("\nLoading real expression...")
    expr = parse_depmap_expression(str(d / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"), kinase_panel, tnbc_ids)
    print("Loading real copy number...")
    cnv = parse_depmap_cnv(str(d / "OmicsCNGeneWGS.csv"), kinase_panel, tnbc_ids)
    print("Loading real damaging mutations...")
    mut = parse_depmap_damaging_mutations(str(d / "OmicsSomaticMutationsMatrixDamaging.csv"), kinase_panel, tnbc_ids)

    print(f"\nLoading real dependency label (source={label_source})...")
    if label_source == "dependency_probability":
        dep = parse_depmap_dependency(str(d / "CRISPRGeneDependency.csv"), kinase_panel, tnbc_ids)
    elif label_source == "chronos":
        dep = parse_depmap_wide_matrix(str(d / "CRISPRGeneEffect.csv"), kinase_panel, tnbc_ids, value_name="dependency_prob")
    else:
        raise ValueError(f"Unknown label_source {label_source!r}; use 'dependency_probability' or 'chronos'.")

    pan_essential_df = None
    if common_essentials_path:
        from depmap_supplemental_loader import flag_pan_essential_kinases
        pan_essential_df = flag_pan_essential_kinases(kinase_panel, common_essentials_path)

    real_ids_in_expr = set(expr["cell_line_id"])
    real_ids_in_dep = set(dep["cell_line_id"])
    print(f"\nReal cell lines with expression data: {len(real_ids_in_expr)}/{len(tnbc_ids)}")
    print(f"Real cell lines with dependency-label data: {len(real_ids_in_dep)}/{len(tnbc_ids)}")

    feature_df = build_feature_matrix(expr, cnv, mut, dep, pan_essential_df=pan_essential_df)
    final_n_lines = feature_df["cell_line_id"].nunique()
    print(f"\nFINAL real cohort actually used for model fitting: {final_n_lines} cell lines "
          f"(started from {len(tnbc_ids)} in the real TNBC list).")
    if final_n_lines != len(tnbc_ids):
        dropped_ids = set(tnbc_ids) - set(feature_df["cell_line_id"].unique())
        print(f"  {len(tnbc_ids) - final_n_lines} real ID(s) dropped due to missing feature/label data: "
              f"{sorted(dropped_ids)}")

    print("\nRunning real gene-identity baseline (grouped CV)...")
    baseline = compute_gene_identity_baseline(feature_df)
    print(f"  Gene-identity R^2 = {baseline['r2']:.4f} (Spearman rho = {baseline['spearman_r']:.4f})")

    print("\nRunning real per-sample multi-omic model (grouped CV)...")
    multiomic = train_and_evaluate(feature_df)
    print(f"  Multi-omic model R^2 = {multiomic['r2']:.4f} (Spearman rho = {multiomic['spearman_r']:.4f})")

    return {
        "label_source": label_source,
        "n_starting_tnbc_ids": len(tnbc_ids),
        "n_final_lines_used": final_n_lines,
        "n_kinase_genes": len(kinase_panel),
        "gene_identity_r2": baseline["r2"],
        "multiomic_model_r2": multiomic["r2"],
        "gene_identity_spearman": baseline["spearman_r"],
        "multiomic_spearman": multiomic["spearman_r"],
    }


# =====================================================================
# SMOKE TEST -- synthetic files matching every real confirmed schema,
# verifies the full pipeline runs end-to-end and produces a real,
# non-trivial gene-identity-vs-multiomic comparison before touching
# actual DepMap downloads.
# =====================================================================

def _run_smoke_test():
    import numpy as np
    import tempfile, os

    tmpdir = tempfile.mkdtemp()
    rng = np.random.default_rng(7)
    kinases = ["EGFR", "ERBB2", "PTK2", "FGFR1", "SRC"]
    cell_lines = [f"ACH-{i:06d}" for i in range(12)]

    def wide_file(path, value_range, seed, gene_dependent=False):
        r = np.random.default_rng(seed)
        if gene_dependent:
            # plant a real gene-identity effect: each gene has its own fixed
            # baseline level, mostly independent of which cell line it's in --
            # this is what compute_gene_identity_baseline() should pick up on.
            gene_base = {g: r.uniform(*value_range) for g in kinases}
            data = np.array([[gene_base[g] + r.normal(0, 0.03) for g in kinases] for _ in cell_lines])
        else:
            data = r.uniform(*value_range, size=(len(cell_lines), len(kinases)))
        cols = [f"{g} ({1000+i})" for i, g in enumerate(kinases)]
        pd.DataFrame(data, index=cell_lines, columns=cols).rename_axis("ModelID").to_csv(path)

    wide_file(f"{tmpdir}/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", (0, 10), 1)
    wide_file(f"{tmpdir}/OmicsCNGeneWGS.csv", (-1, 2), 2)
    wide_file(f"{tmpdir}/CRISPRGeneDependency.csv", (0, 1), 3, gene_dependent=True)

    mut_rows = []
    for cl in cell_lines:
        row = {"ModelID": cl, "SequencingID": "S1", "ModelConditionID": "MC1", "IsDefaultEntryForModel": "Yes"}
        for g in kinases:
            row[f"{g} (1000)"] = int(rng.integers(0, 2))
        mut_rows.append(row)
    pd.DataFrame(mut_rows).to_csv(f"{tmpdir}/OmicsSomaticMutationsMatrixDamaging.csv", index=False)

    with open(f"{tmpdir}/kinase_panel.txt", "w") as f:
        f.write("\n".join(kinases))
    with open(f"{tmpdir}/tnbc_ids.txt", "w") as f:
        f.write("\n".join(cell_lines))

    loader_dir = str(Path(__file__).parent)
    result = run_real_comparison(tmpdir, f"{tmpdir}/kinase_panel.txt", f"{tmpdir}/tnbc_ids.txt", loader_dir)

    print("\n" + "=" * 70)
    print("SMOKE TEST RESULT:", result)
    assert result["n_final_lines_used"] == len(cell_lines), "no lines should have been dropped in this clean synthetic case"
    assert result["gene_identity_r2"] > 0.3, "planted gene-identity signal should be recovered"
    print("PASSED: full real orchestration runs end-to-end on synthetic data matching every real confirmed schema.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        _run_smoke_test()
    else:
        ap = argparse.ArgumentParser()
        ap.add_argument("depmap_dir", help="Directory containing the real CRISPRGeneDependency.csv, "
                         "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv, OmicsCNGeneWGS.csv, "
                         "OmicsSomaticMutationsMatrixDamaging.csv")
        ap.add_argument("kinase_panel_path")
        ap.add_argument("tnbc_cohort_path", help="Real tnbc_model_ids.txt (34 real IDs)")
        ap.add_argument("loader_module_dir", help="Directory containing depmap_multiomic_loader.py, "
                         "depmap_supplemental_loader.py, dependency_ml_model.py")
        ap.add_argument("--label-source", choices=["dependency_probability", "chronos"], default="dependency_probability")
        ap.add_argument("--common-essentials-path", default=None)
        args = ap.parse_args()

        result = run_real_comparison(
            args.depmap_dir, args.kinase_panel_path, args.tnbc_cohort_path, args.loader_module_dir,
            label_source=args.label_source, common_essentials_path=args.common_essentials_path,
        )
        print("\n" + "=" * 70)
        print("FINAL REAL RESULT:", result)
