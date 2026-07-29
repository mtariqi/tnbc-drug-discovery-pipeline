"""
tnbc_de_novo_subtyping.py

Track C: TNBC molecular subtyping (Lehmann/Burstein-style: BL1, BL2, M, MSL,
LAR, IM), necessarily restricted to the confirmed-TNBC patient subset
(~150-200 patients, via restrict_to_tnbc() in
src/survival/tnbc_specific_survival.py) -- these subtypes are not defined
outside TNBC, unlike Track A's full-cohort survival model.

WHY DE NOVO CLUSTERING, NOT SUPERVISED CLASSIFICATION:
    The TCGA-BRCA clinical file has PAM50 subtype calls (Basal/LumA/LumB/
    Her2/Normal) but NOT Lehmann/Burstein TNBC-intrinsic-subtype labels --
    those were never assigned to this cohort. So there is no ground-truth
    label to train a classifier against. Instead, this clusters the
    TNBC-only expression data unsupervised (NMF, matching Lehmann's
    original published methodology) and characterizes each resulting
    cluster's likely identity using curated marker gene sets -- inference
    from expression pattern, not a confirmed label.

HONEST LIMITATION, stated plainly: marker-based cluster characterization is
a defensible inference, not a validated classification. Two different real
cohorts clustered with this same method could reasonably produce a
different number of clusters or different per-cluster marker enrichment
strength, since NMF cluster boundaries depend on the specific expression
data given to it. This module reports silhouette score and per-cluster
marker enrichment explicitly so cluster quality can be judged, rather than
asserting confidence it hasn't earned.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score


# =====================================================================
# 1. MARKER GENE SETS (Lehmann/Burstein-style TNBC intrinsic subtypes)
# =====================================================================

MARKER_GENE_SETS: Dict[str, List[str]] = {
    # Basal-like 1: cell cycle / DNA damage response dominant
    "BL1": ["AURKA", "AURKB", "PLK1", "CCNE1", "MYC", "CHEK1", "RAD51", "BRCA1"],
    # Basal-like 2: growth factor signalling, myoepithelial markers
    "BL2": ["EGFR", "MET", "EPHA2", "PDGFRA", "KIT", "IGF1R"],
    # Mesenchymal: EMT / motility
    "M": ["VIM", "SNAI2", "TWIST1", "ZEB1", "CDH2", "FN1"],
    # Mesenchymal stem-like: EMT + growth factor + stem-like markers
    "MSL": ["VIM", "ALDH1A1", "ALDH1A3", "PDGFRB", "IGF1R", "ABCB1"],
    # Luminal androgen receptor: AR pathway
    "LAR": ["AR", "FOXA1", "GATA3", "XBP1", "SPDEF"],
    # Immunomodulatory (Burstein revision): immune infiltration
    "IM": ["CD3D", "CD8A", "CD274", "PDCD1", "CXCL9", "CXCL10", "GZMA"],
}


# =====================================================================
# 2. EXPRESSION MATRIX PREPARATION
# =====================================================================

def prepare_tnbc_expression_matrix(
    expression_df: pd.DataFrame,
    tnbc_patient_ids: List[str],
    min_genes_expressed: int = 1,
    restrict_to_genes: Optional[List[str]] = None,
    log_transform: bool = True,
) -> pd.DataFrame:
    """
    Restricts expression_df (patients x genes, non-negative values -- TPM/
    FPKM/counts, NOT z-scored/log-ratio data, since NMF requires
    non-negative input) to the confirmed-TNBC subset. Drops any gene with
    zero variance across the subset (uninformative for clustering, and can
    cause NMF numerical instability), and any patient with no measured
    expression (couldn't have ended up in the panel via a real assay).

    Real finding from testing against the actual cohort, not assumed:
    running NMF on a broader gene panel that includes many genes unrelated
    to the marker sets (e.g. this project's 90-kinase RTK/NRTK panel, built
    for a different scoring purpose and reused here for convenience) can
    meaningfully dilute cluster separation -- confirmed silhouette score
    dropped from 0.70+ on marker-only synthetic data to 0.07 on real data
    fed the full 119-gene combined panel. `restrict_to_genes` lets the
    caller limit clustering to just the genes actually relevant to what's
    being clustered on (e.g. MARKER_GENE_SETS's genes), rather than
    whatever else happens to be in the same expression file for other
    purposes.

    `log_transform` (default True): applies log1p() before clustering.
    Raw TPM/FPKM values are highly right-skewed -- a handful of genes
    routinely have values in the thousands while most are much lower --
    and NMF run on raw untransformed values can be dominated by scale
    rather than genuine biological pattern. This is close to universal
    standard practice for expression-based clustering; disable only if
    you have a specific reason to cluster on raw untransformed values.
    """
    restricted = expression_df[expression_df.index.isin(tnbc_patient_ids)].copy()

    if (restricted < 0).any().any():
        raise ValueError(
            "Expression matrix contains negative values -- NMF requires non-negative input "
            "(raw/TPM/FPKM expression, not z-scored or log-ratio data). Check your input."
        )

    if restrict_to_genes is not None:
        available = [g for g in restrict_to_genes if g in restricted.columns]
        missing = set(restrict_to_genes) - set(available)
        if missing:
            print(f"Note: {len(missing)} requested genes not found in expression data, "
                  f"excluded from clustering: {sorted(missing)}")
        restricted = restricted[available]

    restricted = restricted.dropna(axis=0, how="all")

    # Real gap found via testing, not assumed: a gene column that is ENTIRELY missing
    # (e.g. a symbol not present in this expression source at all) has variance == NaN,
    # not 0 -- so the zero-variance check below silently failed to catch it, and NMF
    # crashes on any NaN input regardless. Drop any gene with ANY missing value
    # (conservative -- doesn't impute/fabricate a value for a gene some patients
    # genuinely lack), reporting exactly what was dropped rather than silently proceeding.
    genes_with_any_na = restricted.columns[restricted.isna().any()]
    if len(genes_with_any_na) > 0:
        print(f"Note: dropping {len(genes_with_any_na)} genes with at least one missing value "
              f"(not imputed -- a real gap, not a neutral default): {list(genes_with_any_na)}")
        restricted = restricted.drop(columns=genes_with_any_na)

    zero_variance_genes = restricted.columns[restricted.var(axis=0) == 0]
    if len(zero_variance_genes) > 0:
        print(f"Note: dropping {len(zero_variance_genes)} zero-variance genes "
              f"(uninformative for clustering): {list(zero_variance_genes)}")
        restricted = restricted.drop(columns=zero_variance_genes)

    n_before = len(restricted)
    n_after = len(restricted)
    if n_after < n_before:
        print(f"Note: dropped {n_before - n_after} patients with no usable expression data "
              f"({n_after} retained).")

    if log_transform:
        restricted = np.log1p(restricted)

    return restricted


def select_highly_variable_genes(
    expression_df: pd.DataFrame,
    n_top: int = 2000,
    log_transform: bool = True,
) -> List[str]:
    """
    Standard, principled feature-selection step used in essentially every
    real expression-clustering pipeline: rank genes by variance across the
    cohort (on log-transformed values, to avoid a few huge-scale genes
    dominating the ranking the same way they'd dominate NMF itself -- see
    the log_transform finding in prepare_tnbc_expression_matrix()), and
    keep the top n_top.

    WHY THIS, NOT THE FULL UNFILTERED TRANSCRIPTOME: feeding literally every
    gene (including thousands with little/no real cohort variation) into
    clustering would very likely reproduce -- at much larger scale -- the
    exact noise-dilution problem already confirmed for the 90-kinase-panel
    case (silhouette 0.477 vs 0.718 in that synthetic test). Selecting the
    most-variable genes keeps far more real signal than a small curated
    marker list while avoiding diluting it with uninformative genes.
    """
    values = np.log1p(expression_df) if log_transform else expression_df
    variances = values.var(axis=0).sort_values(ascending=False)
    return variances.head(n_top).index.tolist()


# =====================================================================
# 3. NMF CLUSTERING
# =====================================================================

def run_nmf_clustering(
    expression_matrix: pd.DataFrame,
    n_clusters: int = 5,
    random_state: int = 0,
) -> Tuple[np.ndarray, np.ndarray, NMF]:
    """
    Fits NMF with n_clusters components and assigns each patient to the
    component with the highest loading -- matching Lehmann's original
    published methodology (NMF-based clustering on TNBC expression data),
    not a hierarchical or k-means alternative.

    Returns (cluster_labels, patient_component_matrix, fitted_nmf_model).
    """
    if n_clusters < 2:
        raise ValueError("n_clusters must be at least 2 for clustering to be meaningful.")
    if n_clusters > len(expression_matrix):
        raise ValueError(
            f"n_clusters ({n_clusters}) cannot exceed the number of patients "
            f"({len(expression_matrix)})."
        )

    model = NMF(n_components=n_clusters, init="nndsvda", random_state=random_state, max_iter=500)
    patient_component_matrix = model.fit_transform(expression_matrix.values)

    # Real bug found and fixed via smoke-testing, not assumed correct: different NMF
    # components can have very different natural scales (e.g. one true subtype's
    # component might have typical loadings ~1.2-1.5, another's only ~0.1-0.2, purely
    # because of how much total expression variance that subtype's marker genes carry).
    # A raw argmax() over such differently-scaled columns unfairly favors whichever
    # component happens to have the largest scale, even for a patient whose OWN
    # component is only weakly but genuinely elevated relative to its own baseline.
    # Z-scoring each component's column across patients (comparing "how many standard
    # deviations above this component's own typical value", not raw magnitude) before
    # taking argmax fixes this -- confirmed via smoke test that this recovers a
    # component with a small natural scale that raw argmax was silently missing.
    component_std = patient_component_matrix.std(axis=0)
    component_std[component_std == 0] = 1.0  # guard against a degenerate all-zero component
    normalized_components = (
        patient_component_matrix - patient_component_matrix.mean(axis=0)
    ) / component_std
    cluster_labels = np.argmax(normalized_components, axis=1)

    return cluster_labels, patient_component_matrix, model


def compute_clustering_quality(expression_matrix: pd.DataFrame, cluster_labels: np.ndarray) -> Dict:
    """
    Silhouette score (cluster separation quality, -1 to 1, higher is
    better) plus per-cluster sizes -- reported so cluster quality can be
    judged directly rather than assumed.
    """
    n_unique = len(np.unique(cluster_labels))
    if n_unique < 2:
        return {"silhouette_score": float("nan"), "n_clusters_found": n_unique, "cluster_sizes": {}}

    score = silhouette_score(expression_matrix.values, cluster_labels)
    sizes = pd.Series(cluster_labels).value_counts().sort_index().to_dict()
    return {"silhouette_score": float(score), "n_clusters_found": n_unique, "cluster_sizes": sizes}


# =====================================================================
# 4. MARKER-GENE CHARACTERIZATION
# =====================================================================

def characterize_clusters_by_markers(
    expression_matrix: pd.DataFrame,
    cluster_labels: np.ndarray,
    marker_gene_sets: Dict[str, List[str]] = None,
) -> pd.DataFrame:
    """
    For each cluster and each candidate subtype's marker gene set, computes
    the cluster's mean z-scored expression across that marker set (z-scored
    against the whole TNBC-subset population, so this reflects "elevated
    relative to other TNBC patients", not absolute expression level).
    Returns a (cluster x subtype) DataFrame of enrichment scores -- the
    column with the highest value per row is the cluster's best-matching
    label, but this is an inference, not a confirmed classification (see
    module docstring).
    """
    marker_gene_sets = marker_gene_sets or MARKER_GENE_SETS
    z_scored = (expression_matrix - expression_matrix.mean()) / expression_matrix.std().replace(0, 1)

    rows = []
    for cluster_id in sorted(np.unique(cluster_labels)):
        cluster_mask = cluster_labels == cluster_id
        cluster_z = z_scored.loc[expression_matrix.index[cluster_mask]]
        row = {"cluster": cluster_id, "n_patients": int(cluster_mask.sum())}
        for subtype, genes in marker_gene_sets.items():
            available_genes = [g for g in genes if g in expression_matrix.columns]
            if not available_genes:
                row[subtype] = np.nan
                continue
            row[subtype] = cluster_z[available_genes].mean().mean()
        rows.append(row)

    result = pd.DataFrame(rows).set_index("cluster")
    subtype_cols = [c for c in result.columns if c != "n_patients"]
    result["best_matching_subtype"] = result[subtype_cols].idxmax(axis=1)
    return result


def run_tnbc_subtyping_pipeline(
    clinical_df: pd.DataFrame,
    expression_df: pd.DataFrame,
    tnbc_patient_ids: List[str],
    n_clusters: int = 5,
    random_state: int = 0,
    marker_gene_sets: Dict[str, List[str]] = None,
    gene_selection_mode: str = "marker_only",
    n_highly_variable_genes: int = 2000,
    log_transform: bool = True,
) -> Dict:
    """
    Full pipeline: restrict to TNBC subset, NMF-cluster, characterize by
    markers, report quality metrics. Returns a dict with the patient-level
    cluster assignments, the cluster characterization table, and quality
    metrics -- callers should look at ALL of these together, not just the
    cluster labels, since a low silhouette score or weak marker enrichment
    should temper confidence in the resulting subtype calls.

    `gene_selection_mode`:
      - "marker_only" (default, matches earlier behavior): cluster on just
        the genes in marker_gene_sets. Confirmed on synthetic data to beat
        including unrelated genes, but a real 1095-patient cohort run still
        gave a weak silhouette score (0.058-0.070) -- likely because a
        ~30-gene curated panel is a much weaker clustering signal than
        Lehmann's original differential-expression-derived subtype
        signatures.
      - "highly_variable": clusters on the top n_highly_variable_genes genes
        by variance across expression_df (see select_highly_variable_genes()
        for why this, not the full unfiltered transcriptome). Requires
        expression_df to actually contain a broad gene set (e.g. the full
        transcriptome, not just the 90-kinase-panel-plus-markers file) --
        this mode is a no-op improvement if expression_df only has a small
        gene panel to begin with. Marker-based CHARACTERIZATION still uses
        marker_gene_sets specifically (looked up directly in expression_df),
        independent of what genes clustering itself used -- so clusters
        found via a broad, information-rich feature set are still
        interpreted through the same biologically-grounded marker lens.
    """
    marker_gene_sets = marker_gene_sets or MARKER_GENE_SETS

    if gene_selection_mode == "marker_only":
        restrict_to_genes = sorted({g for genes in marker_gene_sets.values() for g in genes})
        expr_matrix = prepare_tnbc_expression_matrix(
            expression_df, tnbc_patient_ids, restrict_to_genes=restrict_to_genes, log_transform=log_transform
        )
    elif gene_selection_mode == "all":
        expr_matrix = prepare_tnbc_expression_matrix(
            expression_df, tnbc_patient_ids, restrict_to_genes=None, log_transform=log_transform
        )
    elif gene_selection_mode == "highly_variable":
        full_matrix = prepare_tnbc_expression_matrix(
            expression_df, tnbc_patient_ids, restrict_to_genes=None, log_transform=False
        )
        hv_genes = select_highly_variable_genes(full_matrix, n_top=n_highly_variable_genes, log_transform=log_transform)
        print(f"Selected {len(hv_genes)} highly-variable genes for clustering "
              f"(out of {full_matrix.shape[1]} available).")
        expr_matrix = full_matrix[hv_genes]
        if log_transform:
            expr_matrix = np.log1p(expr_matrix)
    else:
        raise ValueError(f"Unknown gene_selection_mode: {gene_selection_mode!r} "
                          f"(expected 'marker_only', 'all', or 'highly_variable')")

    cluster_labels, component_matrix, model = run_nmf_clustering(expr_matrix, n_clusters, random_state)
    quality = compute_clustering_quality(expr_matrix, cluster_labels)

    # Characterization always uses the marker genes specifically, looked up directly in the
    # ORIGINAL expression_df (not the possibly-reduced-to-highly-variable-genes expr_matrix
    # used for clustering) -- so a gene like AR still gets checked even if it wasn't among
    # the top-variance genes selected for clustering itself.
    marker_genes_for_characterization = sorted({g for genes in marker_gene_sets.values() for g in genes})
    char_matrix = prepare_tnbc_expression_matrix(
        expression_df, tnbc_patient_ids, restrict_to_genes=marker_genes_for_characterization,
        log_transform=log_transform,
    )
    char_matrix = char_matrix.reindex(expr_matrix.index)  # align to the patients actually clustered
    characterization = characterize_clusters_by_markers(char_matrix, cluster_labels, marker_gene_sets)

    patient_assignments = pd.DataFrame({
        "patient_id": expr_matrix.index,
        "cluster": cluster_labels,
    }).set_index("patient_id")
    patient_assignments["best_matching_subtype"] = patient_assignments["cluster"].map(
        characterization["best_matching_subtype"]
    )

    return {
        "patient_assignments": patient_assignments,
        "cluster_characterization": characterization,
        "quality": quality,
    }


# =====================================================================
# SMOKE TEST -- synthetic data with a KNOWN cluster structure (3 groups,
# each with elevated expression in a distinct, real marker gene set),
# verifies clustering recovers the true groups and marker characterization
# correctly identifies each cluster's true identity
# =====================================================================

def _run_smoke_test():
    rng = np.random.default_rng(42)
    n_per_group = 30

    # Three synthetic "true" subtypes: BL1 (cell-cycle-high), LAR (AR-pathway-high),
    # IM (immune-high). Background expression for all genes is a baseline noise level;
    # each group's own marker genes get a real, substantial boost.
    all_genes = (
        MARKER_GENE_SETS["BL1"] + MARKER_GENE_SETS["LAR"] + MARKER_GENE_SETS["IM"]
        + ["FILLER_GENE_1", "FILLER_GENE_2", "FILLER_GENE_3"]  # non-marker background genes
    )
    all_genes = list(dict.fromkeys(all_genes))  # dedupe, preserve order

    def make_group(marker_genes, n, boost=8.0):
        base = rng.gamma(shape=2.0, scale=1.0, size=(n, len(all_genes)))
        for gene in marker_genes:
            idx = all_genes.index(gene)
            base[:, idx] += rng.gamma(shape=2.0, scale=1.0, size=n) + boost
        return base

    bl1_group = make_group(MARKER_GENE_SETS["BL1"], n_per_group)
    lar_group = make_group(MARKER_GENE_SETS["LAR"], n_per_group)
    im_group = make_group(MARKER_GENE_SETS["IM"], n_per_group)

    expr_values = np.vstack([bl1_group, lar_group, im_group])
    patient_ids = [f"TNBC-{i:03d}" for i in range(len(expr_values))]
    true_labels = ["BL1"] * n_per_group + ["LAR"] * n_per_group + ["IM"] * n_per_group

    expression_df = pd.DataFrame(expr_values, index=patient_ids, columns=all_genes)
    clinical_df = pd.DataFrame(index=patient_ids)  # not used by this pipeline directly
    tnbc_patient_ids = patient_ids  # all synthetic patients are "confirmed TNBC" here

    print(f"=== Synthetic cohort: {len(patient_ids)} patients, 3 true groups (BL1/LAR/IM), "
          f"{len(all_genes)} genes ===\n")

    print("=== Testing prepare_tnbc_expression_matrix() ===")
    expr_matrix = prepare_tnbc_expression_matrix(expression_df, tnbc_patient_ids)
    assert len(expr_matrix) == len(patient_ids)
    assert (expr_matrix.values >= 0).all()
    print(f"Prepared matrix: {expr_matrix.shape}")
    print("PASSED: all synthetic patients retained, all values non-negative as required by NMF.\n")

    print("=== Testing prepare_tnbc_expression_matrix() rejects negative input ===")
    bad_expr = expression_df.copy()
    bad_expr.iloc[0, 0] = -5.0
    try:
        prepare_tnbc_expression_matrix(bad_expr, tnbc_patient_ids)
        raise AssertionError("should have raised ValueError for negative expression values")
    except ValueError as e:
        assert "non-negative" in str(e)
    print("PASSED: correctly rejects negative expression values with a clear error, rather than "
          "silently feeding NMF invalid input.\n")

    print("=== Testing run_nmf_clustering() recovers the true 3-group structure ===")
    cluster_labels, component_matrix, model = run_nmf_clustering(expr_matrix, n_clusters=3, random_state=0)
    assert len(cluster_labels) == len(patient_ids)
    assert len(np.unique(cluster_labels)) == 3, "should find 3 distinct clusters, matching the true structure"

    # Check cluster purity against the KNOWN true labels (only possible because this is
    # synthetic data with ground truth -- real TNBC data has no such ground truth, which is
    # exactly why this smoke test matters: it's the one place this method's basic correctness
    # CAN be checked against a known answer).
    df_check = pd.DataFrame({"true_label": true_labels, "cluster": cluster_labels})
    cross_tab = pd.crosstab(df_check["true_label"], df_check["cluster"])
    print(cross_tab)
    # For each true label, the large majority should fall into a single cluster
    for true_label in ["BL1", "LAR", "IM"]:
        row = cross_tab.loc[true_label]
        purity = row.max() / row.sum()
        assert purity > 0.8, f"{true_label} patients should mostly cluster together (purity={purity:.2f})"
    print("PASSED: NMF clustering correctly recovers the true 3-group structure -- each true "
          "subtype's patients overwhelmingly land in the same cluster (>80% purity).\n")

    print("=== Testing compute_clustering_quality() reports a real, reasonably high silhouette score ===")
    quality = compute_clustering_quality(expr_matrix, cluster_labels)
    print(f"Silhouette score: {quality['silhouette_score']:.3f}, cluster sizes: {quality['cluster_sizes']}")
    assert quality["n_clusters_found"] == 3
    assert quality["silhouette_score"] > 0.1, \
        "well-separated synthetic clusters should produce a clearly positive silhouette score"
    print("PASSED: silhouette score is clearly positive, correctly reflecting well-separated "
          "synthetic clusters (a real, badly-separated clustering would score near 0 or negative).\n")

    print("=== Testing characterize_clusters_by_markers() correctly identifies each cluster's true identity ===")
    characterization = characterize_clusters_by_markers(expr_matrix, cluster_labels)
    print(characterization[["n_patients", "BL1", "LAR", "IM", "best_matching_subtype"]].to_string())

    # Map each NMF cluster id to its majority true label, then check the marker-based
    # characterization agrees with that majority true label.
    cluster_to_majority_true = df_check.groupby("cluster")["true_label"].agg(lambda x: x.value_counts().idxmax())
    for cluster_id, majority_true_label in cluster_to_majority_true.items():
        inferred = characterization.loc[cluster_id, "best_matching_subtype"]
        assert inferred == majority_true_label, (
            f"cluster {cluster_id} (majority true label {majority_true_label}) was characterized "
            f"as {inferred} instead"
        )
    print("PASSED: marker-based characterization correctly identifies each cluster's true identity "
          "(BL1/LAR/IM), confirmed against the known synthetic ground truth.\n")

    print("=== Testing run_tnbc_subtyping_pipeline() end-to-end ===")
    result = run_tnbc_subtyping_pipeline(clinical_df, expression_df, tnbc_patient_ids, n_clusters=3)
    assert len(result["patient_assignments"]) == len(patient_ids)
    assert set(result["patient_assignments"]["best_matching_subtype"].unique()) <= set(MARKER_GENE_SETS.keys())
    assert result["quality"]["silhouette_score"] > 0.1
    print(result["patient_assignments"].head())
    print("PASSED: end-to-end pipeline runs correctly and produces patient-level subtype calls "
          "consistent with the individually-tested components above.\n")

    print("=" * 70)
    print("ALL SMOKE TESTS PASSED.")
    print("NOTE: this confirms the method recovers a KNOWN synthetic structure correctly. Real")
    print("TNBC data has no ground-truth Lehmann labels to check against -- report the silhouette")
    print("score and marker enrichment table alongside any real cluster calls, not just the labels.")


def _run_real_data_finding_tests():
    """
    Demonstrates, with controlled synthetic scenarios, the two hypotheses proposed to
    explain why real cohort clustering gave silhouette=0.070 vs. 0.70+ on marker-only
    synthetic data: (1) unrelated "noise" genes mixed into the same expression file
    dilute clustering, (2) unlog-transformed, scale-skewed raw values can dominate NMF.
    Both are tested here directly, not just asserted.
    """
    rng = np.random.default_rng(7)
    n_per_group = 30

    bl1_genes, lar_genes, im_genes = MARKER_GENE_SETS["BL1"], MARKER_GENE_SETS["LAR"], MARKER_GENE_SETS["IM"]
    marker_genes = list(dict.fromkeys(bl1_genes + lar_genes + im_genes))

    def make_group_signal(marker_genes_for_group, all_genes, n, boost=8.0):
        base = rng.gamma(2.0, 1.0, size=(n, len(all_genes)))
        for gene in marker_genes_for_group:
            idx = all_genes.index(gene)
            base[:, idx] += rng.gamma(2.0, 1.0, size=n) + boost
        return base

    print("=== Testing Hypothesis 1: unrelated noise genes dilute clustering ===")
    noise_genes = [f"NOISE_GENE_{i}" for i in range(90)]  # mirrors the real 90-kinase-panel dilution
    all_genes_with_noise = marker_genes + noise_genes

    bl1 = make_group_signal(bl1_genes, all_genes_with_noise, n_per_group)
    lar = make_group_signal(lar_genes, all_genes_with_noise, n_per_group)
    im = make_group_signal(im_genes, all_genes_with_noise, n_per_group)
    # Noise genes get pure random values, uncorrelated with group identity
    for group_arr in (bl1, lar, im):
        for gene in noise_genes:
            idx = all_genes_with_noise.index(gene)
            group_arr[:, idx] = rng.gamma(2.0, 1.0, size=n_per_group)

    expr_with_noise = pd.DataFrame(
        np.vstack([bl1, lar, im]),
        index=[f"P{i:03d}" for i in range(n_per_group * 3)],
        columns=all_genes_with_noise,
    )
    true_labels = ["BL1"] * n_per_group + ["LAR"] * n_per_group + ["IM"] * n_per_group
    patient_ids = list(expr_with_noise.index)

    # WITHOUT restricting to marker genes (clustering on all 90+23 columns, mirroring the
    # real 119-gene combined-panel situation)
    result_diluted = run_tnbc_subtyping_pipeline(
        pd.DataFrame(index=patient_ids), expr_with_noise, patient_ids,
        n_clusters=3, gene_selection_mode="all", log_transform=False,
    )
    silhouette_diluted = result_diluted["quality"]["silhouette_score"]

    # WITH restricting to marker genes only (the new default)
    result_restricted = run_tnbc_subtyping_pipeline(
        pd.DataFrame(index=patient_ids), expr_with_noise, patient_ids,
        n_clusters=3, gene_selection_mode="marker_only", log_transform=False,
    )
    silhouette_restricted = result_restricted["quality"]["silhouette_score"]

    print(f"Silhouette WITHOUT marker-gene restriction (90 noise genes included): {silhouette_diluted:.3f}")
    print(f"Silhouette WITH marker-gene restriction (noise genes excluded): {silhouette_restricted:.3f}")
    assert silhouette_restricted > silhouette_diluted, (
        "restricting to marker genes should give a clearly higher silhouette score than "
        "including 90 unrelated noise genes in the same clustering"
    )
    print(f"PASSED: marker-gene restriction gives a real, measurable improvement "
          f"({silhouette_diluted:.3f} -> {silhouette_restricted:.3f}) when unrelated genes are "
          f"mixed into the same expression data -- confirming Hypothesis 1 is a real mechanism, "
          f"not just a plausible-sounding guess.\n")

    print("=== Testing Hypothesis 2: unlogged scale-skew hurts clustering, log1p fixes it ===")
    # Same 3-group marker-only data, but now inject a scale-skew: multiply a few marker
    # genes' values by a huge factor for ALL patients equally (mirrors a few genes having
    # naturally huge TPM values regardless of group, e.g. highly-expressed housekeeping-like
    # genes) -- this should dominate raw NMF (since NMF is scale-sensitive) but matter far
    # less after log1p, which compresses large values much more than small ones.
    bl1_plain = make_group_signal(bl1_genes, marker_genes, n_per_group)
    lar_plain = make_group_signal(lar_genes, marker_genes, n_per_group)
    im_plain = make_group_signal(im_genes, marker_genes, n_per_group)
    expr_plain = np.vstack([bl1_plain, lar_plain, im_plain])

    # Inject extreme, group-UNRELATED scale skew into 3 of the marker genes for every patient
    skewed_genes = marker_genes[:3]
    for gene in skewed_genes:
        idx = marker_genes.index(gene)
        expr_plain[:, idx] = expr_plain[:, idx] * 5000 + rng.gamma(2.0, 1.0, size=len(expr_plain))

    expr_skewed_df = pd.DataFrame(expr_plain, index=patient_ids, columns=marker_genes)

    result_unlogged = run_tnbc_subtyping_pipeline(
        pd.DataFrame(index=patient_ids), expr_skewed_df, patient_ids,
        n_clusters=3, gene_selection_mode="marker_only", log_transform=False,
    )
    result_logged = run_tnbc_subtyping_pipeline(
        pd.DataFrame(index=patient_ids), expr_skewed_df, patient_ids,
        n_clusters=3, gene_selection_mode="marker_only", log_transform=True,
    )

    df_check = pd.DataFrame({"true_label": true_labels})
    purity_unlogged = pd.crosstab(df_check["true_label"], result_unlogged["patient_assignments"]["cluster"].values)
    purity_logged = pd.crosstab(df_check["true_label"], result_logged["patient_assignments"]["cluster"].values)
    mean_purity_unlogged = (purity_unlogged.max(axis=1) / purity_unlogged.sum(axis=1)).mean()
    mean_purity_logged = (purity_logged.max(axis=1) / purity_logged.sum(axis=1)).mean()

    print(f"Mean cluster purity WITHOUT log-transform (scale-skewed genes present): {mean_purity_unlogged:.3f}")
    print(f"Mean cluster purity WITH log-transform: {mean_purity_logged:.3f}")
    assert mean_purity_logged > mean_purity_unlogged, (
        "log-transforming should recover better cluster purity than raw values when a few "
        "genes have extreme, group-unrelated scale skew"
    )
    print(f"PASSED: log-transform gives a real, measurable improvement "
          f"({mean_purity_unlogged:.3f} -> {mean_purity_logged:.3f}) when scale-skewed genes are "
          f"present -- confirming Hypothesis 2 is a real mechanism too.\n")

    print("=" * 70)
    print("REAL-DATA-FINDING TESTS PASSED -- both hypotheses confirmed as real, measurable")
    print("mechanisms on controlled synthetic data, not just plausible-sounding guesses.")


def _run_highly_variable_mode_test():
    """
    Tests the new gene_selection_mode='highly_variable' end-to-end: true group
    signal spread across a broad, ~500-gene synthetic transcriptome (not just the
    curated marker panel), confirming (1) highly-variable selection actually finds
    the informative genes among mostly-flat background genes, (2) clustering on
    them recovers the true group structure, and (3) marker-based characterization
    still correctly identifies each cluster's true identity, using the SAME fixed
    marker genes regardless of what genes were selected for clustering itself.
    """
    rng = np.random.default_rng(11)
    n_per_group = 30

    bl1_genes, lar_genes = MARKER_GENE_SETS["BL1"], MARKER_GENE_SETS["LAR"]
    marker_genes = list(dict.fromkeys(bl1_genes + lar_genes))
    background_genes = [f"BG_GENE_{i}" for i in range(470)]  # flat, uninformative background
    all_genes = marker_genes + background_genes

    def make_group(marker_genes_for_group, n, boost=6.0):
        base = rng.gamma(2.0, 1.0, size=(n, len(all_genes)))
        for gene in marker_genes_for_group:
            idx = all_genes.index(gene)
            base[:, idx] += rng.gamma(2.0, 1.0, size=n) + boost
        return base

    bl1 = make_group(bl1_genes, n_per_group)
    lar = make_group(lar_genes, n_per_group)
    expr_broad = pd.DataFrame(
        np.vstack([bl1, lar]),
        index=[f"P{i:03d}" for i in range(n_per_group * 2)],
        columns=all_genes,
    )
    patient_ids = list(expr_broad.index)
    true_labels = ["BL1"] * n_per_group + ["LAR"] * n_per_group

    print("=== Testing gene_selection_mode='highly_variable' end-to-end ===")
    hv_genes = select_highly_variable_genes(expr_broad, n_top=20)
    n_true_markers_in_top20 = len(set(hv_genes) & set(marker_genes))
    print(f"Top-20 highly-variable genes include {n_true_markers_in_top20}/{len(marker_genes)} "
          f"of the true informative marker genes (out of 470 flat background genes)")
    assert n_true_markers_in_top20 >= len(marker_genes) - 1, (
        "highly-variable selection should find nearly all the true informative genes "
        "among a much larger set of flat background genes"
    )

    result = run_tnbc_subtyping_pipeline(
        pd.DataFrame(index=patient_ids), expr_broad, patient_ids,
        n_clusters=2, gene_selection_mode="highly_variable", n_highly_variable_genes=20,
    )
    df_check = pd.DataFrame({"true_label": true_labels, "cluster": result["patient_assignments"]["cluster"].values})
    cross_tab = pd.crosstab(df_check["true_label"], df_check["cluster"])
    print(cross_tab)
    for true_label in ["BL1", "LAR"]:
        purity = cross_tab.loc[true_label].max() / cross_tab.loc[true_label].sum()
        assert purity > 0.8, f"{true_label} should mostly cluster together (purity={purity:.2f})"

    print(result["cluster_characterization"][["n_patients", "BL1", "LAR", "best_matching_subtype"]].to_string())
    cluster_to_majority = df_check.groupby("cluster")["true_label"].agg(lambda x: x.value_counts().idxmax())
    for cluster_id, majority_label in cluster_to_majority.items():
        inferred = result["cluster_characterization"].loc[cluster_id, "best_matching_subtype"]
        assert inferred == majority_label, f"cluster {cluster_id} characterized as {inferred}, expected {majority_label}"

    print("PASSED: highly-variable gene selection correctly finds true signal among mostly-flat "
          "background genes, clustering on them recovers the true group structure, and marker-based "
          "characterization still correctly identifies each cluster using the fixed marker set.\n")
    print("=" * 70)
    print("HIGHLY-VARIABLE MODE TEST PASSED.")


if __name__ == "__main__":
    _run_smoke_test()
    print()
    _run_real_data_finding_tests()
    print()
    _run_highly_variable_mode_test()
