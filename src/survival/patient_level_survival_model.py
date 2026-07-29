"""
Patient-Level Survival / Prognosis Model (Track A)
====================================================

WHAT THIS ADDS THAT tnbc_specific_survival.py DOES NOT DO:
    tnbc_specific_survival.py answers a per-KINASE question: is this one
    kinase's expression associated with survival, tested independently
    for each of the 90 kinases, restricted to the confirmed-TNBC subset
    (~150-200 patients)?

    This module answers a per-PATIENT question instead: given a patient's
    combined kinase-panel expression plus clinical covariates, what is
    their individual predicted risk score? This is a genuine multivariate
    predictive model, not 90 independent univariate association tests,
    and it is designed to run across the FULL ~1095-1098 patient
    TCGA-BRCA cohort (all subtypes), with subtype included as a feature
    -- not as a pre-filter -- following the same design principle your
    Heuristic_vs_ML_Report.docx (Sec 5.1) already recommended for the
    DepMap dependency-prediction extension.

WHY A CUSTOM RIDGE-PENALIZED COX IMPLEMENTATION, NOT scikit-survival:
    scikit-survival / lifelines are not installed in this environment and
    could not be installed (no network access here). Rather than silently
    fall back to something that isn't actually a survival model (e.g.
    plain logistic regression ignoring censoring, which would misuse
    patients who are censored before an event -- a real correctness bug,
    not a stylistic choice), this module implements the real thing: an
    Efron-tie-corrected Cox partial-likelihood with an L2 (ridge) penalty,
    optimized directly. This is a legitimate, standard survival model,
    the same one CoxnetSurvivalAnalysis(l1_ratio=0.0) would give you.

    If scikit-survival becomes available later, swapping it in is a
    one-function change: pass a different `fit_fn`/`predict_fn` pair into
    train_and_evaluate() (the same dependency-injection pattern already
    used for run_survival_pipeline_fn in tnbc_specific_survival.py). True
    elastic net (L1+L2) is the natural next step once that swap is made,
    since a clean L1 path isn't practical to hand-roll here.

WHAT THIS DOES NOT DO (stated plainly, matching this project's existing
limitations sections):
    - No true elastic net (L1) here -- ridge (L2) only, for the reason above.
    - No treatment-response modeling. This TCGA-BRCA clinical file has
      survival/vital-status fields only, no pCR/RECIST/response field --
      confirmed directly, not assumed. Track B is intentionally NOT
      attempted here; see the separate response-data-gap writeup.
    - Feature set is deliberately restricted to the already-CTS-validated
      90-kinase panel plus a handful of clinical covariates, not the full
      transcriptome -- with n=~1000 and dozens of features this is a
      defensible regime; with n=~1000 and ~20,000 genes it would not be.
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


# =====================================================================
# 1. FEATURE MATRIX ASSEMBLY (full cohort, subtype as feature not filter)
# =====================================================================

def build_patient_feature_matrix(
    clinical_df: pd.DataFrame,
    expression_df: pd.DataFrame,
    kinase_panel: Sequence[str],
    subtype_column: str,
    clinical_covariates: Sequence[str] = ("age", "stage_numeric"),
) -> pd.DataFrame:
    """
    Assembles one row per patient across the FULL cohort (not restricted
    to TNBC): kinase-panel expression (already validated via CTS) +
    clinical covariates + one-hot subtype indicators (subtype as a
    feature, not a pre-filter).

    Both inputs are assumed indexed on patient_id, matching the
    convention already used by restrict_to_tnbc() in
    tnbc_specific_survival.py. Patients missing any kinase-panel
    expression value are dropped (not imputed) so that no row silently
    contributes a fabricated feature value to the model -- consistent
    with this project's stated preference to document real gaps rather
    than paper over them.
    """
    missing_kinases = [g for g in kinase_panel if g not in expression_df.columns]
    if missing_kinases:
        print(f"Note: {len(missing_kinases)} of {len(kinase_panel)} requested kinase-panel genes "
              f"not found in the expression data and will be excluded from the feature matrix: "
              f"{missing_kinases}")
    usable_kinases = [g for g in kinase_panel if g in expression_df.columns]

    expr = expression_df[usable_kinases]
    clinical = clinical_df[list(clinical_covariates) + [subtype_column]]

    joined = clinical.join(expr, how="inner")
    n_before = len(joined)
    joined = joined.dropna(subset=list(clinical_covariates) + usable_kinases)
    n_after = len(joined)
    if n_after < n_before:
        print(f"Note: dropped {n_before - n_after} patients with missing clinical or expression "
              f"values ({n_after} patients retained). Not imputed -- a missing value here is a "
              f"real data gap, not a neutral default.")

    subtype_dummies = pd.get_dummies(joined[subtype_column], prefix="subtype", drop_first=True)
    feature_df = pd.concat(
        [joined[list(clinical_covariates)], subtype_dummies, joined[usable_kinases]], axis=1
    )
    return feature_df


# =====================================================================
# 2. RIDGE-PENALIZED COX PROPORTIONAL HAZARDS (Efron ties)
# =====================================================================

def _efron_neg_log_partial_likelihood(
    beta: np.ndarray, X: np.ndarray, durations: np.ndarray, events: np.ndarray, ridge_lambda: float
) -> float:
    """
    Negative Efron-tie-corrected Cox partial log-likelihood plus an L2
    penalty. Efron's correction (rather than the simpler Breslow
    approximation) is used because TCGA follow-up times are recorded in
    whole days, so exact tied event times are common, not a rare edge case.
    """
    order = np.argsort(durations)
    durations, events, X = durations[order], events[order], X[order]
    risk_scores = X @ beta
    exp_scores = np.exp(risk_scores)

    unique_times = np.unique(durations[events == 1])
    log_lik = 0.0
    for t in unique_times:
        event_mask = (durations == t) & (events == 1)
        risk_mask = durations >= t
        d = event_mask.sum()
        sum_risk_events = risk_scores[event_mask].sum()
        sum_exp_risk_set = exp_scores[risk_mask].sum()
        sum_exp_events = exp_scores[event_mask].sum()

        log_lik += sum_risk_events
        for l in range(int(d)):
            log_lik -= np.log(sum_exp_risk_set - (l / d) * sum_exp_events)

    penalty = ridge_lambda * np.sum(beta ** 2)
    return -log_lik + penalty


def fit_ridge_cox(
    X: np.ndarray,
    durations: np.ndarray,
    events: np.ndarray,
    ridge_lambda: float = 1.0,
    standardize: bool = True,
) -> Dict:
    """
    Fits a ridge-penalized Cox model by direct minimization of the
    Efron partial likelihood. Returns fitted coefficients plus the
    feature means/scales used, so predict_risk() applies the identical
    transform to new data.
    """
    from scipy.optimize import minimize

    X = np.asarray(X, dtype=float)
    durations = np.asarray(durations, dtype=float)
    events = np.asarray(events, dtype=float)

    if standardize:
        means = X.mean(axis=0)
        scales = X.std(axis=0)
        scales[scales == 0] = 1.0
        X_fit = (X - means) / scales
    else:
        means = np.zeros(X.shape[1])
        scales = np.ones(X.shape[1])
        X_fit = X

    beta0 = np.zeros(X.shape[1])
    result = minimize(
        _efron_neg_log_partial_likelihood,
        beta0,
        args=(X_fit, durations, events, ridge_lambda),
        method="L-BFGS-B",
    )
    if not result.success:
        print(f"Note: Cox optimizer did not fully converge ({result.message}); "
              f"coefficients may be unstable -- consider raising ridge_lambda.")

    return {"beta": result.x, "means": means, "scales": scales, "converged": result.success}


def predict_risk(model: Dict, X_new: np.ndarray) -> np.ndarray:
    """Linear predictor (log relative hazard) for new patients, using the
    same standardization fitted on the training set."""
    X_new = np.asarray(X_new, dtype=float)
    X_scaled = (X_new - model["means"]) / model["scales"]
    return X_scaled @ model["beta"]


# =====================================================================
# 3. HARRELL'S C-INDEX (no external survival package required)
# =====================================================================

def concordance_index(durations: np.ndarray, events: np.ndarray, risk_scores: np.ndarray) -> float:
    """
    Standard Harrell's concordance index: among all comparable pairs
    (the patient with the shorter duration must have had an event),
    the fraction where the model correctly ranked the shorter-surviving
    patient as higher risk. O(n^2) -- fine at n~1000, not intended for
    much larger cohorts without a faster implementation.
    """
    durations = np.asarray(durations)
    events = np.asarray(events)
    risk_scores = np.asarray(risk_scores)
    n = len(durations)

    concordant, permissible, tied = 0, 0, 0
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if durations[i] < durations[j] and events[i] == 1:
                permissible += 1
                if risk_scores[i] > risk_scores[j]:
                    concordant += 1
                elif risk_scores[i] == risk_scores[j]:
                    tied += 1
    if permissible == 0:
        return float("nan")
    return (concordant + 0.5 * tied) / permissible


# =====================================================================
# 4. KAPLAN-MEIER ESTIMATOR (for risk-tertile stratification plots)
# =====================================================================

def kaplan_meier_estimate(durations: np.ndarray, events: np.ndarray) -> pd.DataFrame:
    """Standard product-limit estimator. Returns a step-function table
    (time, n_at_risk, n_events, survival_prob) suitable for plotting."""
    durations = np.asarray(durations)
    events = np.asarray(events)
    order = np.argsort(durations)
    durations, events = durations[order], events[order]

    unique_times = np.unique(durations)
    rows = []
    survival = 1.0
    for t in unique_times:
        n_at_risk = (durations >= t).sum()
        n_events = ((durations == t) & (events == 1)).sum()
        if n_at_risk > 0:
            survival *= (1 - n_events / n_at_risk)
        rows.append({"time": t, "n_at_risk": n_at_risk, "n_events": n_events, "survival_prob": survival})
    return pd.DataFrame(rows)


def stratify_by_risk(risk_scores: np.ndarray, n_groups: int = 3) -> np.ndarray:
    """Assigns each patient to a risk group by quantile (default:
    tertiles). Returns integer group labels, 0 = lowest predicted risk."""
    quantiles = np.quantile(risk_scores, np.linspace(0, 1, n_groups + 1))
    quantiles[-1] += 1e-9  # include the max value in the top bin
    return np.digitize(risk_scores, quantiles[1:-1])


# =====================================================================
# 5. TRAIN/EVALUATE WRAPPER (K-fold, dependency-injectable model)
# =====================================================================

def train_and_evaluate_survival_model(
    X: pd.DataFrame,
    durations: pd.Series,
    events: pd.Series,
    n_folds: int = 5,
    ridge_lambda: float = 1.0,
    fit_fn: Callable = fit_ridge_cox,
    predict_fn: Callable = predict_risk,
    random_state: int = 0,
) -> Dict:
    """
    K-fold cross-validated evaluation. fit_fn/predict_fn default to the
    ridge Cox implementation above but can be swapped for
    scikit-survival's CoxnetSurvivalAnalysis or RandomSurvivalForest
    (same dependency-injection pattern as run_survival_pipeline_fn in
    tnbc_specific_survival.py) without changing anything else here.

    Splits are patient-level by construction (one row per patient in X),
    consistent with this project's stated leakage-avoidance discipline
    elsewhere (e.g. the GroupKFold-by-cell-line approach in the DepMap
    ML work) -- there is no multi-sample-per-patient case in this
    clinical-covariate feature matrix to additionally guard against.
    """
    from sklearn.model_selection import KFold

    X_arr = X.to_numpy(dtype=float)
    durations_arr = durations.to_numpy(dtype=float)
    events_arr = events.to_numpy(dtype=float)

    kf = KFold(n_splits=n_folds, shuffle=True, random_state=random_state)
    fold_c_indices = []
    all_risk_scores = np.full(len(X_arr), np.nan)

    for fold_i, (train_idx, test_idx) in enumerate(kf.split(X_arr), 1):
        model = fit_fn(X_arr[train_idx], durations_arr[train_idx], events_arr[train_idx], ridge_lambda)
        risk_test = predict_fn(model, X_arr[test_idx])
        all_risk_scores[test_idx] = risk_test

        c_index = concordance_index(durations_arr[test_idx], events_arr[test_idx], risk_test)
        fold_c_indices.append(c_index)
        print(f"  Fold {fold_i}/{n_folds}: held-out C-index = {c_index:.3f} "
              f"(n_test={len(test_idx)}, n_events_test={int(events_arr[test_idx].sum())})")

    return {
        "fold_c_indices": fold_c_indices,
        "mean_c_index": float(np.nanmean(fold_c_indices)),
        "std_c_index": float(np.nanstd(fold_c_indices)),
        "out_of_fold_risk_scores": all_risk_scores,
    }


# =====================================================================
# 6. PERMUTATION IMPORTANCE (SHAP not installed in this environment)
# =====================================================================

def permutation_importance_cox(
    model: Dict, X: np.ndarray, durations: np.ndarray, events: np.ndarray, n_repeats: int = 10,
    random_state: int = 0,
) -> pd.DataFrame:
    """
    Feature importance via permutation: shuffle one feature column at a
    time, measure the drop in C-index. Not SHAP (not installed here, no
    network access to install it), but a legitimate, standard,
    model-agnostic importance measure that needs nothing beyond numpy.
    If shap becomes available later, this is a direct drop-in
    replacement point -- swap this function for a KernelExplainer/
    TreeExplainer call without touching anything upstream.
    """
    rng = np.random.default_rng(random_state)
    baseline_risk = predict_risk(model, X)
    baseline_c = concordance_index(durations, events, baseline_risk)

    rows = []
    for col_idx in range(X.shape[1]):
        drops = []
        for _ in range(n_repeats):
            X_permuted = X.copy()
            X_permuted[:, col_idx] = rng.permutation(X_permuted[:, col_idx])
            permuted_risk = predict_risk(model, X_permuted)
            permuted_c = concordance_index(durations, events, permuted_risk)
            drops.append(baseline_c - permuted_c)
        rows.append({"feature_index": col_idx, "mean_c_index_drop": np.mean(drops), "std": np.std(drops)})

    return pd.DataFrame(rows).sort_values("mean_c_index_drop", ascending=False).reset_index(drop=True)


# =====================================================================
# SMOKE TEST -- synthetic data with a known true risk structure,
# verifies the Cox fit recovers it and every metric behaves correctly
# =====================================================================

def _run_smoke_test():
    rng = np.random.default_rng(42)
    n = 400

    # Two informative covariates (true log-hazard signal) + two pure-noise covariates
    x_signal_1 = rng.normal(0, 1, n)
    x_signal_2 = rng.normal(0, 1, n)
    x_noise_1 = rng.normal(0, 1, n)
    x_noise_2 = rng.normal(0, 1, n)
    X = np.column_stack([x_signal_1, x_signal_2, x_noise_1, x_noise_2])

    true_log_hazard = 0.8 * x_signal_1 + 0.5 * x_signal_2
    baseline_hazard = 0.05
    event_times = rng.exponential(1.0 / (baseline_hazard * np.exp(true_log_hazard)))
    censor_times = rng.exponential(15.0, n)
    durations = np.minimum(event_times, censor_times)
    events = (event_times <= censor_times).astype(float)

    print(f"=== Synthetic cohort: n={n}, event rate={events.mean():.2f} ===\n")

    print("=== Testing fit_ridge_cox() + predict_risk() recovers the true signal direction ===")
    model = fit_ridge_cox(X, durations, events, ridge_lambda=1.0)
    risk = predict_risk(model, X)
    # Standardized coefficients for the two true signal features should be positive and
    # substantially larger in magnitude than the two pure-noise features.
    assert model["beta"][0] > 0 and model["beta"][1] > 0, "true signal coefficients should be positive"
    assert abs(model["beta"][0]) > abs(model["beta"][2]), "signal feature should outweigh a noise feature"
    assert abs(model["beta"][1]) > abs(model["beta"][3]), "signal feature should outweigh a noise feature"
    print(f"Fitted beta: {np.round(model['beta'], 3)} (features 0,1 = true signal; 2,3 = pure noise)")
    print("PASSED: signal coefficients are positive and larger in magnitude than noise coefficients.\n")

    print("=== Testing concordance_index() distinguishes a real model from random risk scores ===")
    c_real = concordance_index(durations, events, risk)
    c_random = concordance_index(durations, events, rng.permutation(risk))
    print(f"C-index (real model): {c_real:.3f} | C-index (shuffled/random risk): {c_random:.3f}")
    assert c_real > 0.65, "a real fitted model on this strong synthetic signal should clear 0.65"
    assert abs(c_random - 0.5) < 0.1, "randomly shuffled risk scores should be close to the 0.5 no-skill baseline"
    print("PASSED: C-index correctly separates a real fitted model from random risk scores.\n")

    print("=== Testing ridge_lambda actually shrinks coefficients ===")
    model_weak = fit_ridge_cox(X, durations, events, ridge_lambda=0.01)
    model_strong = fit_ridge_cox(X, durations, events, ridge_lambda=50.0)
    norm_weak = np.linalg.norm(model_weak["beta"])
    norm_strong = np.linalg.norm(model_strong["beta"])
    print(f"||beta|| at ridge_lambda=0.01: {norm_weak:.3f} | at ridge_lambda=50.0: {norm_strong:.3f}")
    assert norm_strong < norm_weak, "a much larger ridge penalty should shrink the coefficient vector"
    print("PASSED: stronger ridge penalty correctly produces smaller-magnitude coefficients.\n")

    print("=== Testing kaplan_meier_estimate() produces a monotonically non-increasing curve ===")
    km = kaplan_meier_estimate(durations, events)
    assert (km["survival_prob"].diff().dropna() <= 1e-9).all(), "KM survival probability must never increase"
    assert km["survival_prob"].iloc[0] <= 1.0
    print(f"KM curve has {len(km)} steps, final survival probability = {km['survival_prob'].iloc[-1]:.3f}")
    print("PASSED: KM curve is monotonically non-increasing, as required.\n")

    print("=== Testing stratify_by_risk() produces correctly ordered, roughly balanced tertiles ===")
    groups = stratify_by_risk(risk, n_groups=3)
    assert set(np.unique(groups)) == {0, 1, 2}
    mean_risk_by_group = [risk[groups == g].mean() for g in [0, 1, 2]]
    assert mean_risk_by_group[0] < mean_risk_by_group[1] < mean_risk_by_group[2], \
        "tertile 0 should have the lowest mean risk score, tertile 2 the highest"
    counts = np.bincount(groups)
    assert max(counts) - min(counts) <= n // 10, "tertile sizes should be roughly balanced"
    print(f"Tertile sizes: {counts}, mean risk by tertile: {np.round(mean_risk_by_group, 3)}")
    print("PASSED: risk tertiles are correctly ordered by mean risk and roughly balanced in size.\n")

    print("=== Testing train_and_evaluate_survival_model() end-to-end with K-fold CV ===")
    X_df = pd.DataFrame(X, columns=["signal_1", "signal_2", "noise_1", "noise_2"])
    cv_result = train_and_evaluate_survival_model(
        X_df, pd.Series(durations), pd.Series(events), n_folds=4, ridge_lambda=1.0
    )
    print(f"\nMean held-out C-index across folds: {cv_result['mean_c_index']:.3f} (+/- {cv_result['std_c_index']:.3f})")
    assert cv_result["mean_c_index"] > 0.6, "held-out CV C-index should still show real signal, not just in-sample fit"
    assert not np.isnan(cv_result["out_of_fold_risk_scores"]).any(), "every patient should get an out-of-fold risk score"
    print("PASSED: cross-validated C-index shows real held-out predictive signal, and every patient")
    print("received an out-of-fold risk score with no leakage-through-omission gaps.\n")

    print("=== Testing permutation_importance_cox() ranks true signal features above noise ===")
    importance = permutation_importance_cox(model, X, durations, events, n_repeats=15)
    top_two_features = set(importance.head(2)["feature_index"])
    assert top_two_features == {0, 1}, \
        f"the two true signal features (index 0, 1) should be the top-2 by importance, got {top_two_features}"
    print(importance)
    print("PASSED: permutation importance correctly ranks the two true signal features above the two")
    print("pure-noise features, using only numpy -- no SHAP dependency required.\n")

    print("=== Testing build_patient_feature_matrix() drops missing rows and one-hot-encodes subtype ===")
    clinical_df = pd.DataFrame({
        "age": [45, 60, 55, 70],
        "stage_numeric": [2, 3, 1, np.nan],  # last patient has a real missing value
        "subtype": ["Basal", "LumA", "Basal", "Her2"],
    }, index=["P1", "P2", "P3", "P4"])
    expression_df = pd.DataFrame({
        "EGFR": [5.1, 4.8, 6.0, 5.5], "ERBB2": [3.2, 7.9, 3.0, 8.1],
    }, index=["P1", "P2", "P3", "P4"])
    feat = build_patient_feature_matrix(clinical_df, expression_df, ["EGFR", "ERBB2"], "subtype")
    assert len(feat) == 3, "the one patient with a real missing stage value should be dropped, not imputed"
    assert "P4" not in feat.index
    assert "subtype_Her2" not in feat.columns or "P4" not in feat.index  # Her2 patient was dropped anyway
    assert any(c.startswith("subtype_") for c in feat.columns), "subtype should be one-hot encoded as a feature"
    print(f"Feature matrix shape: {feat.shape}, columns: {feat.columns.tolist()}")
    print("PASSED: missing values are dropped (not silently imputed), and subtype is correctly")
    print("included as a one-hot feature rather than used to pre-filter the cohort.\n")

    print("=" * 70)
    print("ALL SMOKE TESTS PASSED.")


if __name__ == "__main__":
    _run_smoke_test()
