from __future__ import annotations

import numpy as np


def rmse(observed: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.sqrt(np.mean((observed - predicted) ** 2)))


def pearson(observed: np.ndarray, predicted: np.ndarray) -> float:
    x, y = observed.ravel(), predicted.ravel()
    if np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def de_metrics(observed: np.ndarray, predicted: np.ndarray, top_k: int) -> dict[str, float]:
    observed_mean, predicted_mean = observed.mean(0), predicted.mean(0)
    k = min(top_k, observed.shape[1])
    obs_top = set(np.argsort(np.abs(observed_mean))[-k:])
    pred_top = set(np.argsort(np.abs(predicted_mean))[-k:])
    union = obs_top | pred_top
    sign = np.mean(np.sign(observed_mean[list(union)]) == np.sign(predicted_mean[list(union)]))
    return {
        "de_pearson": pearson(observed_mean, predicted_mean),
        f"de_top_{k}_overlap": len(obs_top & pred_top) / k,
        "de_sign_agreement": float(sign),
    }


def _pairwise_sq_dist(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.maximum((x * x).sum(1)[:, None] + (y * y).sum(1)[None, :] - 2 * x @ y.T, 0)


def energy_distance(observed: np.ndarray, predicted: np.ndarray) -> float:
    xy = np.sqrt(_pairwise_sq_dist(observed, predicted) + 1e-12).mean()
    xx = np.sqrt(_pairwise_sq_dist(observed, observed) + 1e-12).mean()
    yy = np.sqrt(_pairwise_sq_dist(predicted, predicted) + 1e-12).mean()
    return float(2 * xy - xx - yy)


def rbf_mmd(observed: np.ndarray, predicted: np.ndarray) -> float:
    joined = np.concatenate((observed, predicted))
    distances = _pairwise_sq_dist(joined, joined)
    positive = distances[distances > 0]
    bandwidth = float(np.median(positive)) if len(positive) else 1.0
    kernel = lambda a, b: np.exp(-_pairwise_sq_dist(a, b) / (2 * bandwidth + 1e-12))
    return float(kernel(observed, observed).mean() + kernel(predicted, predicted).mean()
                 - 2 * kernel(observed, predicted).mean())


def evaluate_arrays(observed: np.ndarray, predicted: np.ndarray, top_k: int) -> dict[str, float]:
    metrics = {
        "rmse": rmse(observed, predicted),
        "pearson": pearson(observed, predicted),
        "energy_distance": energy_distance(observed, predicted),
        "rbf_mmd": rbf_mmd(observed, predicted),
    }
    metrics.update(de_metrics(observed, predicted, top_k))
    return metrics

