"""Shared helpers for the ciprofloxacin resistance notebooks (03 modelling, 04 SHAP).

Feature names follow the format produced by notebook 02:
  - chromosomal amino-acid substitutions, e.g. ``gyrA_S83L`` or ``marR_Q12stop``
  - acquired resistance gene families, e.g. ``acq_qnrS`` or ``acq_aac(6')-Ib-cr``
"""
import re

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin
from sklearn.metrics import (balanced_accuracy_score, matthews_corrcoef,
                             recall_score, roc_auc_score)

# ---------------------------------------------------------------------------
# Feature categories
# ---------------------------------------------------------------------------

# Quinolone resistance-determining regions (QRDRs), E. coli K-12 amino-acid numbering.
# The ranges are deliberately slightly wider than the classic definitions so that
# substitutions reported just outside them (e.g. parE L416F) are included.
QRDR = {
    "gyrA": (67, 106),
    "gyrB": (426, 464),
    "parC": (56, 108),
    "parE": (416, 470),
}
EFFLUX_REGULATORS = ("marR", "acrR", "soxR")

SUBSTITUTION = re.compile(r"^(gyrA|gyrB|parC|parE|marR|acrR|soxR)_([A-Z])(\d+)([A-Z]|stop)$")

CATEGORY_LABELS = {
    "target_qrdr": "Target-site QRDR mutation",
    "acquired": "Acquired resistance gene",
    "efflux_truncation": "Efflux regulator truncation",
    "other": "Other substitution (possible lineage marker)",
}


def feature_category(name):
    """Classify a feature as a known resistance mechanism or 'other'.

    Known mechanisms are QRDR substitutions in the drug targets, acquired quinolone
    resistance genes, and premature stop codons in efflux regulators. Every other
    substitution is labelled 'other': these are mostly lineage-associated
    polymorphisms with no established role in resistance.
    """
    if name.startswith("acq_"):
        return "acquired"
    match = SUBSTITUTION.match(name)
    if not match:
        return "other"
    gene, _, position, alt = match.groups()
    if gene in QRDR:
        start, end = QRDR[gene]
        return "target_qrdr" if start <= int(position) <= end else "other"
    if gene in EFFLUX_REGULATORS and alt == "stop":
        return "efflux_truncation"
    return "other"


def known_mechanism_columns(columns):
    return [c for c in columns if feature_category(c) != "other"]


# ---------------------------------------------------------------------------
# Pipeline components
# ---------------------------------------------------------------------------

class MinCountFilter(BaseEstimator, TransformerMixin):
    """Drop binary features present (or absent) in fewer than ``min_count`` isolates.

    Used inside the model pipelines so the filter is learned from the training
    folds only, avoiding information leakage from the test folds.
    """

    def __init__(self, min_count=3):
        self.min_count = min_count

    def fit(self, X, y=None):
        X = pd.DataFrame(X)
        present = X.sum(axis=0)
        self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        self.keep_ = ((present >= self.min_count) &
                      (present <= len(X) - self.min_count)).to_numpy()
        return self

    def transform(self, X):
        X = pd.DataFrame(X, columns=self.feature_names_in_)
        return X.loc[:, self.keep_]

    def get_feature_names_out(self, input_features=None):
        return self.feature_names_in_[self.keep_]


class KnownMechanismRule(BaseEstimator, ClassifierMixin):
    """Rule-based baseline: count the known resistance mechanisms in each isolate.

    An isolate is predicted resistant when it carries at least ``threshold`` known
    mechanisms. The threshold (1-4) is chosen on the training data to maximise
    balanced accuracy, so the rule is fitted in the same way as the ML models.
    """

    def fit(self, X, y):
        X = pd.DataFrame(X)
        self.classes_ = np.array([0, 1])
        self.columns_ = known_mechanism_columns(X.columns)
        scores = self._count(X)
        self.threshold_ = max(range(1, 5), key=lambda t: balanced_accuracy_score(y, scores >= t))
        return self

    def _count(self, X):
        return pd.DataFrame(X)[self.columns_].sum(axis=1).to_numpy()

    def predict_proba(self, X):
        # Monotonic mapping of the count, with count == threshold giving 0.5
        p = np.clip(self._count(X) / (2 * self.threshold_), 0, 1)
        return np.column_stack([1 - p, p])

    def predict(self, X):
        return (self._count(X) >= self.threshold_).astype(int)


# ---------------------------------------------------------------------------
# Lineage clusters
# ---------------------------------------------------------------------------

def lineage_clusters(distances, threshold):
    """Cluster isolates by average-linkage hierarchical clustering of Mash distances.

    ``threshold`` is the Mash distance at which the tree is cut (0.01 is roughly
    99% average nucleotide identity).
    """
    d = distances.to_numpy(dtype=float)
    d = (d + d.T) / 2
    np.fill_diagonal(d, 0)
    tree = linkage(squareform(d, checks=False), method="average")
    return pd.Series(fcluster(tree, t=threshold, criterion="distance"),
                     index=distances.index, name="clade")


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

METRICS = ["AUROC", "balanced_accuracy", "MCC", "sensitivity", "specificity"]


def evaluate(y_true, prob, threshold=0.5):
    """AUROC, balanced accuracy, MCC, sensitivity and specificity for one set of predictions."""
    y_true = np.asarray(y_true)
    pred = (np.asarray(prob) >= threshold).astype(int)
    both_classes = len(np.unique(y_true)) == 2
    return {
        "AUROC": roc_auc_score(y_true, prob) if both_classes else np.nan,
        "balanced_accuracy": balanced_accuracy_score(y_true, pred) if both_classes else np.nan,
        "MCC": matthews_corrcoef(y_true, pred),
        "sensitivity": recall_score(y_true, pred, pos_label=1, zero_division=np.nan),
        "specificity": recall_score(y_true, pred, pos_label=0, zero_division=np.nan),
    }


def bootstrap_ci(y_true, prob, n_boot=1000, seed=42, alpha=0.05):
    """Point estimates with percentile bootstrap confidence intervals."""
    y_true, prob = np.asarray(y_true), np.asarray(prob)
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), len(y_true))
        if len(np.unique(y_true[idx])) == 2:
            samples.append(evaluate(y_true[idx], prob[idx]))
    samples = pd.DataFrame(samples)
    point = evaluate(y_true, prob)
    return pd.DataFrame({
        "estimate": point,
        "ci_low": samples.quantile(alpha / 2),
        "ci_high": samples.quantile(1 - alpha / 2),
    }).loc[METRICS]
