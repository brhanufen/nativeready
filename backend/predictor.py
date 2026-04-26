"""Predictor for NativeReady.

Loads the trained scikit-learn model from disk. Uses the saved feature
ordering to guarantee features are passed to the model in the same order
as during training.

If the model file is missing or fails to load, falls back to a transparent
heuristic so the API still functions during early development.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional

import numpy as np

from features import extract_features

_HERE = os.path.dirname(os.path.abspath(__file__))
_MODEL_DIR = os.environ.get(
    "NATIVEREADY_MODEL_DIR",
    os.path.normpath(os.path.join(_HERE, "..", "model")),
)
MODEL_PATH = os.path.join(_MODEL_DIR, "model.joblib")
FEATURE_NAMES_PATH = os.path.join(_MODEL_DIR, "feature_names.json")
TRAINED_VERSION = "0.1-baseline"
HEURISTIC_VERSION = "0.1-baseline-heuristic"

_model_bundle: Optional[Dict[str, Any]] = None
_feature_names: Optional[List[str]] = None
_loaded = False


def _try_load() -> bool:
    """Load model + feature names. Returns True on success."""
    global _model_bundle, _feature_names, _loaded
    if _loaded:
        return _model_bundle is not None
    _loaded = True
    try:
        if not (os.path.exists(MODEL_PATH) and os.path.exists(FEATURE_NAMES_PATH)):
            return False
        import joblib
        _model_bundle = joblib.load(MODEL_PATH)
        _feature_names = json.loads(open(FEATURE_NAMES_PATH).read())
        return True
    except Exception as e:
        print(f"[predictor] Failed to load model: {e}")
        _model_bundle = None
        _feature_names = None
        return False


def _label_for_score(score: int) -> str:
    if score >= 80:
        return "Excellent"
    if score >= 65:
        return "Good"
    if score >= 50:
        return "Fair"
    if score >= 35:
        return "Poor"
    return "Unsuitable"


def _risk_level(value: float, low_ok: tuple, medium_ok: tuple) -> str:
    lo, hi = low_ok
    if lo <= value <= hi:
        return "low"
    mlo, mhi = medium_ok
    if mlo <= value <= mhi:
        return "medium"
    return "high"


def _build_risk_factors(features: Dict[str, float]) -> List[Dict[str, str]]:
    length = features["length"]
    mw = features["molecular_weight_kda"]
    gravy = features["gravy"]
    pi = features["isoelectric_point"]
    instability = features["instability_index"]
    pct_cys = features["pct_cysteine"]

    return [
        {
            "name": "Sequence length",
            "value": f"{int(length)} amino acids",
            "risk_level": _risk_level(length, (50, 800), (20, 1500)),
        },
        {
            "name": "Molecular weight",
            "value": f"{mw:.1f} kDa",
            "risk_level": _risk_level(mw, (5, 100), (2, 200)),
        },
        {
            "name": "Hydrophobicity (GRAVY)",
            "value": f"{gravy:+.2f}",
            "risk_level": _risk_level(gravy, (-1.0, 0.2), (-1.5, 0.5)),
        },
        {
            "name": "Isoelectric point",
            "value": f"{pi:.1f}",
            "risk_level": _risk_level(pi, (5.5, 8.5), (4.0, 10.0)),
        },
        {
            "name": "Instability index",
            "value": f"{instability:.1f}",
            "risk_level": _risk_level(instability, (0, 40), (0, 60)),
        },
        {
            "name": "Cysteine content",
            "value": f"{pct_cys:.1f}%",
            "risk_level": "low" if pct_cys >= 1.0 else "medium",
        },
    ]


def _recommendations(features: Dict[str, float], score: int) -> List[str]:
    recs: List[str] = []
    pi = features["isoelectric_point"]
    gravy = features["gravy"]
    mw = features["molecular_weight_kda"]
    length = features["length"]

    buffer_ph = max(5.0, min(8.5, round(pi, 1)))
    recs.append(f"Recommended buffer: 200 mM ammonium acetate, pH {buffer_ph:.1f}")

    if gravy > 0.3:
        recs.append("Hydrophobic sequence — consider detergent screening (DDM, LMNG) before native MS")
    if mw > 150:
        recs.append("Large complex — verify oligomeric state by SEC-MALS prior to MS")
    if length > 1000:
        recs.append("Long sequence — denaturing MS first to confirm sequence assignment is recommended")
    if features["pct_cysteine"] < 0.5:
        recs.append("Low cysteine content — limited disulfide stabilization; sample handling matters")
    if score < 50:
        recs.append("Low predicted suitability — consider construct optimization, removing flexible termini, or chaperone co-expression")
    if score >= 80:
        recs.append("Predicted suitable — standard native MS conditions should work")

    return recs


def _model_score(features: Dict[str, float]) -> Dict[str, Any]:
    """Use the trained model to score the sequence."""
    assert _model_bundle is not None and _feature_names is not None
    model = _model_bundle["model"]
    X = np.array([[features[name] for name in _feature_names]], dtype=float)
    proba = model.predict_proba(X)[0]
    p_suitable = float(proba[1])  # class 1 = native MS suitable
    score_int = int(round(100 * p_suitable))
    # CI based on training-time CV std (~5 percentage points)
    ci_half = 9
    return {
        "suitability_score": score_int,
        "confidence_interval": {
            "lower": max(0, score_int - ci_half),
            "upper": min(100, score_int + ci_half),
        },
        "model_version": TRAINED_VERSION,
    }


def _heuristic_score(features: Dict[str, float]) -> Dict[str, Any]:
    """Transparent fallback used only if model fails to load."""
    score = 70.0
    if 50 <= features["length"] <= 500:
        score += 10
    else:
        score -= 5
    if features["gravy"] > 0.5:
        score -= 15
    if abs(features["isoelectric_point"] - 7.0) <= 1.0:
        score += 5
    cys_count = (features["pct_cysteine"] / 100.0) * features["length"]
    if cys_count >= 2:
        score += min(5.0, cys_count * 0.5)

    score = max(0.0, min(100.0, score))
    score_int = int(round(score))
    return {
        "suitability_score": score_int,
        "confidence_interval": {
            "lower": max(0, score_int - 12),
            "upper": min(100, score_int + 12),
        },
        "model_version": HEURISTIC_VERSION,
    }


def predict(sequence: str) -> Dict[str, Any]:
    """Run prediction on a sequence and return the full response dict."""
    features = extract_features(sequence)

    if _try_load():
        core = _model_score(features)
    else:
        core = _heuristic_score(features)

    score_int = core["suitability_score"]
    return {
        "suitability_score": score_int,
        "suitability_label": _label_for_score(score_int),
        "confidence_interval": core["confidence_interval"],
        "risk_factors": _build_risk_factors(features),
        "recommendations": _recommendations(features, score_int),
        "model_version": core["model_version"],
    }
