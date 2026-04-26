"""Predictor v2 — uses calibrated probabilities and OOD detection.

Drop-in replacement for predictor.py once model_v2.joblib exists.
Gracefully falls back to v1 model if v2 is missing.
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional

import numpy as np

from features import extract_features

MODEL_PATH_V2 = "/Users/bfentaw2/startup/nativeready/model/model_v2.joblib"
MODEL_PATH_V1 = "/Users/bfentaw2/startup/nativeready/model/model.joblib"
FEATURE_NAMES_PATH = "/Users/bfentaw2/startup/nativeready/model/feature_names.json"
TRAINED_VERSION = "0.2-calibrated"
FALLBACK_VERSION = "0.1-baseline"

_bundle: Optional[Dict[str, Any]] = None
_feature_names: Optional[List[str]] = None
_version: str = ""
_loaded = False


def _try_load() -> bool:
    global _bundle, _feature_names, _version, _loaded
    if _loaded:
        return _bundle is not None
    _loaded = True
    try:
        import joblib
        if os.path.exists(MODEL_PATH_V2):
            _bundle = joblib.load(MODEL_PATH_V2)
            _version = TRAINED_VERSION
        elif os.path.exists(MODEL_PATH_V1):
            _bundle = joblib.load(MODEL_PATH_V1)
            _version = FALLBACK_VERSION
        else:
            return False
        if isinstance(_bundle, dict) and "feature_names" in _bundle:
            _feature_names = _bundle["feature_names"]
        else:
            _feature_names = json.loads(open(FEATURE_NAMES_PATH).read())
        return True
    except Exception as e:
        print(f"[predictor] Failed to load model: {e}")
        _bundle = None
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
        {"name": "Sequence length", "value": f"{int(length)} amino acids",
         "risk_level": _risk_level(length, (50, 800), (20, 1500))},
        {"name": "Molecular weight", "value": f"{mw:.1f} kDa",
         "risk_level": _risk_level(mw, (5, 100), (2, 200))},
        {"name": "Hydrophobicity (GRAVY)", "value": f"{gravy:+.2f}",
         "risk_level": _risk_level(gravy, (-1.0, 0.2), (-1.5, 0.5))},
        {"name": "Isoelectric point", "value": f"{pi:.1f}",
         "risk_level": _risk_level(pi, (5.5, 8.5), (4.0, 10.0))},
        {"name": "Instability index", "value": f"{instability:.1f}",
         "risk_level": _risk_level(instability, (0, 40), (0, 60))},
        {"name": "Cysteine content", "value": f"{pct_cys:.1f}%",
         "risk_level": "low" if pct_cys >= 1.0 else "medium"},
    ]


def _recommendations(features: Dict[str, float], score: int, ood: bool) -> List[str]:
    """Per-risk-factor specific recommendations."""
    recs: List[str] = []
    pi = features["isoelectric_point"]
    gravy = features["gravy"]
    mw = features["molecular_weight_kda"]
    length = features["length"]
    pct_cys = features["pct_cysteine"]
    pct_charged = features["pct_charged"]

    # Buffer recommendation
    buffer_ph = max(5.0, min(8.5, round(pi, 1)))
    recs.append(f"Recommended buffer: 200 mM ammonium acetate, pH {buffer_ph:.1f}")

    # Hydrophobicity-specific
    if gravy > 0.5:
        recs.append("Highly hydrophobic — likely membrane protein. Use detergent screening (DDM, LMNG, OG) and titrate detergent concentration before native MS.")
    elif gravy > 0.3:
        recs.append("Moderately hydrophobic — consider mild detergent or amphipol if standard conditions yield poor signal.")

    # Size-specific
    if mw > 500:
        recs.append("Very large complex — use ultra-high-mass spectrometer (Orbitrap UHMR or Q-IM-ToF with extended mass range). Verify oligomeric state by SEC-MALS first.")
    elif mw > 150:
        recs.append("Large complex — verify oligomeric state by SEC-MALS prior to MS. Use high-mass mode if available.")

    # Length-specific
    if length < 50:
        recs.append("Small peptide — use low collision energies and smaller capillary tip. Standard intact-protein conditions should work, but the model's confidence on small peptides is limited.")
    elif length > 1500:
        recs.append("Long sequence — denaturing MS first to confirm sequence assignment is recommended before attempting native conditions.")

    # Stability cues
    if pct_cys < 0.5 and length > 100:
        recs.append("Low cysteine content — limited disulfide stabilization. Sample handling and buffer choice matter more.")
    elif pct_cys >= 2.0:
        recs.append("Cysteine-rich — disulfide bonds likely present. Use non-reducing conditions for native MS.")

    # pI-specific
    if pi > 9.0:
        recs.append("Very basic protein — consider negative-ion mode if positive-ion gives poor signal.")
    elif pi < 5.0:
        recs.append("Very acidic protein — positive-ion ESI may be challenging; try formic-acid-free buffer.")

    # OOD warning
    if ood:
        recs.append("Note: this sequence is unusual relative to the training data. Treat the score with extra caution and validate empirically.")

    # Score-based actions
    if score < 50:
        recs.append("Low predicted suitability — consider construct optimization, removing flexible termini, or chaperone co-expression to improve stability.")
    if score >= 80:
        recs.append("Standard native MS conditions should work. Start with the recommended buffer and standard collision energy ramps.")

    return recs


def _model_score(features: Dict[str, float]):
    assert _bundle is not None and _feature_names is not None
    model = _bundle["model"]
    X = np.array([[features[name] for name in _feature_names]], dtype=float)
    proba = model.predict_proba(X)[0]
    p = float(proba[1])
    score_int = int(round(100 * p))

    # OOD check (only available in v2 bundles)
    ood = False
    ood_distance = None
    if "ood_scaler" in _bundle and "ood_knn" in _bundle:
        scaler = _bundle["ood_scaler"]
        knn = _bundle["ood_knn"]
        threshold = _bundle["ood_threshold"]
        Xs = scaler.transform(X)
        dist, _ = knn.kneighbors(Xs)
        ood_distance = float(dist[0][-1])
        ood = ood_distance > threshold

    # Wider CI when OOD
    ci_half = 15 if ood else 9
    return {
        "suitability_score": score_int,
        "confidence_interval": {
            "lower": max(0, score_int - ci_half),
            "upper": min(100, score_int + ci_half),
        },
        "model_version": _version,
        "ood": ood,
        "ood_distance": ood_distance,
    }


def predict(sequence: str) -> Dict[str, Any]:
    features = extract_features(sequence)
    if not _try_load():
        raise RuntimeError("Model file not found.")
    core = _model_score(features)
    score_int = core["suitability_score"]
    risk_factors = _build_risk_factors(features)
    recs = _recommendations(features, score_int, core["ood"])

    response = {
        "suitability_score": score_int,
        "suitability_label": _label_for_score(score_int),
        "confidence_interval": core["confidence_interval"],
        "risk_factors": risk_factors,
        "recommendations": recs,
        "model_version": core["model_version"],
    }
    if core["ood"]:
        response["warning"] = "out_of_distribution"
        response["warning_message"] = "This sequence is unusual relative to the training data. Treat the score with extra caution."
    return response
