"""Predictor v3 — combined ESM-2 + BioPython model.

Architecture (matches train_v3_compare.py):
  1. Compute 36-dim BioPython physicochemical features
  2. Compute 1280-dim ESM-2 mean-pooled embedding
  3. esm_scaler.transform -> pca.transform -> 256-dim
  4. bio_scaler.transform -> 36-dim
  5. Concatenate -> 292-dim feature vector
  6. Calibrated RF predict_proba

OOD detector (separate bundle ood_detector_v3.joblib):
  - StandardScaler + KNN(k=5, Euclidean) on BioPython feature space
  - Threshold = 95th percentile within-training nearest-neighbor distance

Lazy-loads ESM-2 on first prediction to keep startup fast. Uses CPU; on Railway
hobby tier the first prediction takes ~3-7 seconds (model warm-up), subsequent
predictions are 2-4 seconds depending on sequence length.
"""
from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

import numpy as np

from features import extract_features

_HERE = os.path.dirname(os.path.abspath(__file__))
_MODEL_DIR = os.environ.get(
    "NATIVEREADY_MODEL_DIR",
    os.path.normpath(os.path.join(_HERE, "..", "model")),
)
MODEL_PATH_V3 = os.path.join(_MODEL_DIR, "model_v3.joblib")
OOD_PATH_V3 = os.path.join(_MODEL_DIR, "ood_detector_v3.joblib")
MODEL_PATH_V2 = os.path.join(_MODEL_DIR, "model_v2.joblib")  # graceful fallback
TRAINED_VERSION_V3 = "0.3-esm2-combined"
FALLBACK_VERSION_V2 = "0.2-calibrated"

ESM_MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
ESM_MAX_LEN = 1022  # 1024 model max minus 2 special tokens

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

_v3_bundle: Optional[Dict[str, Any]] = None
_v3_ood: Optional[Dict[str, Any]] = None
_v2_bundle: Optional[Dict[str, Any]] = None  # fallback
_esm_tokenizer = None
_esm_model = None
_esm_device = None
_loaded_classifier = False
_version: str = ""


def _try_load_classifier() -> bool:
    """Load v3 model bundle (preferred) or fall back to v2 if v3 missing."""
    global _v3_bundle, _v3_ood, _v2_bundle, _loaded_classifier, _version
    if _loaded_classifier:
        return _v3_bundle is not None or _v2_bundle is not None
    _loaded_classifier = True
    try:
        import joblib
        if os.path.exists(MODEL_PATH_V3):
            _v3_bundle = joblib.load(MODEL_PATH_V3)
            _version = TRAINED_VERSION_V3
            if os.path.exists(OOD_PATH_V3):
                _v3_ood = joblib.load(OOD_PATH_V3)
            print(f"[predictor_v3] Loaded v3 model from {MODEL_PATH_V3}")
            return True
        elif os.path.exists(MODEL_PATH_V2):
            _v2_bundle = joblib.load(MODEL_PATH_V2)
            _version = FALLBACK_VERSION_V2
            print(f"[predictor_v3] v3 missing, fell back to v2 from {MODEL_PATH_V2}")
            return True
        return False
    except Exception as e:
        print(f"[predictor_v3] Failed to load classifier bundle: {e}")
        return False


def _try_load_esm():
    """Lazy-load ESM-2 (heavy: ~2.5 GB weights)."""
    global _esm_tokenizer, _esm_model, _esm_device
    if _esm_model is not None:
        return True
    try:
        import torch
        from transformers import AutoTokenizer, AutoModel
        print(f"[predictor_v3] Loading ESM-2 model {ESM_MODEL_NAME} (one-time, ~30s)...")
        _esm_tokenizer = AutoTokenizer.from_pretrained(ESM_MODEL_NAME)
        _esm_model = AutoModel.from_pretrained(ESM_MODEL_NAME)
        _esm_model.eval()
        if torch.backends.mps.is_available():
            _esm_device = torch.device("mps")
        elif torch.cuda.is_available():
            _esm_device = torch.device("cuda")
        else:
            _esm_device = torch.device("cpu")
        _esm_model = _esm_model.to(_esm_device)
        print(f"[predictor_v3] ESM-2 ready on {_esm_device}")
        return True
    except Exception as e:
        print(f"[predictor_v3] Failed to load ESM-2: {e}")
        return False


def _clean_seq(seq: str) -> str:
    seq = seq.upper().replace(" ", "").replace("\n", "")
    cleaned = []
    for c in seq:
        if c in STANDARD_AA:
            cleaned.append(c)
        elif c in "BZJU":
            cleaned.append({"B": "N", "Z": "Q", "J": "L", "U": "C"}[c])
    return "".join(cleaned)


def _esm_embed(sequence: str) -> np.ndarray:
    """Compute mean-pooled ESM-2 embedding for one sequence (1280-dim)."""
    import torch
    seq = _clean_seq(sequence)[:ESM_MAX_LEN]
    with torch.no_grad():
        inputs = _esm_tokenizer(seq, return_tensors="pt",
                                truncation=True, max_length=ESM_MAX_LEN + 2).to(_esm_device)
        outputs = _esm_model(**inputs)
        hidden = outputs.last_hidden_state[0, 1:-1, :]  # drop CLS and EOS
        emb = hidden.mean(dim=0).cpu().numpy().astype(np.float32)
    return emb


def _label_for_score(score: int) -> str:
    if score >= 80: return "Excellent"
    if score >= 65: return "Good"
    if score >= 50: return "Fair"
    if score >= 35: return "Poor"
    return "Unsuitable"


def _risk_level(value: float, low_ok, medium_ok) -> str:
    lo, hi = low_ok
    if lo <= value <= hi:
        return "low"
    mlo, mhi = medium_ok
    if mlo <= value <= mhi:
        return "medium"
    return "high"


def _build_risk_factors(features: Dict[str, float]) -> List[Dict[str, str]]:
    return [
        {"name": "Sequence length", "value": f"{int(features['length'])} amino acids",
         "risk_level": _risk_level(features['length'], (50, 800), (20, 1500))},
        {"name": "Molecular weight", "value": f"{features['molecular_weight_kda']:.1f} kDa",
         "risk_level": _risk_level(features['molecular_weight_kda'], (5, 100), (2, 200))},
        {"name": "Hydrophobicity (GRAVY)", "value": f"{features['gravy']:+.2f}",
         "risk_level": _risk_level(features['gravy'], (-1.0, 0.2), (-1.5, 0.5))},
        {"name": "Isoelectric point", "value": f"{features['isoelectric_point']:.1f}",
         "risk_level": _risk_level(features['isoelectric_point'], (5.5, 8.5), (4.0, 10.0))},
        {"name": "Instability index", "value": f"{features['instability_index']:.1f}",
         "risk_level": _risk_level(features['instability_index'], (0, 40), (0, 60))},
        {"name": "Cysteine content", "value": f"{features['pct_cysteine']:.1f}%",
         "risk_level": "low" if features['pct_cysteine'] >= 1.0 else "medium"},
    ]


def _recommendations(features: Dict[str, float], score: int, ood: bool) -> List[str]:
    """Same recommendation logic as v2 (physicochemical-based, model-agnostic)."""
    recs: List[str] = []
    pi, gravy, mw = features['isoelectric_point'], features['gravy'], features['molecular_weight_kda']
    length, pct_cys = features['length'], features['pct_cysteine']
    buffer_ph = max(5.0, min(8.5, round(pi, 1)))
    recs.append(f"Recommended buffer: 200 mM ammonium acetate, pH {buffer_ph:.1f}")
    if gravy > 0.5:
        recs.append("Highly hydrophobic — likely membrane protein. Use detergent screening (DDM, LMNG, OG) and titrate detergent concentration before native MS.")
    elif gravy > 0.3:
        recs.append("Moderately hydrophobic — consider mild detergent or amphipol if standard conditions yield poor signal.")
    if mw > 500:
        recs.append("Very large complex — use ultra-high-mass spectrometer (Orbitrap UHMR or Q-IM-ToF with extended mass range). Verify oligomeric state by SEC-MALS first.")
    elif mw > 150:
        recs.append("Large complex — verify oligomeric state by SEC-MALS prior to MS. Use high-mass mode if available.")
    if length < 50:
        recs.append("Small peptide — use low collision energies and smaller capillary tip. Standard intact-protein conditions should work.")
    elif length > 1500:
        recs.append("Long sequence — denaturing MS first to confirm sequence assignment is recommended before attempting native conditions.")
    if pct_cys < 0.5 and length > 100:
        recs.append("Low cysteine content — limited disulfide stabilization.")
    elif pct_cys >= 2.0:
        recs.append("Cysteine-rich — disulfide bonds likely present. Use non-reducing conditions.")
    if pi > 9.0:
        recs.append("Very basic protein — consider negative-ion mode if positive-ion gives poor signal.")
    elif pi < 5.0:
        recs.append("Very acidic protein — positive-ion ESI may be challenging.")
    if ood:
        recs.append("Note: this sequence is unusual relative to the training data. Treat the score with extra caution and validate empirically.")
    if score < 50:
        recs.append("Low predicted suitability — consider construct optimization, removing flexible termini, or chaperone co-expression.")
    if score >= 80:
        recs.append("Standard native MS conditions should work. Start with the recommended buffer and standard collision energy ramps.")
    return recs


def _v3_score(features: Dict[str, float], sequence: str) -> Dict[str, Any]:
    """Run the v3 combined model (BioPython + PCA-256 ESM-2)."""
    assert _v3_bundle is not None
    if not _try_load_esm():
        raise RuntimeError("ESM-2 unavailable; cannot run v0.3 model.")
    # Features
    bio_vec = np.array([[features[k] for k in [
        "length", "log_length", "molecular_weight_kda", "isoelectric_point",
        "instability_index", "gravy", "aromaticity",
        "secondary_structure_helix", "secondary_structure_turn", "secondary_structure_sheet",
        *[f"pct_{aa}" for aa in STANDARD_AA],
        "pct_cysteine", "pct_proline", "pct_charged", "pct_hydrophobic",
        "pct_aromatic", "pct_polar"
    ]]], dtype=float)
    esm_vec = _esm_embed(sequence).reshape(1, -1)
    # Apply saved transforms
    esm_scaled = _v3_bundle["esm_scaler"].transform(esm_vec)
    esm_pca = _v3_bundle["pca"].transform(esm_scaled)
    bio_scaled = _v3_bundle["bio_scaler"].transform(bio_vec)
    X = np.hstack([esm_pca, bio_scaled])
    proba = _v3_bundle["model"].predict_proba(X)[0]
    p = float(proba[1])
    score_int = int(round(100 * p))
    # OOD using v3 detector on raw BioPython features
    ood, ood_distance = False, None
    if _v3_ood is not None:
        Xs = _v3_ood["ood_scaler"].transform(bio_vec)
        dist, _ = _v3_ood["ood_knn"].kneighbors(Xs)
        ood_distance = float(dist[0][-1])
        ood = ood_distance > _v3_ood["ood_threshold"]
    ci_half = 15 if ood else 9
    return {
        "suitability_score": score_int,
        "confidence_interval": {"lower": max(0, score_int - ci_half),
                                 "upper": min(100, score_int + ci_half)},
        "model_version": _version,
        "ood": ood,
        "ood_distance": ood_distance,
    }


def _v2_score(features: Dict[str, float]) -> Dict[str, Any]:
    """Fallback: v2 BioPython-only model if v3/ESM not available."""
    assert _v2_bundle is not None
    import json
    feat_names = _v2_bundle.get("feature_names")
    if not feat_names:
        feat_names = json.load(open(os.path.join(_MODEL_DIR, "feature_names.json")))
    X = np.array([[features[name] for name in feat_names]], dtype=float)
    proba = _v2_bundle["model"].predict_proba(X)[0]
    p = float(proba[1])
    score_int = int(round(100 * p))
    ood, ood_distance = False, None
    if "ood_scaler" in _v2_bundle and "ood_knn" in _v2_bundle:
        Xs = _v2_bundle["ood_scaler"].transform(X)
        dist, _ = _v2_bundle["ood_knn"].kneighbors(Xs)
        ood_distance = float(dist[0][-1])
        ood = ood_distance > _v2_bundle["ood_threshold"]
    ci_half = 15 if ood else 9
    return {
        "suitability_score": score_int,
        "confidence_interval": {"lower": max(0, score_int - ci_half),
                                 "upper": min(100, score_int + ci_half)},
        "model_version": _version,
        "ood": ood,
        "ood_distance": ood_distance,
    }


def predict(sequence: str) -> Dict[str, Any]:
    if not _try_load_classifier():
        raise RuntimeError("No model file found (looked for v3 then v2).")
    features = extract_features(sequence)
    if _v3_bundle is not None:
        try:
            core = _v3_score(features, sequence)
        except Exception as e:
            print(f"[predictor_v3] v3 inference failed ({e}), falling back to v2")
            if _v2_bundle is None:
                # No fallback; re-raise
                raise
            core = _v2_score(features)
    else:
        core = _v2_score(features)
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
