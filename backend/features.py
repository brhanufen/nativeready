"""Feature extraction for NativeReady — must match training pipeline exactly.

Same logic as /Users/bfentaw2/startup/nativeready/model/train.py to guarantee
that features computed at inference are identical to those at training time.
"""
from __future__ import annotations

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
AMBIGUOUS_MAP = {"B": "N", "Z": "Q", "J": "L", "U": "C"}


def clean_sequence(seq: str) -> str:
    """Strip whitespace, uppercase, replace ambiguous codes, drop unknowns."""
    seq = seq.upper().replace(" ", "").replace("\n", "").replace("\r", "").replace("\t", "")
    cleaned = []
    for c in seq:
        if c in STANDARD_AA:
            cleaned.append(c)
        elif c in AMBIGUOUS_MAP:
            cleaned.append(AMBIGUOUS_MAP[c])
        # X and other unknowns dropped
    return "".join(cleaned)


def extract_features(sequence: str) -> dict:
    """Compute the exact feature dict that the trained model expects.

    Returns a dict keyed by feature name. The caller must order them
    according to feature_names.json before passing to the model.
    """
    seq = clean_sequence(sequence)
    if len(seq) < 5:
        raise ValueError("Sequence too short to analyze (need at least 5 residues after cleaning)")

    pa = ProteinAnalysis(seq)
    aa_pct = pa.amino_acids_percent

    feats = {
        "length": float(len(seq)),
        "log_length": float(np.log1p(len(seq))),
        "molecular_weight_kda": float(pa.molecular_weight()) / 1000.0,
        "isoelectric_point": float(pa.isoelectric_point()),
        "instability_index": float(pa.instability_index()),
        "gravy": float(pa.gravy()),
        "aromaticity": float(pa.aromaticity()),
    }
    helix, turn, sheet = pa.secondary_structure_fraction()
    feats["secondary_structure_helix"] = float(helix)
    feats["secondary_structure_turn"] = float(turn)
    feats["secondary_structure_sheet"] = float(sheet)

    # Per-amino-acid composition (20 features)
    for aa in STANDARD_AA:
        feats[f"pct_{aa}"] = float(aa_pct.get(aa, 0.0))

    # Engineered grouped features (must match train.py exactly)
    feats["pct_cysteine"] = float(aa_pct.get("C", 0.0))
    feats["pct_proline"] = float(aa_pct.get("P", 0.0))
    feats["pct_charged"] = float(sum(aa_pct.get(c, 0.0) for c in "DEKR"))
    feats["pct_hydrophobic"] = float(sum(aa_pct.get(c, 0.0) for c in "AVILMFW"))
    feats["pct_aromatic"] = float(sum(aa_pct.get(c, 0.0) for c in "FWY"))
    feats["pct_polar"] = float(sum(aa_pct.get(c, 0.0) for c in "STNQ"))

    return feats
