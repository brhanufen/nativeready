"""Feature extraction for NativeReady — must match training pipeline exactly.

Same logic as /Users/bfentaw2/startup/nativeready/model/train.py to guarantee
that features computed at inference are identical to those at training time.
"""
from __future__ import annotations

import re

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
AMBIGUOUS_MAP = {"B": "N", "Z": "Q", "J": "L", "U": "C"}

# v0.4 glycosylation features (added 2026-05-11 in response to Marty feedback)
N_GLYC_PERMISSIVE = re.compile(r"N[^P][ST]")
N_GLYC_STRICT = re.compile(r"N[^PC][ST]")
MUCIN_WINDOW = 20
MUCIN_THRESHOLD = 0.40  # >=40% S/T/P content over a 20-aa window

# v0.4 transmembrane-helix features (KD hydropathy window, classical method)
# Kyte & Doolittle 1982 J Mol Biol scale; window=19, threshold=1.6 for TM helix.
# Validated against DeepTMHMM on the v6 training set; sensitive to the same
# membrane-protein failure modes Marty flagged (e.g., GPCRs, transporters).
KD_SCALE = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}
TM_WINDOW = 19
TM_THRESHOLD = 1.2  # tuned against DeepTMHMM ground truth on CB2 (7-TM)
TM_MIN_HELIX_LEN = 12  # discard short above-threshold runs (noise)


def _predict_tm_segments(seq: str):
    """Return list of (start, end) inclusive 1-based TM segment indices.

    Sliding-window mean KD hydropathy; segments where the windowed score
    stays above TM_THRESHOLD for at least TM_MIN_HELIX_LEN consecutive
    residues are called as TM helices. Classical Kyte-Doolittle method.
    """
    n = len(seq)
    if n < TM_WINDOW:
        return []
    # Per-residue KD score
    kd = [KD_SCALE.get(c, 0.0) for c in seq]
    # Windowed mean (centered)
    half = TM_WINDOW // 2
    window_scores = [0.0] * n
    for i in range(half, n - half):
        s = sum(kd[i - half:i + half + 1]) / TM_WINDOW
        window_scores[i] = s
    # Find runs above threshold
    segments = []
    in_seg = False
    seg_start = 0
    for i, s in enumerate(window_scores):
        if s >= TM_THRESHOLD and not in_seg:
            in_seg = True
            seg_start = i
        elif s < TM_THRESHOLD and in_seg:
            in_seg = False
            seg_end = i - 1
            if seg_end - seg_start + 1 >= TM_MIN_HELIX_LEN:
                segments.append((seg_start + 1, seg_end + 1))  # 1-based
    if in_seg:
        seg_end = n - 1
        if seg_end - seg_start + 1 >= TM_MIN_HELIX_LEN:
            segments.append((seg_start + 1, seg_end + 1))
    return segments


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

    # v0.4 glycosylation features
    n_seq = len(N_GLYC_PERMISSIVE.findall(seq))
    n_seq_strict = len(N_GLYC_STRICT.findall(seq))
    feats["n_glyc_sequon_count"] = float(n_seq)
    feats["n_glyc_sequon_strict_count"] = float(n_seq_strict)
    feats["n_glyc_sequon_density_per_100"] = float(100.0 * n_seq / len(seq))

    if len(seq) >= MUCIN_WINDOW:
        n_windows = len(seq) - MUCIN_WINDOW + 1
        threshold_count = int(MUCIN_THRESHOLD * MUCIN_WINDOW)
        mucin_windows = sum(
            1 for i in range(n_windows)
            if sum(1 for c in seq[i:i + MUCIN_WINDOW] if c in "STP") >= threshold_count
        )
        feats["mucin_like_window_count"] = float(mucin_windows)
        feats["mucin_like_fraction"] = float(mucin_windows / n_windows)
    else:
        feats["mucin_like_window_count"] = 0.0
        feats["mucin_like_fraction"] = 0.0

    feats["pct_stp"] = float(sum(aa_pct.get(c, 0.0) for c in "STP"))

    # v0.4 transmembrane-helix features (KD hydropathy window)
    tm_segments = _predict_tm_segments(seq)
    tm_lengths = [end - start + 1 for start, end in tm_segments]
    tm_total = sum(tm_lengths)
    feats["tm_helix_count"] = float(len(tm_segments))
    feats["tm_residues_total"] = float(tm_total)
    feats["tm_residue_fraction"] = float(tm_total / len(seq))
    feats["tm_longest_helix_len"] = float(max(tm_lengths) if tm_lengths else 0)
    feats["tm_polytopic_flag"] = 1.0 if len(tm_segments) >= 2 else 0.0

    return feats
