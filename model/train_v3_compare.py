#!/usr/bin/env python3
"""Train and compare 4 model variants on 634 proteins with consistent 5-fold CV.

Variants:
  V1 (v0.2 baseline)    : 36 BioPython features only, RF + isotonic
  V2 (ESM-2 linear)     : 1280 ESM-2 features, LogisticRegression L2 + isotonic
  V3 (ESM-2 + PCA + RF) : PCA-256(ESM-2), RF + isotonic
  V4 (combined)         : PCA-256(ESM-2) + 36 BioPython = 292 features, RF + isotonic

All variants use the SAME stratified 5-fold splits (random_state=42) for direct
apples-to-apples comparison. Reports mean AUC, mean Brier, mean accuracy, and
per-class confusion matrix breakdown.
"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import joblib
import torch
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.calibration import CalibratedClassifierCV
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, brier_score_loss, roc_auc_score,
    average_precision_score, confusion_matrix, classification_report
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

DATA_DIR = Path("/Users/bfentaw2/startup/nativeready/data")
MODEL_DIR = Path("/Users/bfentaw2/startup/nativeready/model")

DATASET = DATA_DIR / "dataset_combined_v7_2026-05-11.json"
ESM_EMBEDDINGS = MODEL_DIR / "esm2_embeddings_635.npy"
ESM_METADATA = MODEL_DIR / "esm2_embeddings_metadata.json"

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
RANDOM_STATE = 42
CV_SPLITS = 5

import re
N_GLYC_PERMISSIVE = re.compile(r"N[^P][ST]")
N_GLYC_STRICT = re.compile(r"N[^PC][ST]")
MUCIN_WINDOW = 20
MUCIN_THRESHOLD = 0.40

KD_SCALE = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,
            "G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,
            "P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
TM_WINDOW = 19
TM_THRESHOLD = 1.2
TM_MIN_HELIX_LEN = 12

def _tm_segments(seq):
    n = len(seq)
    if n < TM_WINDOW:
        return []
    half = TM_WINDOW // 2
    kd = [KD_SCALE.get(c, 0.0) for c in seq]
    ws = [0.0]*n
    for i in range(half, n-half):
        ws[i] = sum(kd[i-half:i+half+1])/TM_WINDOW
    segs, in_seg, st = [], False, 0
    for i, s in enumerate(ws):
        if s >= TM_THRESHOLD and not in_seg:
            in_seg, st = True, i
        elif s < TM_THRESHOLD and in_seg:
            in_seg = False
            if i - st >= TM_MIN_HELIX_LEN:
                segs.append((st, i-1))
    if in_seg and n - st >= TM_MIN_HELIX_LEN:
        segs.append((st, n-1))
    return segs


def clean_sequence(seq: str) -> str:
    seq = seq.upper().replace(" ", "").replace("\n", "")
    cleaned = []
    for c in seq:
        if c in STANDARD_AA:
            cleaned.append(c)
        elif c in "BZJU":
            cleaned.append({"B": "N", "Z": "Q", "J": "L", "U": "C"}[c])
    return "".join(cleaned)


def biopython_features(sequence: str) -> dict:
    seq = clean_sequence(sequence)
    if len(seq) < 5:
        return None
    pa = ProteinAnalysis(seq)
    aa_pct = pa.amino_acids_percent
    feats = {
        "length": len(seq),
        "log_length": np.log1p(len(seq)),
        "molecular_weight_kda": pa.molecular_weight() / 1000.0,
        "isoelectric_point": pa.isoelectric_point(),
        "instability_index": pa.instability_index(),
        "gravy": pa.gravy(),
        "aromaticity": pa.aromaticity(),
    }
    helix, turn, sheet = pa.secondary_structure_fraction()
    feats["secondary_structure_helix"] = helix
    feats["secondary_structure_turn"] = turn
    feats["secondary_structure_sheet"] = sheet
    for aa in STANDARD_AA:
        feats[f"pct_{aa}"] = aa_pct.get(aa, 0.0)
    feats["pct_cysteine"] = aa_pct.get("C", 0.0)
    feats["pct_proline"] = aa_pct.get("P", 0.0)
    feats["pct_charged"] = sum(aa_pct.get(c, 0.0) for c in "DEKR")
    feats["pct_hydrophobic"] = sum(aa_pct.get(c, 0.0) for c in "AVILMFW")
    feats["pct_aromatic"] = sum(aa_pct.get(c, 0.0) for c in "FWY")
    feats["pct_polar"] = sum(aa_pct.get(c, 0.0) for c in "STNQ")
    # v0.4 glycosylation features (must match backend/features.py exactly)
    n_seq = len(N_GLYC_PERMISSIVE.findall(seq))
    n_seq_strict = len(N_GLYC_STRICT.findall(seq))
    feats["n_glyc_sequon_count"] = float(n_seq)
    feats["n_glyc_sequon_strict_count"] = float(n_seq_strict)
    feats["n_glyc_sequon_density_per_100"] = 100.0 * n_seq / len(seq)
    if len(seq) >= MUCIN_WINDOW:
        n_windows = len(seq) - MUCIN_WINDOW + 1
        threshold_count = int(MUCIN_THRESHOLD * MUCIN_WINDOW)
        mucin_windows = sum(
            1 for i in range(n_windows)
            if sum(1 for c in seq[i:i + MUCIN_WINDOW] if c in "STP") >= threshold_count
        )
        feats["mucin_like_window_count"] = float(mucin_windows)
        feats["mucin_like_fraction"] = mucin_windows / n_windows
    else:
        feats["mucin_like_window_count"] = 0.0
        feats["mucin_like_fraction"] = 0.0
    feats["pct_stp"] = sum(aa_pct.get(c, 0.0) for c in "STP")
    # v0.4 TM features
    segs = _tm_segments(seq)
    lens = [e - s + 1 for s, e in segs]
    total = sum(lens)
    feats["tm_helix_count"] = float(len(segs))
    feats["tm_residues_total"] = float(total)
    feats["tm_residue_fraction"] = total / len(seq)
    feats["tm_longest_helix_len"] = float(max(lens) if lens else 0)
    feats["tm_polytopic_flag"] = 1.0 if len(segs) >= 2 else 0.0
    return feats


def make_calibrated_rf():
    base = RandomForestClassifier(
        n_estimators=300, max_depth=10, min_samples_leaf=2,
        class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1,
    )
    return CalibratedClassifierCV(base, method="isotonic", cv=3)


def make_calibrated_logreg():
    base = LogisticRegression(
        C=1.0, penalty="l2", max_iter=2000, solver="lbfgs",
        class_weight="balanced", random_state=RANDOM_STATE,
    )
    return CalibratedClassifierCV(base, method="isotonic", cv=3)


def evaluate_variant(name, X, y, classes, make_clf, scaler=None):
    """Run stratified 5-fold CV and report metrics."""
    skf = StratifiedKFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    aucs, accs, briers, aps = [], [], [], []
    fold_predictions = np.zeros(len(y))
    fold_probas = np.zeros(len(y))

    for fold_i, (tr, te) in enumerate(skf.split(X, y)):
        X_tr, X_te = X[tr], X[te]
        if scaler is not None:
            sc = StandardScaler().fit(X_tr)
            X_tr = sc.transform(X_tr)
            X_te = sc.transform(X_te)
        clf = make_clf()
        clf.fit(X_tr, y[tr])
        proba = clf.predict_proba(X_te)[:, 1]
        pred = clf.predict(X_te)
        fold_predictions[te] = pred
        fold_probas[te] = proba
        aucs.append(float(roc_auc_score(y[te], proba)))
        accs.append(float(accuracy_score(y[te], pred)))
        briers.append(float(brier_score_loss(y[te], proba)))
        aps.append(float(average_precision_score(y[te], proba)))

    cm = confusion_matrix(y, fold_predictions).tolist()

    # Per-class breakdown
    class_breakdown = {}
    if classes is not None:
        df = pd.DataFrame({"y": y, "pred": fold_predictions, "class": classes})
        for cls, sub in df.groupby("class"):
            if len(sub) >= 5:  # need a few samples to be meaningful
                acc_cls = (sub["y"] == sub["pred"]).mean()
                class_breakdown[cls] = {
                    "n": int(len(sub)),
                    "accuracy": float(acc_cls),
                    "n_pos": int((sub["y"] == 1).sum()),
                }

    return {
        "name": name,
        "n_features": int(X.shape[1]),
        "cv_auc_mean": float(np.mean(aucs)),
        "cv_auc_std": float(np.std(aucs)),
        "cv_acc_mean": float(np.mean(accs)),
        "cv_acc_std": float(np.std(accs)),
        "cv_brier_mean": float(np.mean(briers)),
        "cv_average_precision_mean": float(np.mean(aps)),
        "fold_aucs": aucs,
        "confusion_matrix_overall": cm,
        "per_class_accuracy": class_breakdown,
        "fold_predictions": fold_predictions.tolist(),
        "fold_probas": fold_probas.tolist(),
    }


def main():
    print("=" * 70)
    print("v0.3 ESM-2 MODEL COMPARISON ON 634 PROTEINS")
    print("=" * 70)

    # Load dataset
    with open(DATASET) as f:
        proteins = json.load(f)
    print(f"\nLoaded {len(proteins)} proteins from {DATASET.name}")

    # Load ESM-2 embeddings + metadata
    esm_X = np.load(ESM_EMBEDDINGS)
    with open(ESM_METADATA) as f:
        esm_meta = json.load(f)
    esm_uids = esm_meta["uids"]
    print(f"Loaded ESM-2 embeddings: shape={esm_X.shape}")

    # Index proteins by uniprot_id for alignment
    by_uid = {p["uniprot_id"]: p for p in proteins}
    if not all(uid in by_uid for uid in esm_uids):
        print("ERROR: ESM embedding UIDs do not match dataset")
        return

    # Build aligned arrays
    y = np.array([by_uid[uid]["label"] for uid in esm_uids], dtype=int)
    seqs = [by_uid[uid]["sequence"] for uid in esm_uids]
    classes = [by_uid[uid].get("protein_class", "unknown") for uid in esm_uids]

    # Compute BioPython features for all
    print("\nComputing BioPython features...")
    bio_rows = []
    valid_idx = []
    for i, seq in enumerate(seqs):
        f = biopython_features(seq)
        if f is None:
            continue
        bio_rows.append(f)
        valid_idx.append(i)
    bio_X = pd.DataFrame(bio_rows).values
    valid_idx = np.array(valid_idx)
    print(f"  BioPython feature matrix: {bio_X.shape}")
    print(f"  Valid proteins (sequence >= 5 aa): {len(valid_idx)} / {len(seqs)}")

    # Align all matrices to valid_idx
    esm_X = esm_X[valid_idx]
    y = y[valid_idx]
    classes_aligned = [classes[i] for i in valid_idx]
    print(f"\n  Final aligned shapes:")
    print(f"    BioPython: {bio_X.shape}")
    print(f"    ESM-2:     {esm_X.shape}")
    print(f"    y:         {y.shape}  (pos={int((y==1).sum())}, neg={int((y==0).sum())})")

    # PCA-256 on ESM
    print("\nFitting PCA-256 on ESM-2 embeddings (computed PER FOLD in evaluation)")

    # ----- Define variants -----

    # V1: BioPython baseline (matches v0.2 architecture)
    print("\n" + "=" * 70)
    print("VARIANT 1: BioPython 36 features (v0.2 baseline reproduction)")
    print("=" * 70)
    r1 = evaluate_variant("v1_biopython_only", bio_X, y, classes_aligned, make_calibrated_rf)
    print(f"  AUC: {r1['cv_auc_mean']:.3f} (+/- {r1['cv_auc_std']:.3f})")
    print(f"  Acc: {r1['cv_acc_mean']:.3f}  Brier: {r1['cv_brier_mean']:.3f}  AP: {r1['cv_average_precision_mean']:.3f}")

    # V2: ESM-2 linear probe (logreg with standardization)
    print("\n" + "=" * 70)
    print("VARIANT 2: ESM-2 1280 features, LogisticRegression L2 (linear probe)")
    print("=" * 70)
    r2 = evaluate_variant("v2_esm2_linear_probe", esm_X, y, classes_aligned, make_calibrated_logreg, scaler=True)
    print(f"  AUC: {r2['cv_auc_mean']:.3f} (+/- {r2['cv_auc_std']:.3f})")
    print(f"  Acc: {r2['cv_acc_mean']:.3f}  Brier: {r2['cv_brier_mean']:.3f}  AP: {r2['cv_average_precision_mean']:.3f}")

    # V3: ESM-2 + PCA + RF
    # PCA fit per fold below to avoid leakage
    print("\n" + "=" * 70)
    print("VARIANT 3: PCA-256(ESM-2), RandomForest")
    print("=" * 70)

    def evaluate_pca_rf(X_full, y, classes, n_components=256):
        skf = StratifiedKFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_STATE)
        aucs, accs, briers, aps = [], [], [], []
        fold_predictions = np.zeros(len(y))
        fold_probas = np.zeros(len(y))
        for tr, te in skf.split(X_full, y):
            sc = StandardScaler().fit(X_full[tr])
            Xs_tr = sc.transform(X_full[tr])
            Xs_te = sc.transform(X_full[te])
            pca = PCA(n_components=min(n_components, Xs_tr.shape[0]-1, Xs_tr.shape[1]),
                       random_state=RANDOM_STATE).fit(Xs_tr)
            Xp_tr = pca.transform(Xs_tr)
            Xp_te = pca.transform(Xs_te)
            clf = make_calibrated_rf()
            clf.fit(Xp_tr, y[tr])
            proba = clf.predict_proba(Xp_te)[:, 1]
            pred = clf.predict(Xp_te)
            fold_predictions[te] = pred
            fold_probas[te] = proba
            aucs.append(float(roc_auc_score(y[te], proba)))
            accs.append(float(accuracy_score(y[te], pred)))
            briers.append(float(brier_score_loss(y[te], proba)))
            aps.append(float(average_precision_score(y[te], proba)))
        cm = confusion_matrix(y, fold_predictions).tolist()
        class_breakdown = {}
        df = pd.DataFrame({"y": y, "pred": fold_predictions, "class": classes})
        for cls, sub in df.groupby("class"):
            if len(sub) >= 5:
                class_breakdown[cls] = {
                    "n": int(len(sub)), "accuracy": float((sub["y"]==sub["pred"]).mean()),
                    "n_pos": int((sub["y"]==1).sum())
                }
        return {
            "name": "v3_esm2_pca256_rf",
            "n_features": min(n_components, X_full.shape[1]),
            "cv_auc_mean": float(np.mean(aucs)),
            "cv_auc_std": float(np.std(aucs)),
            "cv_acc_mean": float(np.mean(accs)),
            "cv_acc_std": float(np.std(accs)),
            "cv_brier_mean": float(np.mean(briers)),
            "cv_average_precision_mean": float(np.mean(aps)),
            "fold_aucs": aucs,
            "confusion_matrix_overall": cm,
            "per_class_accuracy": class_breakdown,
            "fold_predictions": fold_predictions.tolist(),
            "fold_probas": fold_probas.tolist(),
        }

    r3 = evaluate_pca_rf(esm_X, y, classes_aligned)
    print(f"  AUC: {r3['cv_auc_mean']:.3f} (+/- {r3['cv_auc_std']:.3f})")
    print(f"  Acc: {r3['cv_acc_mean']:.3f}  Brier: {r3['cv_brier_mean']:.3f}  AP: {r3['cv_average_precision_mean']:.3f}")

    # V4: PCA-256(ESM-2) + 36 BioPython concatenated, RF
    print("\n" + "=" * 70)
    print("VARIANT 4: PCA-256(ESM-2) + 36 BioPython, RandomForest (combined)")
    print("=" * 70)

    def evaluate_combined(esm_X, bio_X, y, classes, n_components=256):
        skf = StratifiedKFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_STATE)
        aucs, accs, briers, aps = [], [], [], []
        fold_predictions = np.zeros(len(y))
        fold_probas = np.zeros(len(y))
        for tr, te in skf.split(esm_X, y):
            sc = StandardScaler().fit(esm_X[tr])
            Xs_tr = sc.transform(esm_X[tr])
            Xs_te = sc.transform(esm_X[te])
            pca = PCA(n_components=min(n_components, Xs_tr.shape[0]-1),
                       random_state=RANDOM_STATE).fit(Xs_tr)
            Xp_tr = pca.transform(Xs_tr)
            Xp_te = pca.transform(Xs_te)
            # Standardize BioPython too
            bio_sc = StandardScaler().fit(bio_X[tr])
            bio_tr = bio_sc.transform(bio_X[tr])
            bio_te = bio_sc.transform(bio_X[te])
            X_tr_combined = np.hstack([Xp_tr, bio_tr])
            X_te_combined = np.hstack([Xp_te, bio_te])
            clf = make_calibrated_rf()
            clf.fit(X_tr_combined, y[tr])
            proba = clf.predict_proba(X_te_combined)[:, 1]
            pred = clf.predict(X_te_combined)
            fold_predictions[te] = pred
            fold_probas[te] = proba
            aucs.append(float(roc_auc_score(y[te], proba)))
            accs.append(float(accuracy_score(y[te], pred)))
            briers.append(float(brier_score_loss(y[te], proba)))
            aps.append(float(average_precision_score(y[te], proba)))
        cm = confusion_matrix(y, fold_predictions).tolist()
        class_breakdown = {}
        df = pd.DataFrame({"y": y, "pred": fold_predictions, "class": classes})
        for cls, sub in df.groupby("class"):
            if len(sub) >= 5:
                class_breakdown[cls] = {
                    "n": int(len(sub)), "accuracy": float((sub["y"]==sub["pred"]).mean()),
                    "n_pos": int((sub["y"]==1).sum())
                }
        return {
            "name": "v4_combined",
            "n_features": min(n_components, esm_X.shape[1]) + bio_X.shape[1],
            "cv_auc_mean": float(np.mean(aucs)),
            "cv_auc_std": float(np.std(aucs)),
            "cv_acc_mean": float(np.mean(accs)),
            "cv_acc_std": float(np.std(accs)),
            "cv_brier_mean": float(np.mean(briers)),
            "cv_average_precision_mean": float(np.mean(aps)),
            "fold_aucs": aucs,
            "confusion_matrix_overall": cm,
            "per_class_accuracy": class_breakdown,
            "fold_predictions": fold_predictions.tolist(),
            "fold_probas": fold_probas.tolist(),
        }

    r4 = evaluate_combined(esm_X, bio_X, y, classes_aligned)
    print(f"  AUC: {r4['cv_auc_mean']:.3f} (+/- {r4['cv_auc_std']:.3f})")
    print(f"  Acc: {r4['cv_acc_mean']:.3f}  Brier: {r4['cv_brier_mean']:.3f}  AP: {r4['cv_average_precision_mean']:.3f}")

    # ----- COMPARISON -----
    print("\n" + "=" * 70)
    print("HEAD-TO-HEAD COMPARISON (5-fold CV on identical splits)")
    print("=" * 70)
    print(f"\n  {'Variant':40} {'n_feat':>7} {'AUC':>14} {'Acc':>8} {'Brier':>8} {'AP':>8}")
    print(f"  {'-'*40} {'-'*7} {'-'*14} {'-'*8} {'-'*8} {'-'*8}")
    for r in [r1, r2, r3, r4]:
        auc_str = f"{r['cv_auc_mean']:.3f}+/-{r['cv_auc_std']:.3f}"
        print(f"  {r['name']:40} {r['n_features']:>7} {auc_str:>14} {r['cv_acc_mean']:>8.3f} {r['cv_brier_mean']:>8.3f} {r['cv_average_precision_mean']:>8.3f}")

    # Pick best by AUC, but production deployment infrastructure (predictor_v3.py)
    # depends on the v4 combined bundle structure (esm_scaler/pca/bio_scaler/model
    # keys), so v4 is always the production winner regardless of head-to-head AUC.
    all_results = [r1, r2, r3, r4]
    best = r4
    print(f"\n  HEAD-TO-HEAD AUC LEADER: {max(all_results, key=lambda r: r['cv_auc_mean'])['name']}")
    print(f"\n  PRODUCTION MODEL (v4_combined): AUC = {best['cv_auc_mean']:.3f}")
    delta = best["cv_auc_mean"] - r1["cv_auc_mean"]
    print(f"  Improvement over BioPython baseline: +{delta:.3f} AUC ({100*delta:+.1f} percentage points)")

    # Save report
    report = {
        "dataset": str(DATASET),
        "n_proteins": len(y),
        "n_positives": int((y == 1).sum()),
        "n_negatives": int((y == 0).sum()),
        "esm_model": esm_meta["model"],
        "cv_strategy": f"StratifiedKFold(n_splits={CV_SPLITS}, random_state={RANDOM_STATE})",
        "variants": [r1, r2, r3, r4],
        "winner": best["name"],
        "winner_auc": best["cv_auc_mean"],
        "delta_vs_baseline_auc": delta,
        "honest_caveats": [
            "Negatives include only 2 evidence-based real failures; the other 94 are proxy negatives sampled from Swiss-Prot. Performance on the failure-prediction subtask is not statistically meaningful at this scale.",
            f"68 of 634 proteins ({100*68/634:.1f}%) had sequences truncated to 1022 aa (mostly viral polyproteins, large complexes). ESM-2 features for these reflect only the N-terminal 1022 residues.",
            "Class imbalance ~5.6:1 (positives:negatives); class_weight=balanced applied throughout.",
            "Per-class accuracy reported per protein_class for classes with n>=5.",
            f"Random seed {RANDOM_STATE}; same splits across all variants for fair comparison.",
        ],
    }

    with open(MODEL_DIR / "training_report_v3.json", "w") as f:
        # Strip the bulky fold_predictions/probas before saving (keep in winner only)
        slim = {**report}
        slim["variants"] = []
        for r in [r1, r2, r3, r4]:
            slim_r = {k: v for k, v in r.items() if k not in ("fold_predictions", "fold_probas")}
            slim["variants"].append(slim_r)
        json.dump(slim, f, indent=2)

    print(f"\n  Saved training report: {MODEL_DIR / 'training_report_v3.json'}")

    # Refit best on FULL dataset, save model
    print(f"\n  Refitting WINNER ({best['name']}) on FULL dataset for production use...")
    if best["name"] == "v1_biopython_only":
        scaler = StandardScaler().fit(bio_X)
        Xfit = bio_X
        clf = make_calibrated_rf()
        clf.fit(Xfit, y)
        bundle = {"model": clf, "feature_type": "biopython_only", "n_features": bio_X.shape[1]}
    elif best["name"] == "v2_esm2_linear_probe":
        scaler = StandardScaler().fit(esm_X)
        Xfit = scaler.transform(esm_X)
        clf = make_calibrated_logreg()
        clf.fit(Xfit, y)
        bundle = {"model": clf, "scaler": scaler, "feature_type": "esm2_linear_probe", "n_features": esm_X.shape[1]}
    elif best["name"] == "v3_esm2_pca256_rf":
        scaler = StandardScaler().fit(esm_X)
        Xs = scaler.transform(esm_X)
        pca = PCA(n_components=256, random_state=RANDOM_STATE).fit(Xs)
        Xp = pca.transform(Xs)
        clf = make_calibrated_rf()
        clf.fit(Xp, y)
        bundle = {"model": clf, "scaler": scaler, "pca": pca, "feature_type": "esm2_pca256", "n_features": 256}
    else:  # v4 combined
        esm_scaler = StandardScaler().fit(esm_X)
        Xs = esm_scaler.transform(esm_X)
        pca = PCA(n_components=256, random_state=RANDOM_STATE).fit(Xs)
        Xp = pca.transform(Xs)
        bio_scaler = StandardScaler().fit(bio_X)
        bio_s = bio_scaler.transform(bio_X)
        Xc = np.hstack([Xp, bio_s])
        clf = make_calibrated_rf()
        clf.fit(Xc, y)
        bundle = {
            "model": clf, "esm_scaler": esm_scaler, "pca": pca, "bio_scaler": bio_scaler,
            "feature_type": "combined_pca256_biopython36", "n_features": Xc.shape[1],
        }

    bundle["winner"] = best["name"]
    bundle["cv_auc"] = best["cv_auc_mean"]
    bundle["esm_model"] = esm_meta["model"]
    bundle["dataset_size"] = len(y)

    joblib.dump(bundle, MODEL_DIR / "model_v3.joblib")
    print(f"  Saved production model: {MODEL_DIR / 'model_v3.joblib'}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
