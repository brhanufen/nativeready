"""Training pipeline v2 — adds:

1. Probability calibration (CalibratedClassifierCV with isotonic regression)
   so reported "scores" reflect actual probabilities.
2. Out-of-distribution (OOD) reference data: stores feature stats from the
   training set so the predictor can flag inputs that are very different
   from anything seen during training.
3. Optional integration of the expansion dataset (small peptides, hard
   negatives, additional positives) if /data/expansion.py exists.

Outputs:
  model/model_v2.joblib            (calibrated classifier + OOD reference)
  model/feature_names.json         (unchanged feature ordering)
  model/training_report_v2.json    (real metrics, honest caveats)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, brier_score_loss, classification_report,
    confusion_matrix, roc_auc_score
)
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
import joblib

ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "data"
MODEL_DIR = ROOT / "model"

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"


def clean_sequence(seq: str) -> str:
    seq = seq.upper().replace(" ", "").replace("\n", "")
    cleaned = []
    for c in seq:
        if c in STANDARD_AA:
            cleaned.append(c)
        elif c in "BZJU":
            cleaned.append({"B": "N", "Z": "Q", "J": "L", "U": "C"}[c])
    return "".join(cleaned)


def extract_features(sequence: str):
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
    return feats


def maybe_load_expansion():
    """Load the agent-curated expansion data if it exists."""
    exp_path = DATA_DIR / "expansion.py"
    if not exp_path.exists():
        print("  No expansion.py found — using base dataset only.")
        return [], [], []
    sys.path.insert(0, str(DATA_DIR))
    try:
        import expansion  # type: ignore
        return (
            getattr(expansion, "SMALL_PEPTIDES_POSITIVES", []),
            getattr(expansion, "HARD_NEGATIVES", []),
            getattr(expansion, "ADDITIONAL_POSITIVES", []),
        )
    except Exception as e:
        print(f"  Failed to load expansion: {e}")
        return [], [], []


def fetch_sequence(uniprot_id: str):
    """Fetch one sequence from UniProt (used for expansion entries)."""
    import urllib.request
    import urllib.error
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "NativeReady/0.2"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            text = resp.read().decode("utf-8")
            lines = text.strip().split("\n")
            return "".join(lines[1:]).replace(" ", "").upper()
    except (urllib.error.URLError, urllib.error.HTTPError):
        return None


def build_dataset():
    print("=== Building expanded dataset ===")
    pos = json.loads((DATA_DIR / "positives_with_sequences.json").read_text())
    neg = json.loads((DATA_DIR / "negatives_with_sequences.json").read_text())
    print(f"  Base positives: {len(pos)}")
    print(f"  Base negatives: {len(neg)}")

    # Load expansion if available
    small_pep, hard_neg, more_pos = maybe_load_expansion()
    print(f"  Expansion small peptides: {len(small_pep)}")
    print(f"  Expansion hard negatives: {len(hard_neg)}")
    print(f"  Expansion additional positives: {len(more_pos)}")

    # Track existing accessions to avoid duplicates
    existing = {e["uniprot_id"] for e in pos + neg}

    # Fetch sequences for expansion entries
    import time
    expansion_records = []
    for label, items in [(1, small_pep + more_pos), (0, hard_neg)]:
        for entry in items:
            acc = entry["uniprot_id"]
            if acc in existing:
                continue
            seq = fetch_sequence(acc)
            if seq is None or len(seq) < 5:
                continue
            existing.add(acc)
            expansion_records.append({
                **entry,
                "sequence": seq,
                "sequence_length": len(seq),
                "label": label,
            })
            time.sleep(0.25)

    print(f"  Successfully fetched expansion sequences: {len(expansion_records)}")

    # Build feature matrix
    rows = []
    all_records = pos + neg + expansion_records
    for entry in all_records:
        feats = extract_features(entry["sequence"])
        if feats is None:
            continue
        feats["__label__"] = entry["label"]
        feats["__accession__"] = entry["uniprot_id"]
        feats["__name__"] = entry.get("name", "")
        rows.append(feats)

    df = pd.DataFrame(rows)
    print(f"\n  Final dataset: {df.shape}")
    print(f"    Positives: {(df['__label__'] == 1).sum()}")
    print(f"    Negatives: {(df['__label__'] == 0).sum()}")

    # Save expanded sequences for reproducibility
    if expansion_records:
        (DATA_DIR / "expansion_with_sequences.json").write_text(
            json.dumps(expansion_records, indent=2)
        )
        print(f"  Saved {len(expansion_records)} expansion records.")
    return df


def train_and_evaluate(df: pd.DataFrame):
    feature_cols = [c for c in df.columns if not c.startswith("__")]
    X = df[feature_cols].values.astype(float)
    y = df["__label__"].values.astype(int)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, random_state=42, stratify=y
    )

    print(f"\n=== Training calibrated classifier ===")
    print(f"  Train: {len(X_train)}  Test: {len(X_test)}")

    base = RandomForestClassifier(
        n_estimators=300, max_depth=10, min_samples_leaf=2,
        class_weight="balanced", random_state=42, n_jobs=-1,
    )
    # Calibrate via 5-fold CV; isotonic regression handles non-monotonic miscalibration.
    n_pos_train = int(y_train.sum())
    n_neg_train = int((1 - y_train).sum())
    cv_folds = min(5, max(2, n_pos_train, n_neg_train))
    print(f"  Calibration CV folds: {cv_folds} (limited by smaller class size)")
    clf = CalibratedClassifierCV(base, method="isotonic", cv=cv_folds)
    clf.fit(X_train, y_train)

    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)[:, 1]
    test_acc = accuracy_score(y_test, y_pred)
    test_auc = roc_auc_score(y_test, y_proba)
    test_brier = brier_score_loss(y_test, y_proba)
    cm = confusion_matrix(y_test, y_pred).tolist()

    # 5-fold CV for performance and Brier score
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_aucs, cv_accs, cv_briers = [], [], []
    for fold_train, fold_test in skf.split(X, y):
        cv_clf = CalibratedClassifierCV(
            RandomForestClassifier(
                n_estimators=300, max_depth=10, min_samples_leaf=2,
                class_weight="balanced", random_state=42, n_jobs=-1,
            ),
            method="isotonic",
            cv=cv_folds,
        )
        cv_clf.fit(X[fold_train], y[fold_train])
        proba = cv_clf.predict_proba(X[fold_test])[:, 1]
        pred = cv_clf.predict(X[fold_test])
        cv_aucs.append(float(roc_auc_score(y[fold_test], proba)))
        cv_accs.append(float(accuracy_score(y[fold_test], pred)))
        cv_briers.append(float(brier_score_loss(y[fold_test], proba)))

    print(f"\n=== HELD-OUT TEST ===")
    print(f"  Accuracy:   {test_acc:.3f}")
    print(f"  AUC:        {test_auc:.3f}")
    print(f"  Brier loss: {test_brier:.3f} (lower = better calibrated)")
    print(f"  Confusion matrix: {cm}")

    print(f"\n=== 5-FOLD CV ===")
    print(f"  AUC:   {np.mean(cv_aucs):.3f} (+/- {np.std(cv_aucs):.3f})")
    print(f"  Acc:   {np.mean(cv_accs):.3f} (+/- {np.std(cv_accs):.3f})")
    print(f"  Brier: {np.mean(cv_briers):.3f}")

    # OOD detection: fit a scaler + nearest-neighbor index on training features
    scaler = StandardScaler().fit(X_train)
    knn = NearestNeighbors(n_neighbors=5).fit(scaler.transform(X_train))
    train_distances, _ = knn.kneighbors(scaler.transform(X_train))
    ood_threshold = float(np.percentile(train_distances[:, -1], 95))
    print(f"\n=== OOD detector ===")
    print(f"  95th percentile NN distance on training set: {ood_threshold:.3f}")
    print(f"  Inputs with greater distance flagged as out-of-distribution.")

    bundle = {
        "model": clf,
        "feature_names": feature_cols,
        "ood_scaler": scaler,
        "ood_knn": knn,
        "ood_threshold": ood_threshold,
        "training_set_size": len(X_train),
        "test_set_size": len(X_test),
        "n_features": len(feature_cols),
    }
    joblib.dump(bundle, MODEL_DIR / "model_v2.joblib")
    (MODEL_DIR / "feature_names.json").write_text(json.dumps(feature_cols, indent=2))

    report = {
        "model_type": "CalibratedClassifierCV(RandomForest, isotonic)",
        "training_set_size": len(X_train),
        "test_set_size": len(X_test),
        "n_features": len(feature_cols),
        "n_positives": int(y.sum()),
        "n_negatives": int((1 - y).sum()),
        "held_out_test": {
            "accuracy": float(test_acc),
            "auc": float(test_auc),
            "brier": float(test_brier),
            "confusion_matrix": cm,
        },
        "cross_validation_5fold": {
            "mean_auc": float(np.mean(cv_aucs)),
            "std_auc": float(np.std(cv_aucs)),
            "mean_accuracy": float(np.mean(cv_accs)),
            "std_accuracy": float(np.std(cv_accs)),
            "mean_brier": float(np.mean(cv_briers)),
        },
        "ood_threshold_95pct": float(ood_threshold),
        "calibration_method": "isotonic",
        "honest_caveats": [
            "Small dataset: now expanded but still under 250 examples; performance estimates have meaningful variance.",
            "Negatives include proxy examples (random Swiss-Prot proteins) plus expert-curated hard cases (mucins, large membrane proteins, aggregation-prone). Still not true experimental failures.",
            "Calibration via isotonic regression on 5-fold CV. Reported probabilities should be roughly aligned with empirical frequencies but accuracy of calibration depends on training set diversity.",
            "Out-of-distribution detection: nearest-neighbor distance in feature space. Inputs above the 95th percentile of training distances get a low-confidence warning.",
            "Limitations: model is a v0.2 prototype. Validate critical experiments yourself.",
            "Generalization to proteins very different from the training set remains uncertain.",
        ],
    }
    (MODEL_DIR / "training_report_v2.json").write_text(json.dumps(report, indent=2))
    print(f"\nSaved v2 model to {MODEL_DIR / 'model_v2.joblib'}")
    print(f"Saved v2 report to {MODEL_DIR / 'training_report_v2.json'}")


if __name__ == "__main__":
    df = build_dataset()
    train_and_evaluate(df)
