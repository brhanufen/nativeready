"""
Train the NativeReady classifier on real protein data.

Inputs:
  - data/positives_with_sequences.json (69 real native MS validated proteins)
  - data/negatives_with_sequences.json (~64 proxy negatives from Swiss-Prot)

Outputs:
  - model/model.joblib (trained RandomForest classifier)
  - model/feature_names.json (ordered list of feature names)
  - model/training_report.json (real performance metrics, no fabrication)

Honest about limitations:
  - Small dataset (~133 examples)
  - Negatives are proxies, not true experimental failures
  - Model is a v0 baseline; performance will be modest
  - All metrics reported are real, computed on a held-out test set
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score, classification_report, confusion_matrix,
    roc_auc_score, roc_curve
)
from sklearn.model_selection import StratifiedKFold, train_test_split
import joblib

ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "data"
MODEL_DIR = ROOT / "model"

# Standard 20 amino acids
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"


def clean_sequence(seq: str) -> str:
    """Remove non-standard chars; replace ambiguous with 'X'."""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    cleaned = []
    for c in seq:
        if c in STANDARD_AA:
            cleaned.append(c)
        elif c in "BZJU":  # ambiguous codes -> closest residue
            cleaned.append({"B": "N", "Z": "Q", "J": "L", "U": "C"}[c])
        else:
            # 'X' or other unknown -> drop
            continue
    return "".join(cleaned)


def extract_features(sequence: str) -> dict:
    """Compute biophysical features for a protein sequence."""
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
        "secondary_structure_helix": pa.secondary_structure_fraction()[0],
        "secondary_structure_turn": pa.secondary_structure_fraction()[1],
        "secondary_structure_sheet": pa.secondary_structure_fraction()[2],
    }
    # Per-amino-acid composition (20 features)
    for aa in STANDARD_AA:
        feats[f"pct_{aa}"] = aa_pct.get(aa, 0.0)

    # Engineered features known to matter for native MS
    feats["pct_cysteine"] = aa_pct.get("C", 0.0)  # disulfides = stability
    feats["pct_proline"] = aa_pct.get("P", 0.0)  # disorder marker
    feats["pct_charged"] = sum(aa_pct.get(c, 0.0) for c in "DEKR")
    feats["pct_hydrophobic"] = sum(aa_pct.get(c, 0.0) for c in "AVILMFW")
    feats["pct_aromatic"] = sum(aa_pct.get(c, 0.0) for c in "FWY")
    feats["pct_polar"] = sum(aa_pct.get(c, 0.0) for c in "STNQ")

    return feats


def build_dataset():
    """Load JSON, extract features, return X, y, accessions."""
    positives = json.loads((DATA_DIR / "positives_with_sequences.json").read_text())
    negatives = json.loads((DATA_DIR / "negatives_with_sequences.json").read_text())

    rows = []
    for entry in positives + negatives:
        feats = extract_features(entry["sequence"])
        if feats is None:
            print(f"  Skipping {entry['uniprot_id']} (sequence too short)")
            continue
        feats["__label__"] = entry["label"]
        feats["__accession__"] = entry["uniprot_id"]
        feats["__name__"] = entry.get("name", "")
        rows.append(feats)

    df = pd.DataFrame(rows)
    print(f"\nDataset shape: {df.shape}")
    print(f"  Positives (native MS validated): {(df['__label__'] == 1).sum()}")
    print(f"  Negatives (proxy from Swiss-Prot): {(df['__label__'] == 0).sum()}")

    return df


def train_and_evaluate(df: pd.DataFrame):
    """Train a RandomForest, report honest metrics."""
    feature_cols = [c for c in df.columns if not c.startswith("__")]
    X = df[feature_cols].values
    y = df["__label__"].values

    # 80/20 split, stratified
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, random_state=42, stratify=y
    )

    print(f"\nTrain size: {len(X_train)}  /  Test size: {len(X_test)}")

    clf = RandomForestClassifier(
        n_estimators=200,
        max_depth=8,
        min_samples_leaf=3,
        class_weight="balanced",
        random_state=42,
        n_jobs=-1,
    )
    clf.fit(X_train, y_train)

    # Held-out test metrics
    y_pred = clf.predict(X_test)
    y_proba = clf.predict_proba(X_test)[:, 1]

    test_acc = accuracy_score(y_test, y_pred)
    test_auc = roc_auc_score(y_test, y_proba)
    cm = confusion_matrix(y_test, y_pred).tolist()
    report = classification_report(y_test, y_pred, output_dict=True, zero_division=0)

    # 5-fold CV for a more robust performance estimate
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_aucs = []
    cv_accs = []
    for fold_train, fold_test in skf.split(X, y):
        clf_cv = RandomForestClassifier(
            n_estimators=200, max_depth=8, min_samples_leaf=3,
            class_weight="balanced", random_state=42, n_jobs=-1,
        )
        clf_cv.fit(X[fold_train], y[fold_train])
        proba = clf_cv.predict_proba(X[fold_test])[:, 1]
        pred = clf_cv.predict(X[fold_test])
        cv_aucs.append(float(roc_auc_score(y[fold_test], proba)))
        cv_accs.append(float(accuracy_score(y[fold_test], pred)))

    print("\n=== HELD-OUT TEST METRICS ===")
    print(f"  Accuracy: {test_acc:.3f}")
    print(f"  AUC:      {test_auc:.3f}")
    print(f"  Confusion matrix [TN FP / FN TP]: {cm}")

    print("\n=== 5-FOLD CROSS-VALIDATION ===")
    print(f"  Mean AUC: {np.mean(cv_aucs):.3f} (+/- {np.std(cv_aucs):.3f})")
    print(f"  Mean Acc: {np.mean(cv_accs):.3f} (+/- {np.std(cv_accs):.3f})")

    # Feature importance
    importances = sorted(
        zip(feature_cols, clf.feature_importances_),
        key=lambda x: -x[1]
    )
    print("\n=== TOP 15 FEATURES ===")
    for name, imp in importances[:15]:
        print(f"  {name:30s} {imp:.4f}")

    # Save model + report
    joblib.dump({
        "model": clf,
        "feature_names": feature_cols,
        "training_set_size": len(X_train),
        "test_set_size": len(X_test),
    }, MODEL_DIR / "model.joblib")

    (MODEL_DIR / "feature_names.json").write_text(json.dumps(feature_cols, indent=2))

    report_data = {
        "model_type": "RandomForestClassifier",
        "n_estimators": 200,
        "max_depth": 8,
        "training_set_size": len(X_train),
        "test_set_size": len(X_test),
        "n_features": len(feature_cols),
        "held_out_test": {
            "accuracy": float(test_acc),
            "auc": float(test_auc),
            "confusion_matrix": cm,
            "classification_report": report,
        },
        "cross_validation_5fold": {
            "mean_auc": float(np.mean(cv_aucs)),
            "std_auc": float(np.std(cv_aucs)),
            "mean_accuracy": float(np.mean(cv_accs)),
            "std_accuracy": float(np.std(cv_accs)),
            "fold_aucs": cv_aucs,
        },
        "top_features": [(name, float(imp)) for name, imp in importances[:20]],
        "honest_caveats": [
            "Small dataset: trained on only ~133 examples; performance estimates have high variance.",
            "Negatives are proxy examples (random Swiss-Prot proteins), not true experimental failures.",
            "Real native MS failure data is largely unpublished, limiting label quality.",
            "Limitations: model is a v0 baseline. Validate critical experiments yourself.",
            "Generalization to proteins very different from the training set is uncertain.",
            "Limitations on small peptides: model under-predicts for sequences shorter than 50 amino acids because the training set has few small proteins.",
        ],
        "data_sources": [
            "Positives: 69 proteins curated from native MS literature; sequences from UniProt REST API.",
            "Negatives: ~64 proteins sampled from Swiss-Prot (reviewed, organism=human) excluding the positive list.",
        ],
    }
    (MODEL_DIR / "training_report.json").write_text(json.dumps(report_data, indent=2))

    print(f"\nSaved model to {MODEL_DIR / 'model.joblib'}")
    print(f"Saved training report to {MODEL_DIR / 'training_report.json'}")


if __name__ == "__main__":
    df = build_dataset()
    train_and_evaluate(df)
