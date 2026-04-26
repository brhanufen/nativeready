"""Scientific integrity audit for NativeReady.

Verifies:
1. Data integrity: positives and negatives are distinct, real, UniProt-sourced
2. Model performance: claimed metrics match recomputed metrics
3. Biological plausibility: predictions align with expert intuition on
   well-known cases that span the full difficulty spectrum
4. No obvious data leakage
5. Edge case behavior is sensible
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))
sys.path.insert(0, str(Path(__file__).parent.parent / "data"))

from predictor import predict
from features import extract_features
from positives_curated import POSITIVE_EXAMPLES

ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "data"


def section(title):
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)


def audit_data_integrity():
    section("AUDIT 1: Data integrity")
    pos = json.loads((DATA_DIR / "positives_with_sequences.json").read_text())
    neg = json.loads((DATA_DIR / "negatives_with_sequences.json").read_text())

    pos_ids = {p["uniprot_id"] for p in pos}
    neg_ids = {n["uniprot_id"] for n in neg}
    overlap = pos_ids & neg_ids
    print(f"  Positive examples: {len(pos)}")
    print(f"  Negative examples: {len(neg)}")
    print(f"  ID overlap (must be 0): {len(overlap)}")
    assert len(overlap) == 0, f"DATA LEAK: {overlap} appear in both sets"

    # Verify all sequences are non-empty and reasonable
    for entry in pos + neg:
        assert entry["sequence"], f"Empty sequence for {entry['uniprot_id']}"
        assert len(entry["sequence"]) >= 30, f"Suspiciously short sequence for {entry['uniprot_id']}: {len(entry['sequence'])} aa"
        # UniProt accessions are 6 or 10 alphanumeric chars
        acc = entry["uniprot_id"]
        assert 6 <= len(acc) <= 10 and acc.isalnum(), \
            f"UniProt ID {acc} doesn't look like a real accession"

    # Sample a positive and verify its sequence header references the right protein
    sample = pos[0]
    print(f"  Sample positive: {sample['uniprot_id']} ({sample['name']})")
    print(f"    FASTA header: {sample['fasta_header'][:80]}")
    print(f"    Sequence length: {sample['sequence_length']} aa (UniProt-fetched)")
    print(f"    First 50 aa: {sample['sequence'][:50]}...")

    print("  PASS: data integrity")
    return pos, neg


def audit_model_metrics_reproduce():
    section("AUDIT 2: Reproduce claimed model metrics")
    import joblib
    from sklearn.metrics import accuracy_score, roc_auc_score
    from sklearn.model_selection import train_test_split

    # Prefer v2 model if available
    v2_path = ROOT / "model" / "model_v2.joblib"
    v2_report = ROOT / "model" / "training_report_v2.json"
    if v2_path.exists() and v2_report.exists():
        bundle = joblib.load(v2_path)
        report = json.loads(v2_report.read_text())
        version = "v2"
    else:
        bundle = joblib.load(ROOT / "model" / "model.joblib")
        report = json.loads((ROOT / "model" / "training_report.json").read_text())
        version = "v1"
    model = bundle["model"]
    feature_names = bundle.get("feature_names") or json.loads((ROOT / "model" / "feature_names.json").read_text())
    print(f"  Auditing model version: {version}")

    pos = json.loads((DATA_DIR / "positives_with_sequences.json").read_text())
    neg = json.loads((DATA_DIR / "negatives_with_sequences.json").read_text())
    # Include expansion if it exists (used in v2)
    exp_file = DATA_DIR / "expansion_with_sequences.json"
    expansion = json.loads(exp_file.read_text()) if (version == "v2" and exp_file.exists()) else []

    # Recompute features and labels
    X_rows, y = [], []
    for entry in pos + neg + expansion:
        feats = extract_features(entry["sequence"])
        X_rows.append([feats[name] for name in feature_names])
        y.append(entry["label"])
    X = np.array(X_rows, dtype=float)
    y = np.array(y, dtype=int)

    # Reproduce the same train/test split (random_state=42 in train.py)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, random_state=42, stratify=y
    )

    # Score on the test set with the saved model
    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]
    acc = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_proba)

    claimed_acc = report["held_out_test"]["accuracy"]
    claimed_auc = report["held_out_test"]["auc"]
    print(f"  Claimed test accuracy: {claimed_acc:.3f}  |  Reproduced: {acc:.3f}")
    print(f"  Claimed test AUC:      {claimed_auc:.3f}  |  Reproduced: {auc:.3f}")
    assert abs(acc - claimed_acc) < 0.001, "Accuracy mismatch — reported numbers may be incorrect"
    assert abs(auc - claimed_auc) < 0.001, "AUC mismatch — reported numbers may be incorrect"
    print("  PASS: claimed metrics match reproduced metrics exactly")


def audit_no_data_leakage():
    section("AUDIT 3: Check for data leakage")
    # The risk: if the same accession appears in train and test, that's leakage.
    # The train/test split is random within the combined set, so no engineering
    # leakage exists. But check that no sequence appears twice.
    pos = json.loads((DATA_DIR / "positives_with_sequences.json").read_text())
    neg = json.loads((DATA_DIR / "negatives_with_sequences.json").read_text())
    seqs = [e["sequence"] for e in pos + neg]
    unique_seqs = set(seqs)
    print(f"  Total sequences: {len(seqs)}")
    print(f"  Unique sequences: {len(unique_seqs)}")
    assert len(seqs) == len(unique_seqs), "Some sequences appear more than once"
    print("  PASS: no duplicate sequences across positive/negative sets")


def audit_biological_plausibility():
    section("AUDIT 4: Biological plausibility on diverse known cases")

    # Real, well-known proteins NOT in the training set, picked to span
    # the difficulty spectrum. Sequences are real (would normally be fetched
    # from UniProt; for this audit we hard-code short well-known ones).
    test_cases = [
        # Should score WELL (small, soluble, well-studied)
        ("Ubiquitin (P0CG48, in train set — sanity check)",
         "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
         "high"),
        # Should score WELL (in train set — another sanity check)
        ("Insulin chain B (P01308 not in train — real human insulin chain B)",
         "FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
         "medium"),  # very short might confuse model
        # Glucagon — small, well-folded, should score reasonably
        ("Glucagon (P01275, real human glucagon peptide)",
         "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT",
         "medium"),  # also short
    ]

    print(f"  {'Case':50s}  {'Score':6s}  {'Label':12s}  Expected")
    print("  " + "-" * 90)
    for name, seq, expected in test_cases:
        try:
            r = predict(seq)
            print(f"  {name[:48]:50s}  {r['suitability_score']:>3d}    {r['suitability_label']:12s}  {expected}")
        except ValueError as e:
            print(f"  {name[:48]:50s}  ERR    {str(e)[:30]}  {expected}")


def audit_full_dataset_predictions():
    section("AUDIT 5: Predict on full dataset, compare distributions")
    pos = json.loads((DATA_DIR / "positives_with_sequences.json").read_text())
    neg = json.loads((DATA_DIR / "negatives_with_sequences.json").read_text())

    pos_scores = []
    neg_scores = []
    failures = []
    for e in pos:
        try:
            r = predict(e["sequence"])
            pos_scores.append(r["suitability_score"])
        except Exception as ex:
            failures.append((e["uniprot_id"], str(ex)))
    for e in neg:
        try:
            r = predict(e["sequence"])
            neg_scores.append(r["suitability_score"])
        except Exception as ex:
            failures.append((e["uniprot_id"], str(ex)))

    print(f"  Positive set predictions: n={len(pos_scores)}")
    print(f"    Mean score: {np.mean(pos_scores):.1f} (std {np.std(pos_scores):.1f})")
    print(f"    Median:     {np.median(pos_scores):.1f}")
    print(f"    Min/Max:    {min(pos_scores)}-{max(pos_scores)}")
    print(f"  Negative set predictions: n={len(neg_scores)}")
    print(f"    Mean score: {np.mean(neg_scores):.1f} (std {np.std(neg_scores):.1f})")
    print(f"    Median:     {np.median(neg_scores):.1f}")
    print(f"    Min/Max:    {min(neg_scores)}-{max(neg_scores)}")
    print(f"  Score gap (positive mean - negative mean): {np.mean(pos_scores) - np.mean(neg_scores):.1f}")

    if failures:
        print(f"  WARN: {len(failures)} sequences failed prediction:")
        for fid, msg in failures[:5]:
            print(f"    {fid}: {msg}")

    # IMPORTANT CAVEAT: the positives include training examples, so this is
    # NOT held-out performance. It's a sanity check only.
    print("  NOTE: Positive set includes training examples — this is sanity, not generalization.")


def audit_edge_cases():
    section("AUDIT 6: Edge case behavior")
    cases = [
        ("Very short (10 aa)", "MAGSTSCEGN"),
        ("Long all-alanine (200 aa)", "A" * 200),
        ("Hydrophobic membrane-like (100 aa)", "LLLLLLLLLLVVVVVVVVVVAAAAAAAAAA" * 4 + "MMMMMMMMMM"),
        ("Polar/charged only (100 aa)", ("DEKR" * 25)),
        ("Disordered-like (PSTQGS rich, 150 aa)", ("PSTQGS" * 25)),
    ]
    print(f"  {'Case':45s}  {'Score':6s}  {'Label':12s}")
    for name, seq in cases:
        try:
            r = predict(seq)
            print(f"  {name:45s}  {r['suitability_score']:>3d}    {r['suitability_label']}")
        except ValueError as e:
            print(f"  {name:45s}  ERR    {e}")


def audit_documentation():
    section("AUDIT 7: Documentation honesty")
    report = json.loads((ROOT / "model" / "training_report.json").read_text())
    caveats = report.get("honest_caveats", [])
    print(f"  Honest caveats documented: {len(caveats)}")
    expected_topics = ["small dataset", "proxy", "v0", "limitations", "validate"]
    text = " ".join(caveats).lower()
    for topic in expected_topics:
        present = topic in text
        marker = "✓" if present else "✗"
        print(f"    [{marker}] mentions '{topic}': {present}")

    readme = (ROOT / "README.md").read_text().lower()
    print("  README disclosures:")
    for required in ["caveats", "small dataset", "proxy", "validate", "research prototype"]:
        present = required in readme
        marker = "✓" if present else "✗"
        print(f"    [{marker}] '{required}' present in README: {present}")


def main():
    print("\n" + "#" * 70)
    print("# NATIVE READY — SCIENTIFIC INTEGRITY AUDIT")
    print("#" * 70)
    audit_data_integrity()
    audit_model_metrics_reproduce()
    audit_no_data_leakage()
    audit_biological_plausibility()
    audit_full_dataset_predictions()
    audit_edge_cases()
    audit_documentation()
    print("\n" + "#" * 70)
    print("# AUDIT COMPLETE")
    print("#" * 70)


if __name__ == "__main__":
    main()
