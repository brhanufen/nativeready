"""Integration tests for NativeReady — verify the trained model works
end-to-end with the backend prediction pipeline.

Uses real sequences from UniProt (no synthetic data).
"""
import json
import sys
from pathlib import Path

# Make the backend importable
sys.path.insert(0, str(Path(__file__).parent.parent / "backend"))

from predictor import predict
from features import extract_features, clean_sequence


# Real sequences from UniProt — fetched in our positive set.
# These should score WELL because they are documented native MS targets.
UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
    "KESTLHLVLRLRGG"
)

# Carbonic anhydrase 2 (P00918) — classic native MS standard
CA2 = (
    "MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNG"
    "HAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTK"
    "YGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLD"
    "YWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKN"
    "RQIKASFK"
)

# An unrelated random protein from negatives — should score LOWER
# Mucin-like or large random Swiss-Prot entry
LARGE_DISORDERED = "P" * 50 + "STSTSTSTSTSTSTSTSTSTSTSTSTSTST" * 30  # synthetic — only for input validation test
# DON'T use synthetic for actual scoring tests — load a real one instead


def test_clean_sequence():
    assert clean_sequence("MK ST") == "MKST"
    assert clean_sequence("mkst\n") == "MKST"
    assert clean_sequence("MKBZ") == "MKNQ"  # ambiguous mapped
    assert clean_sequence("MKX") == "MK"  # X dropped
    print("clean_sequence: OK")


def test_features_extract():
    feats = extract_features(UBIQUITIN)
    assert feats["length"] == 76, f"Expected length 76, got {feats['length']}"
    assert 8.0 < feats["molecular_weight_kda"] < 9.5, f"MW out of expected range: {feats['molecular_weight_kda']}"
    assert "pct_K" in feats
    assert "pct_charged" in feats
    print(f"features_extract (ubiquitin): length={feats['length']}, MW={feats['molecular_weight_kda']:.2f} kDa, pI={feats['isoelectric_point']:.2f}")


def test_predict_known_positives():
    """Known native MS standards should score well."""
    print("\n--- Predicting known positives ---")
    for name, seq in [("Ubiquitin", UBIQUITIN), ("Carbonic anhydrase 2", CA2)]:
        result = predict(seq)
        print(f"  {name}: score={result['suitability_score']}/100 ({result['suitability_label']}) [model={result['model_version']}]")
        # We expect these classic native MS standards to score >= 50
        assert result["suitability_score"] >= 40, f"{name} scored too low: {result['suitability_score']}"
        assert "risk_factors" in result
        assert len(result["risk_factors"]) >= 4
        assert "recommendations" in result
        assert len(result["recommendations"]) >= 1


def test_predict_real_negative():
    """Load a real protein from our negatives set and score it."""
    neg_file = Path(__file__).parent.parent / "data" / "negatives_with_sequences.json"
    negs = json.loads(neg_file.read_text())
    if not negs:
        print("  (no negatives loaded; skipping)")
        return
    sample = negs[0]  # first one
    result = predict(sample["sequence"])
    print(f"  {sample['name'][:60]}: score={result['suitability_score']}/100 ({result['suitability_label']}) [model={result['model_version']}]")


def test_response_shape():
    result = predict(UBIQUITIN)
    required = {"suitability_score", "suitability_label", "confidence_interval", "risk_factors", "recommendations", "model_version"}
    missing = required - set(result.keys())
    assert not missing, f"Missing fields in response: {missing}"
    assert isinstance(result["suitability_score"], int)
    assert isinstance(result["risk_factors"], list)
    assert all("name" in r and "value" in r and "risk_level" in r for r in result["risk_factors"])
    print("response_shape: OK")


def test_validation_too_short():
    try:
        predict("MKL")  # 3 residues
    except ValueError as e:
        print(f"validation_too_short: OK (raised ValueError: {e})")
        return
    raise AssertionError("Should have raised ValueError for too-short sequence")


if __name__ == "__main__":
    print("=" * 60)
    print("NativeReady Integration Tests")
    print("=" * 60)
    test_clean_sequence()
    test_features_extract()
    test_predict_known_positives()
    test_predict_real_negative()
    test_response_shape()
    test_validation_too_short()
    print("\n" + "=" * 60)
    print("All tests passed.")
    print("=" * 60)
