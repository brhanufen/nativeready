"""Smoke tests against the live NativeReady API.

Run with: python -m pytest tests/test_smoke.py -v

These hit the real production endpoint, so they take 30-60 seconds total
(first call triggers an ESM-2 cold start on the server).
"""
import pytest

from nativeready import Client, predict, parse_fasta
from nativeready.client import _clean_sequence, _extract_uniprot_from_header


CARBONIC_ANHYDRASE_2 = (
    "MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVE"
    "FDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPD"
    "GLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVT"
    "WIVLKEPISVSSEQVLK"
)

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHL"
    "VLRLRGG"
)

CARBONIC_ANHYDRASE_2_FASTA = (
    f">sp|P00918|CAH2_HUMAN Carbonic anhydrase 2\n{CARBONIC_ANHYDRASE_2}\n"
)


# ---------------- Pure unit tests (no network) ----------------

class TestSequenceCleaning:
    def test_strip_fasta_header(self):
        assert _clean_sequence(CARBONIC_ANHYDRASE_2_FASTA) == CARBONIC_ANHYDRASE_2

    def test_uppercase_and_strip_whitespace(self):
        raw = "  mqif vktlt\n gktit\tlevepsd ti  "
        assert _clean_sequence(raw) == "MQIFVKTLTGKTITLEVEPSDTI"

    def test_extract_uniprot_from_swiss_prot_header(self):
        assert _extract_uniprot_from_header(CARBONIC_ANHYDRASE_2_FASTA) == "P00918"

    def test_extract_uniprot_returns_none_for_non_fasta(self):
        assert _extract_uniprot_from_header(CARBONIC_ANHYDRASE_2) is None


class TestFastaParser:
    def test_parse_single_record(self):
        records = parse_fasta(CARBONIC_ANHYDRASE_2_FASTA)
        assert len(records) == 1
        assert records[0]["sequence"] == CARBONIC_ANHYDRASE_2
        assert records[0]["uniprot_id"] == "P00918"

    def test_parse_multiple_records(self):
        text = (
            f">sp|P00918|CAH2_HUMAN Carbonic anhydrase 2\n{CARBONIC_ANHYDRASE_2}\n"
            f">sp|P0CG48|UBC_HUMAN Polyubiquitin\n{UBIQUITIN}\n"
        )
        records = parse_fasta(text)
        assert len(records) == 2
        assert records[0]["uniprot_id"] == "P00918"
        assert records[1]["uniprot_id"] == "P0CG48"


# ---------------- Live API tests (network required) ----------------

@pytest.fixture(scope="session")
def client():
    return Client()


class TestLiveAPI:
    def test_health(self, client):
        h = client.health()
        assert h.get("status") == "ok"

    def test_predict_carbonic_anhydrase(self, client):
        """Canonical native MS calibrant: must score very high."""
        result = client.predict(CARBONIC_ANHYDRASE_2)
        assert result.score >= 80, f"Carbonic anhydrase scored {result.score}"
        assert result.label in ("Excellent", "Good")
        assert result.model_version, "model_version should be set"
        assert "0.3" in result.model_version or "0.2" in result.model_version

    def test_predict_ubiquitin(self, client):
        """Tiny well-folded protein: must score high."""
        result = client.predict(UBIQUITIN)
        assert result.score >= 80
        assert result.label in ("Excellent", "Good")

    def test_predict_with_fasta_header(self, client):
        """Should auto-extract UniProt ID from a FASTA header."""
        result = client.predict(CARBONIC_ANHYDRASE_2_FASTA)
        assert result.uniprot_id == "P00918"
        assert result.score >= 80

    def test_predict_uniprot_lookup(self, client):
        """Fetch from UniProt and predict."""
        result = client.predict_uniprot("P0CG48")  # ubiquitin
        assert result.uniprot_id == "P0CG48"
        assert result.score > 0


if __name__ == "__main__":
    # Allow running as `python tests/test_smoke.py` for a quick check
    import sys
    pytest.main([__file__, "-v"] + sys.argv[1:])
