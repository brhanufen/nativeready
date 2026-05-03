"""Basic usage examples for the NativeReady Python SDK."""
from nativeready import predict, Client


# ==========================================================================
# Example 1: one-liner prediction
# ==========================================================================
print("=" * 60)
print("Example 1: one-liner")
print("=" * 60)

ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
result = predict(ubiquitin)
print(result)
print(f"Score: {result.score}/100 ({result.label})")
print(f"95% CI: {result.confidence_lower}-{result.confidence_upper}")
print()


# ==========================================================================
# Example 2: from a UniProt accession
# ==========================================================================
print("=" * 60)
print("Example 2: from UniProt accession")
print("=" * 60)

client = Client()
ca2 = client.predict_uniprot("P00918")  # Carbonic anhydrase 2
print(f"{ca2.uniprot_id}: {ca2.score}/100 ({ca2.label}), len={ca2.sequence_length}")
print()


# ==========================================================================
# Example 3: batch prediction
# ==========================================================================
print("=" * 60)
print("Example 3: batch prediction (3 well-known native MS standards)")
print("=" * 60)

batch = [
    {"id": "ubiquitin", "sequence": ubiquitin},
    {"id": "lysozyme", "sequence":
     "MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQAT"
     "NRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNG"
     "MNAWVAWRNRCKGTDVQAWIRGCRL"},
    {"id": "cytochrome_c", "sequence":
     "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITW"
     "KEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"},
]

results = client.predict_batch(batch, progress=False)
for r in results:
    print(f"  {r.uniprot_id or '?':16} score={r.score:3d}  {r.label}")
