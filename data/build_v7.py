"""Build dataset_combined_v7_2026-05-11.json from v6, removing CB2 (P34972).

CB2 was previously included in v6 as a confirmed real failure based on a
personal communication. Subsequent clarification indicated the experimental
outcome is undetermined (the protein has not yet succeeded in native MS
but may eventually). The honest position is to remove it from training
and treat it as a held-out test case until the experiment converges.

After removal: 635 records (538 positives + 97 negatives, of which 3
remain evidence-based real failures: insulin P01317, AAV8 VP1 Q8JQF8,
thyroglobulin P01267).
"""
import json
from pathlib import Path

DATA = Path("/Users/bfentaw2/startup/nativeready/data")
SRC = DATA / "dataset_combined_v6_2026-05-11.json"
DST = DATA / "dataset_combined_v7_2026-05-11.json"

with open(SRC) as f:
    records = json.load(f)
print(f"Loaded {len(records)} records from {SRC.name}")

before = len(records)
records = [r for r in records if r["uniprot_id"] != "P34972"]
removed = before - len(records)
assert removed == 1, f"Expected to remove exactly 1 record (CB2 P34972), removed {removed}"
print(f"Removed CB2 P34972 -> {len(records)} records")

pos = sum(1 for r in records if r["label"] == 1)
neg = sum(1 for r in records if r["label"] == 0)
real_fails = [r for r in records if r.get("label") == 0 and r.get("failure_mode")]
print(f"  Positives: {pos}")
print(f"  Negatives: {neg}")
print(f"  Real failures (with failure_mode): {len(real_fails)}")
for r in real_fails:
    print(f"    - {r['uniprot_id']:8s} {r['name'][:50]:50s} mode={r['failure_mode']}")

with open(DST, "w") as f:
    json.dump(records, f, indent=2)
print(f"\nSaved: {DST}")
