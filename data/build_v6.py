"""Build dataset_combined_v6_2026-05-11.json from v5 + CB2 (P34972).

CB2 = cannabinoid receptor 2 (human), 360 aa, 7-TM membrane protein.
Added as the 4th evidence-based real failure following Marty's email
feedback (UT Austin, UniDec author, personal communication 2026-05).

Marty's report: CB2 fails native MS in his lab — none of the tested
detergents preserve the intact complex. Current model scores it 84
(Good), which is wrong, because GRAVY underweights TM-helix
hydrophobicity for membrane proteins. This record is the seed for
the membrane-protein failure mode.
"""
import json
from pathlib import Path

DATA = Path("/Users/bfentaw2/startup/nativeready/data")
SRC = DATA / "dataset_combined_v5_2026-05-03.json"
DST = DATA / "dataset_combined_v6_2026-05-11.json"

CB2_SEQUENCE = (
    "MEECWVTEIANGSKDGLDSNPMKDYMILSGPQKTAVAVLCTLLGLLSALENVAVLYLILS"
    "SHQLRRKPSYLFIGSLAGADFLASVVFACSFVNFHVFHGVDSKAVFLLKIGSVTMTFTAS"
    "VGSLLLTAIDRYLCLRYPPSYKALLTRGRALVTLGIMWVLSALVSYLPLMGWTCCPRPCS"
    "ELFPLIPNDYLLSWLLFIAFLFSGIIYTYGHVLWKAHQHVASLSGHQDRQVPGMARMRLD"
    "VRLAKTLGLVLAVLLICWFPVLALMAHSLATTLSDQVKKAFAFCSMLCLINSMVNPVIYA"
    "LRSGEIRSSAHHCLAHWKKCVRGLGSEAKEEAPRSSVTETEADGKITPWPDSRDLDLSDC"
)
assert len(CB2_SEQUENCE) == 360, f"CB2 length unexpected: {len(CB2_SEQUENCE)}"

cb2_record = {
    "uniprot_id": "P34972",
    "name": "Cannabinoid receptor 2 (human)",
    "protein_class": "membrane protein (7-TM GPCR)",
    "mw_kda": 39.7,
    "notes": (
        "Source: Michael Marty (UT Austin), personal communication 2026-05. "
        "Failure mode: none of the tested detergents preserve the intact "
        "complex in native MS; GPCR/7-TM membrane proteins remain a hard "
        "case for detergent screening. Surfaced model weakness: current "
        "v0.3 scores P34972 at 84 (Good) because GRAVY is a poor proxy "
        "for membrane-helix hydrophobicity. Seed record for the "
        "detergent-incompatible / membrane failure mode."
    ),
    "sequence": CB2_SEQUENCE,
    "sequence_length": 360,
    "fasta_header": "sp|P34972|CNR2_HUMAN Cannabinoid receptor 2 OS=Homo sapiens OX=9606 GN=CNR2 PE=1 SV=1",
    "label": 0,
    "outcome_label": 5,
    "failure_mode": "detergent_incompatible",
    "evidence_paragraph": (
        "Personal communication from Michael Marty (UT Austin), 2026-05: "
        "CB2 has been tested in native MS but no detergent screen preserves "
        "the intact complex. Reported as a representative case of the "
        "broader detergent-affinity prediction problem for membrane proteins."
    ),
    "source_pmcid": None,
    "source_doi": None,
    "source_type": "personal_communication",
}

with open(SRC) as f:
    records = json.load(f)
print(f"Loaded {len(records)} records from {SRC.name}")

if any(r["uniprot_id"] == "P34972" for r in records):
    raise SystemExit("P34972 already present; aborting")

records.append(cb2_record)
print(f"Appended CB2 P34972 -> {len(records)} records")

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
