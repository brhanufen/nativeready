#!/usr/bin/env python3
"""Extract native MS proteins from Zenodo by mining UniProt IDs from record descriptions."""
import json
import re
import time
import warnings
warnings.filterwarnings("ignore")

import requests

DATA_DIR = "/Users/bfentaw2/startup/nativeready/data"
EXISTING_IDS = set(json.load(open("/tmp/existing_uniprot_ids.json")))

try:
    pdb_records = json.load(open(f"{DATA_DIR}/pdb_extracted.json"))
    PDB_IDS = {r["uniprot_id"] for r in pdb_records}
except Exception:
    PDB_IDS = set()

try:
    pride_records = json.load(open(f"{DATA_DIR}/pride_extracted.json"))
    PRIDE_IDS = {r["uniprot_id"] for r in pride_records}
except Exception:
    PRIDE_IDS = set()

ALL_EXISTING = EXISTING_IDS | PDB_IDS | PRIDE_IDS

print(f"Dedup set: existing={len(EXISTING_IDS)} pdb={len(PDB_IDS)} pride={len(PRIDE_IDS)} total={len(ALL_EXISTING)}")

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
ZENODO_URL = "https://zenodo.org/api/records"
UNIPROT_RE = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\b")

QUERIES = [
    '"native mass spectrometry"',
    '"native top-down"',
    '"intact mass spectrometry"'
]


def classify_protein(udata, name_text):
    name_low = name_text.lower()
    keywords = []
    for kw in udata.get("keywords", []):
        if isinstance(kw, dict):
            keywords.append(kw.get("name", "").lower())
        elif isinstance(kw, str):
            keywords.append(kw.lower())
    kw_text = " ".join(keywords)
    full = name_low + " " + kw_text
    if "antibody" in full or "immunoglobulin" in full:
        return "antibody"
    if "kinase" in full:
        return "kinase"
    if "g protein-coupled" in full or "gpcr" in full:
        return "GPCR"
    if "ion channel" in full or "channel" in kw_text:
        return "ion channel"
    if "membrane" in kw_text or "transmembrane" in kw_text:
        return "membrane protein"
    if "capsid" in full or "virus" in full or "viral" in full:
        return "viral capsid"
    if "intrinsically disordered" in full or "disordered" in full:
        return "IDP"
    if "ribosom" in full:
        return "complex"
    if "chaperon" in full:
        return "chaperone"
    if "transferase" in full or "hydrolase" in full or "ligase" in full or "lyase" in full or "isomerase" in full or "oxidoreductase" in full or "enzyme" in kw_text:
        return "enzyme"
    return "protein"


def fetch_uniprot(uid):
    try:
        r = requests.get(UNIPROT_URL.format(uid), timeout=30)
        if r.status_code != 200:
            return None
        return r.json()
    except Exception:
        return None


def get_organism(udata):
    org = udata.get("organism", {})
    sci = org.get("scientificName", "")
    common = org.get("commonName", "")
    return f"{sci} ({common})" if common else sci


def get_name(udata):
    pd = udata.get("proteinDescription", {})
    rec = pd.get("recommendedName", {})
    full = rec.get("fullName", {})
    if isinstance(full, dict):
        return full.get("value", "")
    return ""


extracted = []
seen_uids = set()
records_examined = 0
records_with_uniprot_text = 0
duplicates_against_existing = 0

for q in QUERIES:
    print(f"\nZenodo query: {q}")
    try:
        # Zenodo caps size at 25; we'd paginate but for our query the total is well under 50
        q_plus = q.replace(" ", "+")
        url = f'{ZENODO_URL}?q=%22{q_plus}%22&size=25'
        r = requests.get(url, headers={"User-Agent": "NativeReadyBot/1.0"}, timeout=60)
        if r.status_code != 200:
            print(f"  Zenodo returned {r.status_code}")
            continue
        data = r.json()
    except Exception as e:
        print(f"  Zenodo failed: {e}")
        continue

    hits = data.get("hits", {}).get("hits", [])
    print(f"  Zenodo returned {len(hits)} records")

    for rec in hits:
        records_examined += 1
        meta = rec.get("metadata", {})
        title = meta.get("title", "")
        desc = meta.get("description", "") or ""
        # strip simple HTML
        desc_clean = re.sub(r"<[^>]+>", " ", desc)
        keywords = " ".join(meta.get("keywords", []) or [])
        doi = rec.get("doi") or meta.get("doi") or ""
        text_blob = " ".join([title, desc_clean, keywords])
        matches = UNIPROT_RE.findall(text_blob)
        flat_ids = set()
        for m in matches:
            if isinstance(m, tuple):
                flat_ids.add(m[0])
            else:
                flat_ids.add(m)
        if flat_ids:
            records_with_uniprot_text += 1

        for uid in flat_ids:
            if uid in ALL_EXISTING:
                duplicates_against_existing += 1
                continue
            if uid in seen_uids:
                continue
            udata = fetch_uniprot(uid)
            if not udata:
                continue
            seq = (udata.get("sequence") or {}).get("value", "")
            mw = (udata.get("sequence") or {}).get("molWeight", 0)
            if not seq:
                continue
            name = get_name(udata)
            if not name:
                continue
            organism = get_organism(udata)
            full_name = f"{name} ({organism})" if organism else name
            pclass = classify_protein(udata, name)
            notes = f"Source: Zenodo DOI {doi}; {title}".strip()
            record = {
                "uniprot_id": uid,
                "name": full_name,
                "protein_class": pclass,
                "mw_kda": round(mw / 1000.0, 2) if mw else 0.0,
                "notes": notes,
                "sequence": seq,
                "sequence_length": len(seq),
                "label": 1
            }
            extracted.append(record)
            seen_uids.add(uid)
            time.sleep(0.05)

with open(f"{DATA_DIR}/zenodo_extracted.json", "w") as f:
    json.dump(extracted, f, indent=2)

print("=" * 60)
print(f"ZENODO FINAL: extracted={len(extracted)}")
print(f"  records examined: {records_examined}")
print(f"  records with UniProt IDs in text: {records_with_uniprot_text}")
print(f"  duplicates rejected: {duplicates_against_existing}")
