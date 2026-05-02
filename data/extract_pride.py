#!/usr/bin/env python3
"""Extract native MS proteins from PRIDE Archive by mining UniProt IDs from project descriptions."""
import json
import re
import time
import warnings
warnings.filterwarnings("ignore")

import requests

DATA_DIR = "/Users/bfentaw2/startup/nativeready/data"
EXISTING_IDS = set(json.load(open("/tmp/existing_uniprot_ids.json")))

# Load PDB extracted to also dedup against (avoid double-add)
try:
    pdb_records = json.load(open(f"{DATA_DIR}/pdb_extracted.json"))
    PDB_IDS = {r["uniprot_id"] for r in pdb_records}
    print(f"Loaded {len(PDB_IDS)} PDB-extracted IDs to also dedup against")
except Exception:
    PDB_IDS = set()

ALL_EXISTING = EXISTING_IDS | PDB_IDS

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
PRIDE_SEARCH = "https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects"
# UniProt accession regex
UNIPROT_RE = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\b")

KEYWORDS = [
    "native mass spectrometry",
    "native top-down",
    "native ESI",
    "intact mass spectrometry",
    "HDX-MS"
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
projects_examined = 0
projects_with_uniprot_text = 0
duplicates_against_existing = 0

for keyword in KEYWORDS:
    print(f"\nKeyword: '{keyword}'")
    try:
        r = requests.get(PRIDE_SEARCH, params={"keyword": keyword, "pageSize": 100}, timeout=60)
        if r.status_code != 200:
            print(f"  PRIDE returned {r.status_code}")
            continue
        data = r.json()
    except Exception as e:
        print(f"  PRIDE search failed: {e}")
        continue

    if isinstance(data, list):
        projects = data
    else:
        projects = data.get("_embedded", {}).get("compactprojects", []) or data.get("data", []) or []
    print(f"  PRIDE returned {len(projects)} projects")

    for proj in projects:
        projects_examined += 1
        accession = proj.get("accession") or proj.get("projectAccession") or ""
        title = proj.get("title", "")
        desc = proj.get("projectDescription", "") or proj.get("description", "")
        sample = proj.get("sampleProcessingProtocol", "") or ""
        keywords_field = " ".join(proj.get("keywords", []) or []) if isinstance(proj.get("keywords"), list) else (proj.get("keywords") or "")

        text_blob = " ".join([title, desc, sample, keywords_field])
        matches = UNIPROT_RE.findall(text_blob)
        # findall returns tuples for grouped patterns; flatten
        flat_ids = set()
        for m in matches:
            if isinstance(m, tuple):
                flat_ids.add(m[0])
            else:
                flat_ids.add(m)

        if flat_ids:
            projects_with_uniprot_text += 1

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
            notes = f"Source: PRIDE {accession}; native MS submission; {title}".strip()
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

with open(f"{DATA_DIR}/pride_extracted.json", "w") as f:
    json.dump(extracted, f, indent=2)

print("=" * 60)
print(f"PRIDE FINAL: extracted={len(extracted)}")
print(f"  projects examined: {projects_examined}")
print(f"  projects with extractable UniProt IDs in text: {projects_with_uniprot_text}")
print(f"  duplicates rejected (against existing+PDB): {duplicates_against_existing}")
