#!/usr/bin/env python3
"""Extract native MS proteins from RCSB PDB."""
import json
import time
import sys
import warnings
warnings.filterwarnings("ignore")

import requests

DATA_DIR = "/Users/bfentaw2/startup/nativeready/data"
EXISTING_IDS = set(json.load(open("/tmp/existing_uniprot_ids.json")))
print(f"Loaded {len(EXISTING_IDS)} existing UniProt IDs to dedup against")

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"

# Step 2a: search PDB
search_body = {
    "query": {
        "type": "terminal",
        "service": "full_text",
        "parameters": {"value": "native mass spectrometry"}
    },
    "request_options": {
        "results_content_type": ["experimental"],
        "paginate": {"start": 0, "rows": 1198}
    },
    "return_type": "polymer_entity"
}

print("Searching PDB...")
r = requests.post(SEARCH_URL, json=search_body, timeout=60)
r.raise_for_status()
search_results = r.json()
total = search_results.get("total_count", 0)
hits = search_results.get("result_set", [])
print(f"PDB returned {total} total polymer_entity hits, processing {len(hits)}")

# Save raw search to a side file for traceability
with open(f"{DATA_DIR}/pdb_search_raw.json", "w") as f:
    json.dump(search_results, f, indent=2)


def classify_protein(uniprot_data, name_text):
    """Heuristic protein class from UniProt keywords + name."""
    name_low = name_text.lower()
    keywords = []
    for kw in uniprot_data.get("keywords", []):
        if isinstance(kw, dict):
            kw_name = kw.get("name", "").lower()
            keywords.append(kw_name)
        elif isinstance(kw, str):
            keywords.append(kw.lower())
    kw_text = " ".join(keywords)
    full = name_low + " " + kw_text

    if "antibody" in full or "immunoglobulin" in full or "fab " in full:
        return "antibody"
    if "kinase" in full:
        return "kinase"
    if "g-protein coupled" in full or "g protein-coupled" in full or "gpcr" in full:
        return "GPCR"
    if "ion channel" in full or "channel" in kw_text:
        return "ion channel"
    if "membrane" in kw_text or "transmembrane" in kw_text or "membrane protein" in full:
        return "membrane protein"
    if "capsid" in full or "viral" in full or "virus" in full or "nucleocapsid" in full:
        return "viral capsid"
    if "intrinsically disordered" in full or "disordered" in full:
        return "IDP"
    if "ribosom" in full:
        return "complex"
    if "chaperon" in full:
        return "chaperone"
    if "transferase" in full or "hydrolase" in full or "ligase" in full or "lyase" in full or "isomerase" in full or "oxidoreductase" in full or "enzyme" in kw_text:
        return "enzyme"
    if "transcription factor" in full:
        return "transcription factor"
    return "protein"


def fetch_uniprot(uid):
    try:
        r = requests.get(UNIPROT_URL.format(uid), timeout=30)
        if r.status_code != 200:
            return None
        return r.json()
    except Exception as e:
        print(f"  uniprot fetch fail {uid}: {e}")
        return None


def get_organism(udata):
    org = udata.get("organism", {})
    sci = org.get("scientificName", "")
    common = org.get("commonName", "")
    if common:
        return f"{sci} ({common})"
    return sci


def get_name(udata):
    pd = udata.get("proteinDescription", {})
    rec = pd.get("recommendedName", {})
    full = rec.get("fullName", {})
    if isinstance(full, dict):
        return full.get("value", "")
    return ""


extracted = []
seen_uids = set()
duplicates_against_existing = 0
duplicates_within_pdb = 0
no_uniprot_count = 0
errors = 0

for i, hit in enumerate(hits):
    pe_id = hit.get("identifier")  # e.g. "8UDN_1"
    if not pe_id or "_" not in pe_id:
        continue
    entry_id, entity_id = pe_id.split("_", 1)
    try:
        er = requests.get(ENTITY_URL.format(entry_id, entity_id), timeout=30)
        if er.status_code != 200:
            errors += 1
            continue
        ent = er.json()
    except Exception as e:
        errors += 1
        continue

    container = ent.get("rcsb_polymer_entity_container_identifiers", {}) or {}
    uniprot_ids = container.get("uniprot_ids") or []
    if not uniprot_ids:
        no_uniprot_count += 1
        continue
    uid = uniprot_ids[0]
    if uid in EXISTING_IDS:
        duplicates_against_existing += 1
        continue
    if uid in seen_uids:
        duplicates_within_pdb += 1
        continue

    udata = fetch_uniprot(uid)
    if not udata:
        errors += 1
        continue
    seq = (udata.get("sequence") or {}).get("value", "")
    mw = (udata.get("sequence") or {}).get("molWeight", 0)
    if not seq:
        no_uniprot_count += 1
        continue
    name = get_name(udata)
    if not name:
        # fallback to entity-level name
        name = (ent.get("rcsb_polymer_entity") or {}).get("pdbx_description", "") or uid
    organism = get_organism(udata)
    full_name = f"{name} ({organism})" if organism else name

    pclass = classify_protein(udata, name)

    # Try to get paper title from entity
    cit_title = ""
    refs = ent.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])
    # Get from entry-level citation
    try:
        rr = requests.get(ENTRY_URL.format(entry_id), timeout=30)
        if rr.status_code == 200:
            ej = rr.json()
            citations = ej.get("citation", [])
            if citations:
                cit_title = citations[0].get("title", "")
    except Exception:
        pass

    notes = f"Source: PDB {entry_id}; native MS-validated; {cit_title or 'cited in PDB entry'}"

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

    if (i + 1) % 25 == 0:
        print(f"  processed {i+1}/{len(hits)}; new={len(extracted)} dup_existing={duplicates_against_existing} dup_within={duplicates_within_pdb} no_up={no_uniprot_count} err={errors}")
        # save partial progress
        with open(f"{DATA_DIR}/pdb_extracted.json", "w") as f:
            json.dump(extracted, f, indent=2)
    time.sleep(0.05)

with open(f"{DATA_DIR}/pdb_extracted.json", "w") as f:
    json.dump(extracted, f, indent=2)

print("=" * 60)
print(f"PDB FINAL: extracted={len(extracted)}")
print(f"  duplicates vs existing 232: {duplicates_against_existing}")
print(f"  duplicates within PDB results: {duplicates_within_pdb}")
print(f"  no UniProt / no seq: {no_uniprot_count}")
print(f"  errors: {errors}")
