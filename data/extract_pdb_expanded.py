#!/usr/bin/env python3
"""Expanded PDB extraction with broader native-MS-adjacent queries.

Runs 7 expanded full-text queries against RCSB PDB and dedupes against
the 492 unique UniProts already in the dataset (232 original + 260 from
the May-1 PDB extraction).
"""
import json
import os
import time
import warnings
warnings.filterwarnings("ignore")

import requests

DATA_DIR = "/Users/bfentaw2/startup/nativeready/data"

# Load all existing UniProt IDs from the 4 dataset files
EXISTING_FILES = [
    "positives_with_sequences.json",
    "negatives_with_sequences.json",
    "expansion_with_sequences.json",
    "extracted_new_2026-05-01.json",
]
EXISTING_IDS = set()
for fn in EXISTING_FILES:
    with open(os.path.join(DATA_DIR, fn)) as f:
        for r in json.load(f):
            if r.get("uniprot_id"):
                EXISTING_IDS.add(r["uniprot_id"])
print(f"Loaded {len(EXISTING_IDS)} existing UniProt IDs to dedup against")

# Expanded query list
QUERIES = [
    "native ESI",
    "intact mass spectrometry",
    "ion mobility mass spectrometry",
    "native top-down",
    "non-covalent mass spectrometry",
    "charge detection mass spectrometry",
    "individual ion mass spectrometry",
    "variable temperature electrospray",
]

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{}"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"


def classify_protein(uniprot_data, name_text):
    name_low = name_text.lower()
    keywords = []
    for kw in uniprot_data.get("keywords", []):
        if isinstance(kw, dict):
            keywords.append(kw.get("name", "").lower())
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
    except Exception:
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


def search_pdb(query_text, max_rows=400):
    body = {
        "query": {"type": "terminal", "service": "full_text",
                  "parameters": {"value": query_text}},
        "request_options": {
            "results_content_type": ["experimental"],
            "paginate": {"start": 0, "rows": max_rows},
        },
        "return_type": "polymer_entity",
    }
    r = requests.post(SEARCH_URL, json=body, timeout=60)
    if r.status_code != 200:
        return 0, []
    data = r.json()
    return data.get("total_count", 0), data.get("result_set", [])


def fetch_entry_title(entry_id):
    try:
        rr = requests.get(ENTRY_URL.format(entry_id), timeout=30)
        if rr.status_code != 200:
            return ""
        ej = rr.json()
        cits = ej.get("citation", [])
        if cits:
            return cits[0].get("title", "")
    except Exception:
        pass
    return ""


extracted = []
seen_uids = set()
per_query_yield = {}
duplicates_against_existing = 0
duplicates_within_run = 0

for query in QUERIES:
    print(f"\nQuery: '{query}'")
    total, hits = search_pdb(query)
    print(f"  PDB total_count={total}, processing {len(hits)} polymer entities")
    n_added_this_query = 0
    for i, hit in enumerate(hits):
        pe_id = hit.get("identifier")
        if not pe_id or "_" not in pe_id:
            continue
        entry_id, entity_id = pe_id.split("_", 1)
        try:
            er = requests.get(ENTITY_URL.format(entry_id, entity_id), timeout=30)
            if er.status_code != 200:
                continue
            ent = er.json()
        except Exception:
            continue
        container = ent.get("rcsb_polymer_entity_container_identifiers", {}) or {}
        uniprot_ids = container.get("uniprot_ids") or []
        if not uniprot_ids:
            continue
        uid = uniprot_ids[0]
        if uid in EXISTING_IDS:
            duplicates_against_existing += 1
            continue
        if uid in seen_uids:
            duplicates_within_run += 1
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
            name = (ent.get("rcsb_polymer_entity") or {}).get("pdbx_description", "") or uid
        organism = get_organism(udata)
        full_name = f"{name} ({organism})" if organism else name
        pclass = classify_protein(udata, name)
        cit_title = fetch_entry_title(entry_id)
        notes = f"Source: PDB {entry_id} (query: '{query}'); native MS-adjacent; {cit_title or 'cited in PDB entry'}"

        record = {
            "uniprot_id": uid,
            "name": full_name,
            "protein_class": pclass,
            "mw_kda": round(mw / 1000.0, 2) if mw else 0.0,
            "notes": notes,
            "sequence": seq,
            "sequence_length": len(seq),
            "label": 1,
            "source_query": query,
        }
        extracted.append(record)
        seen_uids.add(uid)
        n_added_this_query += 1

        if (i + 1) % 50 == 0:
            print(f"  ... processed {i+1}/{len(hits)}; query_added={n_added_this_query}; total_run={len(extracted)}")
            # partial save
            with open(f"{DATA_DIR}/extracted_pdb_expanded_2026-05-02.json", "w") as f:
                json.dump(extracted, f, indent=2)
        time.sleep(0.05)

    per_query_yield[query] = n_added_this_query
    print(f"  Done with '{query}': +{n_added_this_query} new records")

# Final save
with open(f"{DATA_DIR}/extracted_pdb_expanded_2026-05-02.json", "w") as f:
    json.dump(extracted, f, indent=2)

print("\n" + "=" * 60)
print("EXPANDED PDB EXTRACTION FINAL SUMMARY")
print("=" * 60)
print(f"Total NEW records added (after dedup): {len(extracted)}")
print(f"Duplicates against existing 492: {duplicates_against_existing}")
print(f"Duplicates within this run: {duplicates_within_run}")
print(f"\nPer-query yield:")
for q, n in per_query_yield.items():
    print(f"  {q:40} +{n}")
print(f"\nProjected total dataset size: {len(EXISTING_IDS)} + {len(extracted)} = {len(EXISTING_IDS) + len(extracted)}")
print(f"\nFile: {DATA_DIR}/extracted_pdb_expanded_2026-05-02.json")
