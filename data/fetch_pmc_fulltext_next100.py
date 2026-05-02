#!/usr/bin/env python3
"""Fetch fulltext XML for papers ranked 31-130 in the EuropePMC queue.

Same logic as fetch_pmc_fulltext.py but starts at offset 30 to skip already-processed.
Saves into the same europepmc_fulltext/ directory.
"""

import json
import os
import csv
import re
import time
import warnings
warnings.filterwarnings("ignore")

import requests

DATA = "/Users/bfentaw2/startup/nativeready/data"
QUEUE = f"{DATA}/europepmc_queue_2026-05-01.json"
OUT_DIR = f"{DATA}/europepmc_fulltext"
INDEX_CSV = f"{DATA}/europepmc_extraction_index_next100.csv"

START_OFFSET = 30
TOP_N = 100

EPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC{pmcid_num}/fullTextXML"
EPMC_ANNOT = ("https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds"
              "?articleIds=PMC%3A{pmcid}&type=Gene_Proteins&format=JSON")

UNIPROT_RE = re.compile(r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b|\b[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}\b")

NATIVE_MS_KEYWORDS = [
    "native mass spectrometry", "native MS", "native ESI", "intact mass",
    "non-covalent mass spectrometry", "ion mobility mass spectrometry",
    "native top-down", "DAR analysis", "Q-Exactive UHMR", "Orbitrap UHMR",
    "ammonium acetate", "supercharging", "charge state distribution",
    "native condition", "non-denaturing"
]
NATIVE_MS_RE = re.compile("|".join(re.escape(k) for k in NATIVE_MS_KEYWORDS), re.IGNORECASE)

FAILURE_KEYWORDS = [
    "did not", "could not", "failed", "no signal", "no peak", "did not yield",
    "did not ionize", "denatured", "aggregat", "insoluble",
    "could not be detected", "did not detect", "unsuccessful",
    "discarded", "rejected", "not observed",
]
FAILURE_RE = re.compile("|".join(re.escape(k) for k in FAILURE_KEYWORDS), re.IGNORECASE)

os.makedirs(OUT_DIR, exist_ok=True)

with open(QUEUE) as f:
    all_papers = json.load(f)

candidates = [p for p in all_papers if p.get("pmcid") and p.get("open_access")]
candidates.sort(key=lambda x: -x["score"])
batch = candidates[START_OFFSET:START_OFFSET + TOP_N]
print(f"Loaded {len(all_papers)} papers; {len(candidates)} OA+PMCID")
print(f"Processing offset {START_OFFSET} to {START_OFFSET + TOP_N - 1} ({len(batch)} papers)")


def extract_supp_file_refs(xml_text):
    refs = []
    for m in re.finditer(r'xlink:href="([^"]+\.(?:xlsx|xls|csv|docx|doc|pdf|txt|zip))"', xml_text, re.IGNORECASE):
        refs.append(m.group(1))
    for m in re.finditer(r'href="([^"]+\.(?:xlsx|xls|csv|docx|doc|pdf|txt|zip))"', xml_text, re.IGNORECASE):
        if m.group(1) not in refs:
            refs.append(m.group(1))
    return refs


def extract_native_ms_context_paragraphs(xml_text, max_paras=8):
    paragraphs = re.findall(r"<p[^>]*>(.*?)</p>", xml_text, re.DOTALL)
    flagged = []
    for p in paragraphs:
        text = re.sub(r"<[^>]+>", " ", p)
        text = re.sub(r"\s+", " ", text).strip()
        if NATIVE_MS_RE.search(text) and len(text) > 50:
            flagged.append(text[:500])
        if len(flagged) >= max_paras:
            break
    return flagged


def extract_failure_paragraphs(xml_text, max_paras=5):
    paragraphs = re.findall(r"<p[^>]*>(.*?)</p>", xml_text, re.DOTALL)
    flagged = []
    for p in paragraphs:
        text = re.sub(r"<[^>]+>", " ", p)
        text = re.sub(r"\s+", " ", text).strip()
        if FAILURE_RE.search(text) and NATIVE_MS_RE.search(text) and len(text) > 50:
            flagged.append(text[:500])
        if len(flagged) >= max_paras:
            break
    return flagged


def find_uniprot_in_text(xml_text):
    plain = re.sub(r"<[^>]+>", " ", xml_text)
    return sorted(set(UNIPROT_RE.findall(plain)))


def get_annotated_proteins(pmcid):
    try:
        r = requests.get(EPMC_ANNOT.format(pmcid=pmcid), timeout=30)
        if r.status_code != 200:
            return []
        data = r.json()
        proteins = set()
        for art in data:
            for ann in art.get("annotations", []):
                if ann.get("type") in ("Gene_Proteins", "gene_protein"):
                    name = ann.get("exact") or ann.get("name") or ""
                    if name:
                        proteins.add(name)
        return sorted(proteins)[:50]
    except Exception:
        return []


index_rows = []
n_ok = 0
n_skip_existing = 0
n_fail = 0

for i, paper in enumerate(batch):
    pmcid_raw = paper.get("pmcid", "")
    pmcid_num = pmcid_raw.replace("PMC", "")
    paper_dir = os.path.join(OUT_DIR, f"PMC{pmcid_num}")

    # Skip if already downloaded
    if os.path.exists(os.path.join(paper_dir, "fulltext.xml")):
        n_skip_existing += 1
        if (i + 1) % 20 == 0:
            print(f"  [{i+1}/{TOP_N}] checked; ok={n_ok} skip_existing={n_skip_existing} fail={n_fail}")
        continue

    os.makedirs(paper_dir, exist_ok=True)
    with open(os.path.join(paper_dir, "metadata.json"), "w") as f:
        json.dump(paper, f, indent=2)

    url = EPMC_FULLTEXT.format(pmcid_num=pmcid_num)
    try:
        r = requests.get(url, timeout=60)
        if r.status_code != 200:
            n_fail += 1
            index_rows.append({
                "pmcid": pmcid_raw, "title": paper.get("title", "")[:80],
                "score": paper.get("score"),
                "fulltext_status": f"http_{r.status_code}",
                "fulltext_size_bytes": 0,
                "supp_file_refs": "", "native_ms_paragraphs": 0,
                "failure_paragraphs": 0, "uniprot_ids_in_text": "",
                "annotated_proteins": "",
                "manual_download_url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid_raw}/",
            })
            continue
        xml_text = r.text
        with open(os.path.join(paper_dir, "fulltext.xml"), "w") as f:
            f.write(xml_text)
    except Exception as e:
        n_fail += 1
        continue

    supp_refs = extract_supp_file_refs(xml_text)
    nms_paras = extract_native_ms_context_paragraphs(xml_text)
    fail_paras = extract_failure_paragraphs(xml_text)
    uids = find_uniprot_in_text(xml_text)
    proteins = get_annotated_proteins(pmcid_raw)

    extract = {
        "pmcid": pmcid_raw,
        "title": paper.get("title"),
        "doi": paper.get("doi"),
        "supp_file_refs": supp_refs,
        "native_ms_paragraphs": nms_paras,
        "failure_paragraphs": fail_paras,
        "uniprot_ids_in_text": uids,
        "annotated_proteins": proteins,
    }
    with open(os.path.join(paper_dir, "extract.json"), "w") as f:
        json.dump(extract, f, indent=2)

    index_rows.append({
        "pmcid": pmcid_raw,
        "title": paper.get("title", "")[:80],
        "score": paper.get("score"),
        "fulltext_status": "ok",
        "fulltext_size_bytes": len(xml_text),
        "supp_file_refs": "; ".join(supp_refs),
        "native_ms_paragraphs": len(nms_paras),
        "failure_paragraphs": len(fail_paras),
        "uniprot_ids_in_text": "; ".join(uids[:10]),
        "annotated_proteins": "; ".join(proteins[:10]),
        "manual_download_url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid_raw}/",
    })
    n_ok += 1

    if (i + 1) % 20 == 0:
        print(f"  [{i+1}/{TOP_N}] ok={n_ok} skip_existing={n_skip_existing} fail={n_fail}")
    time.sleep(0.3)

# Save index
fields = ["pmcid", "score", "fulltext_status", "fulltext_size_bytes",
          "native_ms_paragraphs", "failure_paragraphs", "uniprot_ids_in_text",
          "annotated_proteins", "supp_file_refs", "manual_download_url", "title"]
with open(INDEX_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    for row in index_rows:
        w.writerow(row)

print("\n" + "=" * 70)
print(f"FETCH COMPLETE")
print(f"  Downloaded fulltext: {n_ok}")
print(f"  Skipped (already present): {n_skip_existing}")
print(f"  Failed: {n_fail}")
print(f"  Index: {INDEX_CSV}")

# Summary stats from new ones
ok_rows = [r for r in index_rows if r["fulltext_status"] == "ok"]
if ok_rows:
    print(f"\n  Of {len(ok_rows)} new papers:")
    print(f"    With native MS paragraphs: {sum(1 for r in ok_rows if r['native_ms_paragraphs'] > 0)}")
    print(f"    With failure-flagged paragraphs: {sum(1 for r in ok_rows if r['failure_paragraphs'] > 0)}")
    print(f"    With UniProt IDs in text: {sum(1 for r in ok_rows if r['uniprot_ids_in_text'])}")
    print(f"    With supp file references: {sum(1 for r in ok_rows if r['supp_file_refs'])}")
