#!/usr/bin/env python3
"""Fetch fulltext XML from EuropePMC for top-N native MS papers.

Why this approach:
- The PMC binary supplementary files are protected by a JavaScript proof-of-work
  challenge that breaks scripted downloads.
- EuropePMC's fullTextXML endpoint serves the JATS XML for OA papers without PoW.
- JATS XML often contains:
    - Inline supplementary tables and figure captions
    - References to supplementary files (sec id="supplementary-material" / supplementary-material/media)
    - Protein names mentioned in methods/results text
- We also save a "manual download" CSV listing supp files the user can click
  through in a browser (browser sessions solve the PoW automatically).
"""

import json
import os
import csv
import re
import time
import warnings
warnings.filterwarnings("ignore")

import requests
import xml.etree.ElementTree as ET

DATA = "/Users/bfentaw2/startup/nativeready/data"
QUEUE = f"{DATA}/europepmc_queue_2026-05-01.json"
OUT_DIR = f"{DATA}/europepmc_fulltext"
INDEX_CSV = f"{DATA}/europepmc_extraction_index.csv"

TOP_N = 30  # bumped from 20

EPMC_FULLTEXT = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC{pmcid_num}/fullTextXML"
EPMC_ANNOT = ("https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds"
              "?articleIds=PMC%3A{pmcid}&type=Gene_Proteins&format=JSON")

# UniProt accession regex (canonical form)
UNIPROT_RE = re.compile(r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b|\b[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}\b")

# Native MS context keywords - flag papers/sections that mention these near a protein
NATIVE_MS_KEYWORDS = [
    "native mass spectrometry", "native MS", "native ESI", "intact mass",
    "non-covalent mass spectrometry", "ion mobility mass spectrometry",
    "native top-down", "DAR analysis", "Q-Exactive UHMR", "Orbitrap UHMR",
    "ammonium acetate", "supercharging", "charge state distribution",
    "native condition", "non-denaturing"
]
NATIVE_MS_RE = re.compile("|".join(re.escape(k) for k in NATIVE_MS_KEYWORDS), re.IGNORECASE)

# Failure-mode keywords to flag potential negatives
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
top_papers = candidates[:TOP_N]
print(f"Loaded {len(all_papers)} papers; {len(candidates)} OA+PMCID; processing top {TOP_N}")


def strip_xml_tags(text):
    """Quick text extraction from XML element."""
    if text is None:
        return ""
    return re.sub(r"<[^>]+>", " ", text)


def extract_supp_file_refs(xml_text):
    """Find supplementary file references in the JATS XML."""
    refs = []
    # JATS pattern: <supplementary-material> with <media xlink:href="..." or <self-uri xlink:href="...
    for m in re.finditer(r'xlink:href="([^"]+\.(?:xlsx|xls|csv|docx|doc|pdf|txt|zip))"', xml_text, re.IGNORECASE):
        refs.append(m.group(1))
    # Sometimes hrefs without xlink prefix
    for m in re.finditer(r'href="([^"]+\.(?:xlsx|xls|csv|docx|doc|pdf|txt|zip))"', xml_text, re.IGNORECASE):
        if m.group(1) not in refs:
            refs.append(m.group(1))
    return refs


def extract_native_ms_context_paragraphs(xml_text, max_paras=8):
    """Pull paragraphs that mention native MS keywords - these likely contain
    the actual experimental outcomes (clean spectrum, failure, etc.)."""
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
    """Paragraphs mentioning failure keywords - high-priority for negatives."""
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
    """Find UniProt accessions mentioned in the paper text."""
    plain = re.sub(r"<[^>]+>", " ", xml_text)
    return sorted(set(UNIPROT_RE.findall(plain)))


def get_annotated_proteins(pmcid):
    """Use EuropePMC text-mining annotations API for gene/protein mentions."""
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
        return sorted(proteins)[:50]  # cap to avoid noise
    except Exception:
        return []


index_rows = []

for i, paper in enumerate(top_papers):
    pmcid_raw = paper.get("pmcid", "")
    pmcid_num = pmcid_raw.replace("PMC", "")
    title = paper.get("title", "")[:80]
    print(f"\n[{i+1}/{TOP_N}] PMC{pmcid_num} (score={paper.get('score')})")
    print(f"  {title}")

    paper_dir = os.path.join(OUT_DIR, f"PMC{pmcid_num}")
    os.makedirs(paper_dir, exist_ok=True)

    # Save metadata
    with open(os.path.join(paper_dir, "metadata.json"), "w") as f:
        json.dump(paper, f, indent=2)

    # 1. Fetch fulltext XML
    url = EPMC_FULLTEXT.format(pmcid_num=pmcid_num)
    try:
        r = requests.get(url, timeout=120)
        if r.status_code != 200:
            print(f"  fullTextXML: {r.status_code}; skipping")
            index_rows.append({
                "pmcid": pmcid_raw, "title": title, "score": paper.get("score"),
                "fulltext_status": f"http_{r.status_code}",
                "fulltext_size_bytes": 0, "supp_file_refs": "",
                "native_ms_paragraphs": 0, "failure_paragraphs": 0,
                "uniprot_ids_in_text": "", "annotated_proteins": "",
                "manual_download_url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid_raw}/",
            })
            continue
        xml_text = r.text
        xml_size = len(xml_text)
        print(f"  fullTextXML: {xml_size} bytes")

        with open(os.path.join(paper_dir, "fulltext.xml"), "w") as f:
            f.write(xml_text)
    except Exception as e:
        print(f"  fullTextXML failed: {e}")
        continue

    # 2. Extract supplementary file references
    supp_refs = extract_supp_file_refs(xml_text)
    print(f"  Supp file refs in XML: {len(supp_refs)}")

    # 3. Extract native MS context paragraphs
    nms_paras = extract_native_ms_context_paragraphs(xml_text)
    print(f"  Native MS context paragraphs: {len(nms_paras)}")

    # 4. Extract failure-flagged paragraphs (potential negatives)
    fail_paras = extract_failure_paragraphs(xml_text)
    print(f"  Failure-flagged paragraphs: {len(fail_paras)}")

    # 5. UniProt IDs mentioned
    uids = find_uniprot_in_text(xml_text)
    print(f"  UniProt IDs in text: {len(uids)}")

    # 6. Text-mined protein/gene mentions
    proteins = get_annotated_proteins(pmcid_raw)
    print(f"  Annotated proteins (EuropePMC text mining): {len(proteins)}")

    # Save extracted highlights
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
        "title": title,
        "score": paper.get("score"),
        "fulltext_status": "ok",
        "fulltext_size_bytes": xml_size,
        "supp_file_refs": "; ".join(supp_refs),
        "native_ms_paragraphs": len(nms_paras),
        "failure_paragraphs": len(fail_paras),
        "uniprot_ids_in_text": "; ".join(uids[:10]),
        "annotated_proteins": "; ".join(proteins[:10]),
        "manual_download_url": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid_raw}/",
    })
    time.sleep(0.5)

# Save index
fields = [
    "pmcid", "score", "fulltext_status", "fulltext_size_bytes",
    "native_ms_paragraphs", "failure_paragraphs", "uniprot_ids_in_text",
    "annotated_proteins", "supp_file_refs", "manual_download_url", "title",
]
with open(INDEX_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    for row in index_rows:
        w.writerow(row)

# Summary
print("\n" + "=" * 70)
print("FULLTEXT EXTRACTION SUMMARY")
print("=" * 70)
ok = [r for r in index_rows if r["fulltext_status"] == "ok"]
print(f"Top {TOP_N} papers attempted")
print(f"Papers with fulltext XML retrieved: {len(ok)} / {TOP_N}")
print(f"Papers with native MS paragraphs found: {sum(1 for r in ok if r['native_ms_paragraphs'] > 0)}")
print(f"Papers with failure-flagged paragraphs: {sum(1 for r in ok if r['failure_paragraphs'] > 0)}")
print(f"Papers with UniProt IDs in text: {sum(1 for r in ok if r['uniprot_ids_in_text'])}")
print(f"Papers with EuropePMC-annotated protein mentions: {sum(1 for r in ok if r['annotated_proteins'])}")
print(f"Papers referencing supplementary files: {sum(1 for r in ok if r['supp_file_refs'])}")

print(f"\nIndex CSV: {INDEX_CSV}")
print(f"Per-paper extracts: {OUT_DIR}/PMC*/extract.json")
print(f"Per-paper fulltext XML: {OUT_DIR}/PMC*/fulltext.xml")
