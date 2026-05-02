"""Harvest Europe PMC papers that are likely to contain extractable native MS protein data
in their supplementary tables, ranked for human triage.

Output: a queue (CSV + JSON) of papers with:
- title, authors, year, journal, DOI, PMID, PMCID
- open access flag
- abstract (so the user can quickly skim relevance)
- direct PMC link to the paper (where supplements are accessible)
- relevance score (recency + open access + has PMCID + ADC/antibody mention)

The user (or an undergrad) opens each PMC link, downloads the supplements,
and extracts protein-level data into the schema.

Scope: prioritizes the antibody/ADC/biopharma slice because that's where supplementary
tables most often contain explicit protein lists with native MS outcomes.
"""

import json
import time
import csv
from urllib.parse import quote
import requests

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
OUT_JSON = "/Users/bfentaw2/startup/nativeready/data/europepmc_queue_2026-05-01.json"
OUT_CSV = "/Users/bfentaw2/startup/nativeready/data/europepmc_queue_2026-05-01.csv"

# Targeted queries. Order matters: highest-relevance first.
QUERIES = [
    ('"antibody-drug conjugate" "native mass spectrometry"', "ADC native MS", 5),
    ('ADC DAR "native mass spectrometry"', "ADC DAR analysis", 5),
    ('bispecific "native mass spectrometry"', "Bispecific native MS", 4),
    ('antibody "native mass spectrometry" supplementary', "Antibody native MS w/ supplements", 4),
    ('"intact mass" antibody supplementary', "Intact mass antibody w/ supplements", 3),
    ('"native top-down" proteoform', "Native top-down", 3),
    ('"native mass spectrometry" "supplementary table" protein', "General native MS w/ supp tables", 2),
    ('"membrane protein" "native mass spectrometry"', "Membrane protein native MS", 2),
    ('"viral capsid" "native mass spectrometry"', "Viral capsid native MS", 2),
    ('"protein complex" "native mass spectrometry" stoichiometry', "Protein complex native MS", 2),
]

YEAR_MIN = 2018  # more recent = more likely to have machine-readable supplements

def search(query, page_size=100, max_pages=3):
    """Yield paper dicts from Europe PMC."""
    cursor_mark = "*"
    seen = set()
    for page in range(max_pages):
        params = {
            "query": f"{query} AND PUB_YEAR:[{YEAR_MIN} TO 2026]",
            "resultType": "core",
            "format": "json",
            "pageSize": page_size,
            "cursorMark": cursor_mark,
        }
        try:
            r = requests.get(EUROPEPMC, params=params, timeout=60)
            if r.status_code != 200:
                print(f"  Europe PMC returned {r.status_code} on page {page}")
                return
            data = r.json()
        except Exception as e:
            print(f"  Europe PMC failed: {e}")
            return

        results = data.get("resultList", {}).get("result", []) or []
        next_cursor = data.get("nextCursorMark") or "*"
        if not results:
            return
        for paper in results:
            pid = paper.get("id") or paper.get("pmid") or paper.get("doi") or ""
            if pid in seen:
                continue
            seen.add(pid)
            yield paper
        if next_cursor == cursor_mark:
            return
        cursor_mark = next_cursor


def normalize_paper(paper, query_label, query_priority):
    """Extract the fields we care about into a flat dict."""
    pmid = paper.get("pmid") or ""
    pmcid = paper.get("pmcid") or ""
    doi = paper.get("doi") or ""
    title = paper.get("title", "").rstrip(".")
    journal = paper.get("journalTitle", "")
    year = paper.get("pubYear", "")
    abstract = paper.get("abstractText", "") or ""
    open_access = paper.get("isOpenAccess") in ("Y", "y", True)
    has_supp = paper.get("hasSuppl") in ("Y", "y", True)
    has_pdf = paper.get("hasPDF") in ("Y", "y", True)
    has_text_mined = paper.get("hasTextMinedTerms") in ("Y", "y", True)
    author_str = paper.get("authorString", "")

    # Build PMC link if available (supplements live there)
    if pmcid:
        pmc_url = f"https://europepmc.org/article/PMC/{pmcid.replace('PMC', '')}"
    elif pmid:
        pmc_url = f"https://europepmc.org/article/MED/{pmid}"
    elif doi:
        pmc_url = f"https://doi.org/{doi}"
    else:
        pmc_url = ""

    # Relevance score
    score = query_priority
    if open_access:
        score += 2
    if has_supp:
        score += 3  # explicit supplementary indicator
    if pmcid:
        score += 1
    if has_pdf:
        score += 1
    if has_text_mined:
        score += 1
    try:
        yr = int(year)
        score += min(3, max(0, yr - 2018))  # recency bonus
    except Exception:
        pass

    return {
        "score": score,
        "query_label": query_label,
        "year": year,
        "title": title,
        "journal": journal,
        "authors": author_str,
        "doi": doi,
        "pmid": pmid,
        "pmcid": pmcid,
        "pmc_url": pmc_url,
        "open_access": open_access,
        "has_supplementary": has_supp,
        "has_pdf": has_pdf,
        "has_text_mined": has_text_mined,
        "abstract": abstract[:600],
    }


print("=" * 70)
print("EUROPE PMC NATIVE MS PAPER HARVEST")
print("=" * 70)

all_papers = {}  # key by (pmid or doi) to dedupe
per_query_counts = {}

for query, label, priority in QUERIES:
    print(f"\nQuery: {label}")
    print(f"  EuropePMC q: {query}")
    n_added = 0
    for paper in search(query):
        norm = normalize_paper(paper, label, priority)
        key = norm["pmid"] or norm["doi"] or norm["title"]
        if key in all_papers:
            # If duplicate, keep the higher-scoring entry
            if norm["score"] > all_papers[key]["score"]:
                all_papers[key] = norm
            continue
        all_papers[key] = norm
        n_added += 1
    per_query_counts[label] = n_added
    print(f"  Returned {n_added} unique papers added to queue")
    time.sleep(0.5)

# Sort by score desc
queue = sorted(all_papers.values(), key=lambda x: -x["score"])

# Save JSON (full metadata)
with open(OUT_JSON, "w") as f:
    json.dump(queue, f, indent=2)

# Save CSV (compact triage view)
csv_fields = [
    "score", "year", "query_label", "title", "journal", "authors",
    "open_access", "has_supplementary", "has_pdf", "pmcid", "pmid", "doi", "pmc_url"
]
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
    w.writeheader()
    for p in queue:
        w.writerow(p)

print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print(f"Total unique papers in queue: {len(queue)}")
print(f"  Open access: {sum(1 for p in queue if p['open_access'])}")
print(f"  Has supplementary flag: {sum(1 for p in queue if p['has_supplementary'])}")
print(f"  Has PMCID (supplements likely accessible): {sum(1 for p in queue if p['pmcid'])}")
print(f"\nPer-query breakdown:")
for label, n in per_query_counts.items():
    print(f"  {label:42} {n:5}")

print(f"\nFiles written:")
print(f"  {OUT_JSON}")
print(f"  {OUT_CSV}")

print(f"\nTop 10 papers by relevance score:")
for p in queue[:10]:
    print(f"  [{p['score']}] {p['year']} {p['journal'][:30]:30} {p['title'][:90]}")
