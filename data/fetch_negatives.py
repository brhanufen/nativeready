"""
Fetch proxy negative examples from UniProt — real proteins from the reviewed
Swiss-Prot subset that are NOT in our curated positive list.

Approach (honest about limitations):
- Sample real, reviewed proteins from public UniProt
- Bias the sample toward classes known to be harder for native MS:
  * Very large multi-pass membrane proteins
  * Highly glycosylated mucins / secreted proteins
  * Very long sequences (>1500 aa)
  * Highly disordered proteins
- Mix with randomly sampled smaller proteins for balance
- Exclude any protein in the positive list

This is a proxy: failures are rarely published, so we can't get true negatives
from the literature. By sampling diverse Swiss-Prot proteins not in the
positive set, we get a class of proteins that statistically includes many
hard-to-study examples. Acknowledge this limitation in the model card.
"""
import json
import time
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path
from typing import Optional, List, Dict

from positives_curated import POSITIVE_EXAMPLES

DATA_DIR = Path(__file__).parent
OUTPUT_NEGATIVES = DATA_DIR / "negatives_with_sequences.json"
POSITIVE_IDS = {p["uniprot_id"] for p in POSITIVE_EXAMPLES}

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"


def search_uniprot(query: str, limit: int = 25) -> List[Dict]:
    """Search Swiss-Prot reviewed entries with a query, return JSON."""
    params = {
        "query": query,
        "format": "json",
        "size": str(limit),
        "fields": "accession,id,protein_name,sequence,length,mass,cc_subcellular_location",
    }
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "NativeReady/0.1"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data.get("results", [])
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        print(f"  Search failed for '{query}': {e}")
        return []


def fetch_fasta(accession: str) -> Optional[str]:
    """Get the bare FASTA sequence."""
    url = UNIPROT_FASTA.format(accession)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "NativeReady/0.1"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            text = resp.read().decode("utf-8")
            lines = text.strip().split("\n")
            return "".join(lines[1:]).replace(" ", "").upper()
    except (urllib.error.URLError, urllib.error.HTTPError):
        return None


# Queries designed to enrich for proteins NOT typically targets of native MS
# (very large, very disordered, heavily glycosylated, multi-pass membrane).
# All queries restricted to Swiss-Prot reviewed entries (reviewed:true).
HARD_QUERIES = [
    # Very long proteins (>2000 aa)
    'reviewed:true AND length:[2000 TO 5000] AND organism_id:9606',
    # Multi-pass membrane proteins (likely transmembrane >7)
    'reviewed:true AND ft_transmem:* AND length:[400 TO 1500] AND organism_id:9606',
    # Highly glycosylated mucins
    'reviewed:true AND name:mucin AND organism_id:9606',
    # Intrinsically disordered
    'reviewed:true AND keyword:KW-1185 AND organism_id:9606',
    # Random small proteins for balance
    'reviewed:true AND length:[100 TO 300] AND organism_id:9606',
    # Random medium proteins for balance
    'reviewed:true AND length:[300 TO 600] AND organism_id:9606',
]


def main():
    print("Sampling proxy negative examples from UniProt Swiss-Prot...\n")
    print(f"(Will exclude {len(POSITIVE_IDS)} positive accessions to avoid overlap)\n")

    candidates: List[Dict] = []
    for query in HARD_QUERIES:
        print(f"Query: {query[:70]}{'...' if len(query) > 70 else ''}")
        results = search_uniprot(query, limit=15)
        print(f"  -> {len(results)} results")
        for r in results:
            acc = r.get("primaryAccession")
            if not acc or acc in POSITIVE_IDS:
                continue
            seq_obj = r.get("sequence", {})
            sequence = seq_obj.get("value", "")
            if not sequence or len(sequence) < 30:
                continue
            protein_desc = r.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown")
            candidates.append({
                "uniprot_id": acc,
                "name": protein_desc,
                "protein_class": "proxy negative (Swiss-Prot sample)",
                "mw_kda": seq_obj.get("molWeight", 0) / 1000.0,
                "notes": f"Sampled from query: {query[:50]}",
                "sequence": sequence.upper(),
                "sequence_length": len(sequence),
                "label": 0,
            })
        time.sleep(1)

    # Deduplicate by accession
    seen = set()
    unique = []
    for c in candidates:
        if c["uniprot_id"] in seen:
            continue
        seen.add(c["uniprot_id"])
        unique.append(c)

    # Cap at ~70 to roughly balance positives
    if len(unique) > 70:
        unique = unique[:70]

    OUTPUT_NEGATIVES.write_text(json.dumps(unique, indent=2))
    print(f"\n--- Summary ---")
    print(f"Negative examples collected: {len(unique)}")
    print(f"Saved to: {OUTPUT_NEGATIVES}")


if __name__ == "__main__":
    main()
