"""
Fetch real protein sequences from the UniProt REST API.

UniProt is the canonical public database of protein sequences. All sequences
fetched here are real, peer-reviewed entries. No synthetic data.

API documentation: https://www.uniprot.org/help/api_queries
Endpoint format: https://rest.uniprot.org/uniprotkb/{accession}.fasta
"""
import json
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional, Tuple

from positives_curated import POSITIVE_EXAMPLES

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{}.json"
DATA_DIR = Path(__file__).parent
OUTPUT_FASTA = DATA_DIR / "positives.fasta"
OUTPUT_META = DATA_DIR / "positives_with_sequences.json"


def fetch_fasta(accession: str, retries: int = 3) -> Optional[str]:
    """Fetch FASTA sequence for a UniProt accession."""
    url = UNIPROT_URL.format(accession)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "NativeReady/0.1"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                if resp.status == 200:
                    return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f"  [404] {accession} not found")
                return None
            print(f"  [HTTP {e.code}] {accession} attempt {attempt+1}")
        except (urllib.error.URLError, TimeoutError) as e:
            print(f"  [error] {accession} attempt {attempt+1}: {e}")
        time.sleep(1.5 * (attempt + 1))
    return None


def parse_fasta(fasta: str) -> Tuple[str, str]:
    """Parse a single FASTA record. Returns (header, sequence)."""
    lines = fasta.strip().split("\n")
    header = lines[0].lstrip(">")
    sequence = "".join(lines[1:]).replace(" ", "").upper()
    return header, sequence


def main():
    print(f"Fetching {len(POSITIVE_EXAMPLES)} sequences from UniProt REST API...\n")
    enriched = []
    fasta_lines = []
    successes = 0
    failures = []

    for i, entry in enumerate(POSITIVE_EXAMPLES, 1):
        accession = entry["uniprot_id"]
        print(f"[{i:>2}/{len(POSITIVE_EXAMPLES)}] Fetching {accession} ({entry['name']})...")
        fasta = fetch_fasta(accession)
        if fasta is None:
            failures.append(accession)
            continue
        header, sequence = parse_fasta(fasta)
        enriched.append({
            **entry,
            "sequence": sequence,
            "sequence_length": len(sequence),
            "fasta_header": header,
            "label": 1,  # positive example: known native MS target
        })
        fasta_lines.append(fasta.strip())
        successes += 1
        time.sleep(0.3)  # be polite to UniProt

    OUTPUT_FASTA.write_text("\n\n".join(fasta_lines) + "\n")
    OUTPUT_META.write_text(json.dumps(enriched, indent=2))

    print(f"\n--- Summary ---")
    print(f"Successfully fetched: {successes}/{len(POSITIVE_EXAMPLES)}")
    if failures:
        print(f"Failed: {len(failures)} ({', '.join(failures)})")
    print(f"FASTA saved to: {OUTPUT_FASTA}")
    print(f"Metadata saved to: {OUTPUT_META}")


if __name__ == "__main__":
    main()
