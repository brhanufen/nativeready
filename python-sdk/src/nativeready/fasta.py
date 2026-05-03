"""Minimal FASTA parser. No external dependencies."""
from __future__ import annotations

import os
import re
from typing import Iterable, Iterator, List, Optional

UNIPROT_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$")


def parse_fasta(path_or_text: str) -> List[dict]:
    """Parse a FASTA file or raw FASTA text into a list of records.

    Each record is a dict {'id': str, 'description': str, 'sequence': str}.
    The 'id' field is the first whitespace-delimited token after '>'.

    Parameters
    ----------
    path_or_text : str
        Either a path to a FASTA file or a string of FASTA-formatted text.

    Returns
    -------
    list of dict
    """
    text = _read_text(path_or_text)
    records = []
    current_header: Optional[str] = None
    current_seq: List[str] = []

    for line in text.splitlines():
        if line.startswith(">"):
            if current_header is not None:
                records.append(_finalize(current_header, current_seq))
            current_header = line[1:].strip()
            current_seq = []
        elif line.strip():
            current_seq.append(re.sub(r"\s+", "", line).upper())
    if current_header is not None:
        records.append(_finalize(current_header, current_seq))

    return records


def _read_text(path_or_text: str) -> str:
    """Treat the input as a path if it exists on disk; otherwise as raw text."""
    if path_or_text and len(path_or_text) < 4096 and os.path.exists(path_or_text):
        with open(path_or_text, "r", encoding="utf-8") as f:
            return f.read()
    return path_or_text


def _finalize(header: str, seq_lines: List[str]) -> dict:
    seq = "".join(seq_lines)
    # Try common header conventions:
    #   >sp|P00918|CAH2_HUMAN Carbonic anhydrase 2 OS=Homo sapiens
    #   >P00918 my description
    #   >ubiquitin custom
    rec_id = header.split()[0] if header else ""
    description = header[len(rec_id):].strip() if header else ""
    uniprot_id = _extract_uniprot(rec_id)
    return {
        "id": uniprot_id or rec_id,
        "description": description,
        "sequence": seq,
        "uniprot_id": uniprot_id,
    }


def _extract_uniprot(rec_id: str) -> Optional[str]:
    """Pull a UniProt accession out of common ID formats."""
    if not rec_id:
        return None
    # Format 1: sp|ACCESSION|NAME
    m = re.match(r"^(?:sp|tr)\|([A-Z0-9]+)\|", rec_id)
    if m and UNIPROT_RE.match(m.group(1)):
        return m.group(1)
    # Format 2: bare accession
    if UNIPROT_RE.match(rec_id):
        return rec_id
    return None


def write_fasta(records: Iterable[dict], path: str) -> None:
    """Write a list of records back out as FASTA."""
    with open(path, "w", encoding="utf-8") as f:
        for r in records:
            header = r.get("id", "unknown")
            desc = r.get("description", "")
            seq = r.get("sequence", "")
            full_header = f">{header} {desc}".strip()
            f.write(full_header + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
