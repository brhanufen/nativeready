"""Command-line interface for the NativeReady SDK.

Usage examples:
    nativeready predict --sequence "MQIFVKTL..."
    nativeready predict --uniprot P00918
    nativeready predict --fasta my_proteins.fasta --output results.csv
    nativeready predict --fasta my_proteins.fasta --output results.json
    nativeready health
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import List

from . import __version__
from .client import Client, PredictionResult, NativeReadyError


def _format_table(results: List[PredictionResult], detailed: bool = False) -> str:
    """Pretty-print a small results table to stdout."""
    if not results:
        return "(no results)"
    rows = []
    rows.append(["#", "id", "len", "score", "label", "ci", "model", "ood"])
    for i, r in enumerate(results, 1):
        rows.append([
            str(i),
            (r.uniprot_id or "")[:20],
            str(r.sequence_length or ""),
            str(r.score),
            r.label[:12],
            f"{r.confidence_lower}-{r.confidence_upper}",
            r.model_version[:18],
            "Y" if r.is_ood else "",
        ])
    widths = [max(len(row[c]) for row in rows) for c in range(len(rows[0]))]
    lines = []
    for ri, row in enumerate(rows):
        lines.append("  ".join(cell.ljust(widths[ci]) for ci, cell in enumerate(row)))
        if ri == 0:
            lines.append("  ".join("-" * w for w in widths))
    return "\n".join(lines)


def _save_results(results: List[PredictionResult], path: str) -> None:
    """Write results to CSV or JSON based on file extension."""
    p = Path(path)
    if p.suffix.lower() == ".json":
        with p.open("w", encoding="utf-8") as f:
            json.dump([r.as_dict() for r in results], f, indent=2)
    else:
        # default to CSV
        with p.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow([
                "uniprot_id", "sequence_length", "score", "label",
                "ci_lower", "ci_upper", "is_ood", "model_version",
            ])
            for r in results:
                w.writerow([
                    r.uniprot_id or "",
                    r.sequence_length or "",
                    r.score,
                    r.label,
                    r.confidence_lower,
                    r.confidence_upper,
                    "Y" if r.is_ood else "",
                    r.model_version,
                ])


def cmd_predict(args: argparse.Namespace) -> int:
    client = Client(
        base_url=args.url,
        timeout=args.timeout,
    )

    if args.sequence:
        result = client.predict(args.sequence)
        results = [result]
    elif args.uniprot:
        accessions = [a.strip() for a in args.uniprot.split(",") if a.strip()]
        if len(accessions) == 1:
            results = [client.predict_uniprot(accessions[0])]
        else:
            print(f"Predicting {len(accessions)} UniProt accessions...")
            results = []
            for acc in accessions:
                try:
                    results.append(client.predict_uniprot(acc))
                except NativeReadyError as e:
                    print(f"  {acc}: failed ({e})")
    elif args.fasta:
        print(f"Reading FASTA: {args.fasta}")
        results = client.predict_fasta(args.fasta, progress=not args.quiet)
    else:
        print("Error: must provide --sequence, --uniprot, or --fasta", file=sys.stderr)
        return 2

    print()
    print(_format_table(results))
    print()

    if args.output:
        _save_results(results, args.output)
        print(f"Saved to {args.output}")

    return 0


def cmd_health(args: argparse.Namespace) -> int:
    client = Client(base_url=args.url, timeout=args.timeout)
    try:
        h = client.health()
        print(json.dumps(h, indent=2))
        return 0
    except NativeReadyError as e:
        print(f"Health check failed: {e}", file=sys.stderr)
        return 1


def cmd_stats(args: argparse.Namespace) -> int:
    client = Client(base_url=args.url, timeout=args.timeout)
    try:
        s = client.feedback_stats()
        print(json.dumps(s, indent=2))
        return 0
    except NativeReadyError as e:
        print(f"Stats fetch failed: {e}", file=sys.stderr)
        return 1


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="nativeready",
        description="Predict native mass spectrometry suitability from protein sequences.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    sub = parser.add_subparsers(dest="command", required=True)

    # predict
    p = sub.add_parser("predict", help="Predict suitability for one or many sequences")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--sequence", "-s", help="Raw protein sequence (FASTA header optional)")
    src.add_argument("--uniprot", "-u", help="UniProt accession (or comma-separated list)")
    src.add_argument("--fasta", "-f", help="Path to a FASTA file")
    p.add_argument("--output", "-o", help="Save results to CSV or JSON (extension determines format)")
    p.add_argument("--url", default="https://nativeready-production.up.railway.app",
                   help="Override the API base URL")
    p.add_argument("--timeout", type=int, default=120, help="Request timeout in seconds")
    p.add_argument("--quiet", action="store_true", help="Suppress progress bar")
    p.set_defaults(func=cmd_predict)

    # health
    p2 = sub.add_parser("health", help="Health check the API")
    p2.add_argument("--url", default="https://nativeready-production.up.railway.app")
    p2.add_argument("--timeout", type=int, default=15)
    p2.set_defaults(func=cmd_health)

    # stats
    p3 = sub.add_parser("stats", help="Fetch public aggregate feedback statistics")
    p3.add_argument("--url", default="https://nativeready-production.up.railway.app")
    p3.add_argument("--timeout", type=int, default=15)
    p3.set_defaults(func=cmd_stats)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
