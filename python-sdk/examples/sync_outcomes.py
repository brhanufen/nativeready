"""Sync a lab's experimental outcomes back to NativeReady, once and forever.

This is the design-partner integration that makes the flywheel run without any
per-experiment human effort. You run it on a schedule (a nightly or weekly cron
job). It reads the lab's own results table, finds rows that have not been logged
yet, and reports each real outcome to NativeReady via the SDK.

The idea
--------
A pipeline calls ``predict()`` at design time. Nobody's code naturally calls
``report_outcome()`` weeks later when the wet-lab result comes back. This script
is that missing bridge: point it at wherever the lab records what happened, map a
few columns, and the loop closes itself from then on.

What you adapt for each partner
-------------------------------
1. ``load_outcome_rows()`` - read from wherever the lab keeps results
   (a CSV export, a LIMS/ELN query, a Google Sheet, a database). Only this
   function is lab-specific; everything below is generic.
2. ``COLUMN_MAP`` - name the columns in their table.
3. Nothing else. The dedup + reporting logic is reusable.

Idempotent by design
---------------------
The script keeps a small local file of already-synced row ids, so re-running it
never double-logs. The API also has its own content-signature dedup guard as a
second line of defense, so a duplicate is harmless either way.

Usage
-----
    export NATIVEREADY_URL="https://nativeready-production.up.railway.app"
    python sync_outcomes.py --results lab_results.csv
    python sync_outcomes.py --results lab_results.csv --dry-run   # preview only

Then wire it into cron, e.g. every night at 2am:
    0 2 * * *  cd /path/to && /usr/bin/python3 sync_outcomes.py --results lab_results.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List

from nativeready import Client
from nativeready.client import NativeReadyError

# --------------------------------------------------------------------------
# 1. Map the lab's column names to NativeReady's fields.
#    Left  = the key this script uses.
#    Right = the exact column header in the lab's results table.
#    Only "sequence", "predicted_score", and "outcome" are required; the rest
#    are optional but are what make a "failed" result scientifically useful.
# --------------------------------------------------------------------------
COLUMN_MAP: Dict[str, str] = {
    "row_id": "experiment_id",        # a stable unique id per experiment (for dedup)
    "sequence": "protein_sequence",
    "predicted_score": "nativeready_score",
    "outcome": "result",              # must map to worked / failed / not_tested
    "buffer": "buffer",
    "construct": "construct",
    "expression_system": "expression",
    "instrument": "ms_instrument",
    "resolution": "resolving_power",
    "failure_mode": "failure_reason",
    "model_version": "model_version",
}

# How the lab's own result words map to NativeReady's three outcomes.
# Adapt the left-hand strings to whatever the lab actually writes.
OUTCOME_MAP: Dict[str, str] = {
    "success": "worked", "worked": "worked", "pass": "worked", "yes": "worked",
    "fail": "failed", "failed": "failed", "no signal": "failed", "no": "failed",
    "pending": "not_tested", "planned": "not_tested", "not run": "not_tested",
}

SYNCED_STATE = Path(".nativeready_synced_ids.json")


def load_outcome_rows(results_path: str) -> List[Dict[str, str]]:
    """Read the lab's results table into a list of dict rows.

    ADAPT THIS PER PARTNER. The default reads a CSV. To read from a LIMS/ELN or
    a database instead, replace the body with your query and return a list of
    dicts keyed by the lab's own column names.
    """
    rows: List[Dict[str, str]] = []
    with open(results_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def _load_synced() -> set:
    if SYNCED_STATE.exists():
        try:
            return set(json.loads(SYNCED_STATE.read_text()))
        except (json.JSONDecodeError, OSError):
            return set()
    return set()


def _save_synced(ids: Iterable[str]) -> None:
    SYNCED_STATE.write_text(json.dumps(sorted(ids)))


def _get(row: Dict[str, str], key: str) -> str | None:
    col = COLUMN_MAP.get(key)
    if not col:
        return None
    val = (row.get(col) or "").strip()
    return val or None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Sync lab outcomes to NativeReady.")
    ap.add_argument("--results", required=True, help="Path to the lab's results CSV")
    ap.add_argument("--url", default=os.environ.get(
        "NATIVEREADY_URL", "https://nativeready-production.up.railway.app"))
    ap.add_argument("--dry-run", action="store_true",
                    help="Print what would be sent, send nothing")
    args = ap.parse_args(argv)

    client = Client(base_url=args.url)
    synced = _load_synced()
    rows = load_outcome_rows(args.results)

    sent, skipped, errors = 0, 0, 0
    newly_synced = set(synced)

    for row in rows:
        row_id = _get(row, "row_id")
        if row_id and row_id in synced:
            skipped += 1
            continue

        seq = _get(row, "sequence")
        raw_outcome = (_get(row, "outcome") or "").lower()
        outcome = OUTCOME_MAP.get(raw_outcome)
        score = _get(row, "predicted_score")

        # A row needs a sequence, a recognized outcome, and a score to be usable.
        if not seq or not outcome or not score:
            skipped += 1
            continue

        # Data-quality gate, mirroring the API. A worked/failed outcome must carry
        # the minimum conditions that make it interpretable; a failed one also
        # needs a failure mode. Rows that cannot meet this are skipped (with a
        # reason) rather than sent and rejected. "not_tested" rows are exempt.
        _buf = _get(row, "buffer")
        _instr = _get(row, "instrument")
        _fmode = _get(row, "failure_mode")
        if outcome in ("worked", "failed") and not (_buf or _instr):
            print(f"  skip {row_id or '?'}: {outcome} needs buffer or instrument",
                  file=sys.stderr)
            skipped += 1
            continue
        if outcome == "failed" and not _fmode:
            print(f"  skip {row_id or '?'}: failed needs a failure_mode",
                  file=sys.stderr)
            skipped += 1
            continue

        payload = dict(
            sequence=seq,
            predicted_score=int(float(score)),
            outcome=outcome,
            buffer=_get(row, "buffer"),
            construct=_get(row, "construct"),
            expression_system=_get(row, "expression_system"),
            instrument=_get(row, "instrument"),
            resolution=_get(row, "resolution"),
            failure_mode=_get(row, "failure_mode"),
            model_version=_get(row, "model_version"),
        )

        if args.dry_run:
            print("WOULD SEND:", json.dumps({k: v for k, v in payload.items() if v}))
            sent += 1
            if row_id:
                newly_synced.add(row_id)
            continue

        try:
            client.report_outcome(**payload)
            sent += 1
            if row_id:
                newly_synced.add(row_id)
        except (NativeReadyError, ValueError) as e:
            print(f"  row {row_id or '?'}: {e}", file=sys.stderr)
            errors += 1

    if not args.dry_run:
        _save_synced(newly_synced)

    print(f"\nDone. sent={sent} skipped={skipped} errors={errors} "
          f"(dry-run)" if args.dry_run else
          f"\nDone. sent={sent} skipped={skipped} errors={errors}")
    return 0 if errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
