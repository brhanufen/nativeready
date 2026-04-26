"""NativeReady FastAPI backend.

Endpoints:
  GET  /         -> health check
  POST /predict  -> native-MS suitability prediction for a protein sequence
  GET  /docs     -> Swagger UI (FastAPI default)
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Literal, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

from predictor_v2 import predict as run_prediction

_HERE = Path(__file__).resolve().parent
FEEDBACK_LOG = Path(
    os.environ.get(
        "NATIVEREADY_FEEDBACK_LOG",
        str(_HERE.parent / "data" / "feedback.jsonl"),
    )
)

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYXBZ")
MIN_LEN = 10
MAX_LEN = 5000

app = FastAPI(
    title="NativeReady",
    description=(
        "Predicts whether a protein sequence is likely to give usable native "
        "mass spectrometry data."
    ),
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # tighten in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class PredictRequest(BaseModel):
    sequence: str = Field(
        ...,
        description=(
            "Protein sequence as a FASTA-style string. A leading '>' header "
            "line is allowed and will be stripped."
        ),
    )


def _clean_sequence(raw: str) -> str:
    """Strip FASTA header (if any), whitespace, and uppercase the sequence."""
    if not isinstance(raw, str):
        raise HTTPException(
            status_code=400, detail="Sequence must be a string."
        )
    text = raw.strip()
    if text.startswith(">"):
        # Drop the first header line
        parts = text.split("\n", 1)
        text = parts[1] if len(parts) > 1 else ""
    # Remove all whitespace (spaces, newlines, tabs)
    cleaned = re.sub(r"\s+", "", text).upper()
    return cleaned


def _validate_sequence(seq: str) -> None:
    if len(seq) == 0:
        raise HTTPException(
            status_code=400,
            detail="Empty sequence — please provide a protein sequence.",
        )
    if len(seq) < MIN_LEN:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Sequence too short to analyze "
                f"({len(seq)} residues; minimum is {MIN_LEN})."
            ),
        )
    if len(seq) > MAX_LEN:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Sequence too long for current model "
                f"({len(seq)} residues; maximum is {MAX_LEN})."
            ),
        )
    bad = sorted({c for c in seq if c not in ALLOWED_AA})
    if bad:
        raise HTTPException(
            status_code=400,
            detail=(
                "Sequence contains invalid character(s): "
                f"{', '.join(bad)}. Allowed are the 20 standard amino acids "
                "plus X, B, Z."
            ),
        )


@app.get("/")
def health() -> Dict[str, str]:
    return {"status": "ok", "service": "nativeready"}


@app.post("/predict")
def predict_endpoint(req: PredictRequest) -> Dict[str, Any]:
    seq = _clean_sequence(req.sequence)
    _validate_sequence(seq)
    try:
        return run_prediction(seq)
    except HTTPException:
        raise
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(
            status_code=500, detail=f"Prediction failed: {exc}"
        )


# --------------------------------------------------------------------------
# Feedback endpoint — collects real-world experimental outcomes from users.
# This is the data flywheel that lets the model improve over time.
# --------------------------------------------------------------------------

class FeedbackRequest(BaseModel):
    sequence: str = Field(..., description="Sequence the user tested (will be hashed for privacy).")
    predicted_score: int = Field(..., ge=0, le=100, description="Score the model returned.")
    user_outcome: Literal["worked", "failed", "not_tested"] = Field(
        ..., description="What actually happened in the lab."
    )
    note: Optional[str] = Field(None, max_length=500, description="Optional context (e.g., conditions used).")
    model_version: Optional[str] = Field(None, max_length=64)


@app.post("/feedback")
def feedback_endpoint(req: FeedbackRequest, request: Request) -> Dict[str, Any]:
    """Append a real-world outcome to the feedback log.

    Stores: timestamp, hashed sequence (privacy-preserving), predicted score,
    user-reported outcome, optional note. No personal data, no IP storage,
    no email collected.
    """
    seq = _clean_sequence(req.sequence)
    _validate_sequence(seq)

    # Hash sequence for privacy (don't store raw user proteins)
    sequence_hash = hashlib.sha256(seq.encode("utf-8")).hexdigest()[:16]

    record = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sequence_hash": sequence_hash,
        "sequence_length": len(seq),
        "predicted_score": req.predicted_score,
        "user_outcome": req.user_outcome,
        "note": (req.note or "").strip()[:500] or None,
        "model_version": req.model_version,
    }
    # Append-only JSONL log
    try:
        FEEDBACK_LOG.parent.mkdir(parents=True, exist_ok=True)
        with open(FEEDBACK_LOG, "a", encoding="utf-8") as f:
            f.write(json.dumps(record) + "\n")
    except OSError as exc:
        raise HTTPException(status_code=500, detail=f"Feedback log write failed: {exc}")

    return {
        "status": "received",
        "message": "Thanks. Your input helps the model learn from real experiments.",
    }


@app.get("/feedback/stats")
def feedback_stats() -> Dict[str, Any]:
    """Public summary stats — no individual records exposed."""
    if not FEEDBACK_LOG.exists():
        return {"total_feedback": 0, "outcomes": {}}
    counts = {"worked": 0, "failed": 0, "not_tested": 0}
    total = 0
    try:
        with open(FEEDBACK_LOG, "r", encoding="utf-8") as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    outcome = rec.get("user_outcome")
                    if outcome in counts:
                        counts[outcome] += 1
                    total += 1
                except json.JSONDecodeError:
                    continue
    except OSError:
        pass
    return {"total_feedback": total, "outcomes": counts}
