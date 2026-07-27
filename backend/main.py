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

from predictor_v3 import predict as run_prediction

_HERE = Path(__file__).resolve().parent
FEEDBACK_LOG = Path(
    os.environ.get(
        "NATIVEREADY_FEEDBACK_LOG",
        str(_HERE.parent / "data" / "feedback.jsonl"),
    )
)
PREDICTION_LOG = Path(
    os.environ.get(
        "NATIVEREADY_PREDICTION_LOG",
        str(_HERE.parent / "data" / "predictions.jsonl"),
    )
)
# Shared-secret token for /admin/stats. Set NATIVEREADY_ADMIN_TOKEN in Railway env.
# Strip whitespace defensively — env var systems sometimes append newlines.
ADMIN_TOKEN = os.environ.get("NATIVEREADY_ADMIN_TOKEN", "").strip()

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYXBZ")
MIN_LEN = 10
MAX_LEN = 5000

app = FastAPI(
    title="NativeReady",
    description=(
        "Predicts whether a protein sequence is likely to give usable native "
        "mass spectrometry data."
    ),
    version="0.3.0",
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


def _visitor_id(request: Optional[Request]) -> str:
    """Privacy-preserving distinct-user signal.

    Hashes client IP + User-Agent into a short, non-reversible token so we can
    count *distinct* users without storing any personal data. Respects the
    X-Forwarded-For header set by Railway's proxy. An SDK caller may also send
    an explicit X-NativeReady-Client header to identify itself consistently.
    """
    if request is None:
        return "unknown"
    # Explicit client id (SDK can set this) takes precedence.
    client_hdr = request.headers.get("x-nativeready-client", "").strip()
    if client_hdr:
        basis = "client:" + client_hdr
    else:
        xff = request.headers.get("x-forwarded-for", "")
        ip = xff.split(",")[0].strip() if xff else (
            request.client.host if request.client else "?"
        )
        ua = request.headers.get("user-agent", "")
        basis = f"{ip}|{ua}"
    return hashlib.sha256(basis.encode("utf-8")).hexdigest()[:16]


@app.post("/predict")
def predict_endpoint(req: PredictRequest, request: Request = None) -> Dict[str, Any]:
    seq = _clean_sequence(req.sequence)
    _validate_sequence(seq)
    try:
        result = run_prediction(seq)
    except HTTPException:
        raise
    except Exception as exc:  # pragma: no cover - defensive
        raise HTTPException(
            status_code=500, detail=f"Prediction failed: {exc}"
        )
    # Append-only prediction log for private usage dashboard.
    # Stores hashed sequence + hashed visitor id only — no raw protein data,
    # no IPs, no personal data persisted.
    try:
        PREDICTION_LOG.parent.mkdir(parents=True, exist_ok=True)
        seq_hash = hashlib.sha256(seq.encode("utf-8")).hexdigest()[:16]
        log_record = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "sequence_hash": seq_hash,
            "visitor_id": _visitor_id(request),
            "sequence_length": len(seq),
            "suitability_score": result.get("suitability_score"),
            "suitability_label": result.get("suitability_label"),
            "model_version": result.get("model_version"),
            "ood": bool(result.get("ood", False)),
        }
        with open(PREDICTION_LOG, "a", encoding="utf-8") as f:
            f.write(json.dumps(log_record) + "\n")
    except OSError:
        pass  # logging never blocks the response
    return result


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
    email_for_followup: Optional[str] = Field(
        None,
        max_length=200,
        description="Optional email — opts user in to a single follow-up message in 2-4 weeks asking how the experiment went.",
    )


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
        "email_for_followup": (req.email_for_followup or "").strip()[:200] or None,
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


# --------------------------------------------------------------------------
# Admin stats — private dashboard for the founder. Token-protected.
# Set NATIVEREADY_ADMIN_TOKEN in the Railway environment to enable.
# --------------------------------------------------------------------------

def _read_jsonl(path: Path):
    if not path.exists():
        return
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                try:
                    yield json.loads(line)
                except json.JSONDecodeError:
                    continue
    except OSError:
        return


@app.get("/admin/stats")
def admin_stats(token: str = "") -> Dict[str, Any]:
    """Private usage dashboard. Requires ?token=... matching NATIVEREADY_ADMIN_TOKEN."""
    if not ADMIN_TOKEN:
        # Diagnostic: env var not loaded by container at all.
        raise HTTPException(status_code=503, detail="admin_token_not_configured_on_server")
    if token.strip() != ADMIN_TOKEN:
        raise HTTPException(status_code=403, detail="forbidden")

    # ---- Predictions
    pred_total = 0
    per_day: Dict[str, int] = {}
    score_buckets = {"0-34": 0, "35-49": 0, "50-64": 0, "65-79": 0, "80-100": 0}
    ood_count = 0
    model_versions: Dict[str, int] = {}
    visitors_all: set = set()          # distinct users (hashed), all-time
    visitors_30d: set = set()          # distinct users in the last 30 days
    users_per_day: Dict[str, set] = {}  # distinct users seen each day
    sequences_all: set = set()          # distinct sequences ever scored
    # 30-day cutoff (YYYY-MM-DD string compare is fine for ISO dates)
    from datetime import timedelta
    cutoff = (datetime.now(timezone.utc) - timedelta(days=30)).strftime("%Y-%m-%d")
    for rec in _read_jsonl(PREDICTION_LOG):
        pred_total += 1
        ts = rec.get("timestamp", "")[:10]
        if ts:
            per_day[ts] = per_day.get(ts, 0) + 1
        vid = rec.get("visitor_id")
        if vid:
            visitors_all.add(vid)
            users_per_day.setdefault(ts, set()).add(vid)
            if ts and ts >= cutoff:
                visitors_30d.add(vid)
        sh = rec.get("sequence_hash")
        if sh:
            sequences_all.add(sh)
        score = rec.get("suitability_score")
        if isinstance(score, int):
            if score < 35: score_buckets["0-34"] += 1
            elif score < 50: score_buckets["35-49"] += 1
            elif score < 65: score_buckets["50-64"] += 1
            elif score < 80: score_buckets["65-79"] += 1
            else: score_buckets["80-100"] += 1
        if rec.get("ood"):
            ood_count += 1
        mv = rec.get("model_version") or "unknown"
        model_versions[mv] = model_versions.get(mv, 0) + 1

    # Last 30 days of activity, oldest -> newest (predictions + distinct users)
    sorted_days = sorted(per_day.items())[-30:]
    users_by_day = [(d, len(users_per_day.get(d, ()))) for d, _ in sorted_days]

    # ---- Feedback
    fb_total = 0
    fb_outcomes = {"worked": 0, "failed": 0, "not_tested": 0}
    recent_notes = []
    emails_collected = 0
    for rec in _read_jsonl(FEEDBACK_LOG):
        fb_total += 1
        out = rec.get("user_outcome")
        if out in fb_outcomes:
            fb_outcomes[out] += 1
        if rec.get("email_for_followup"):
            emails_collected += 1
        note = rec.get("note")
        if note:
            recent_notes.append({
                "timestamp": rec.get("timestamp"),
                "outcome": out,
                "predicted_score": rec.get("predicted_score"),
                "note": note,
            })

    return {
        "predictions": {
            "total": pred_total,
            "distinct_users_all_time": len(visitors_all),
            "distinct_users_last_30d": len(visitors_30d),
            "distinct_sequences_all_time": len(sequences_all),
            "per_day_last_30": sorted_days,
            "distinct_users_per_day_last_30": users_by_day,
            "score_distribution": score_buckets,
            "ood_count": ood_count,
            "model_versions": model_versions,
        },
        "feedback": {
            "total": fb_total,
            "outcomes": fb_outcomes,
            "emails_collected_for_followup": emails_collected,
            "recent_notes": recent_notes[-20:],
        },
    }
