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
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field

from predictor_v3 import predict as run_prediction

# Method 1 heterogeneity forward-simulator. Imported defensively: this is a
# purely additive annotation on /predict, so a problem loading it must never
# stop the API from serving model scores.
try:
    from heterogeneity import heterogeneity_report
except Exception as _exc:  # pragma: no cover - defensive
    heterogeneity_report = None
    print(f"[main] heterogeneity module unavailable, annotations disabled: {_exc}")

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


@app.on_event("startup")
def _warm_model() -> None:
    """Warm the model at container startup so the first real request is fast.

    The predictor lazy-loads ESM-2 (~2.5 GB) and the classifier bundle on the
    first prediction, which adds ~60-70s of cold-start latency to whichever
    user happens to hit /predict first. Running one throwaway prediction here
    pays that cost during boot instead. Wrapped defensively: if warm-up fails
    for any reason, the app still starts and the lazy path handles the first
    real request as before.
    """
    import time as _time

    warm_seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    try:
        _t0 = _time.time()
        run_prediction(warm_seq)
        print(f"[startup] model warm-up complete in {_time.time() - _t0:.1f}s")
    except Exception as exc:  # noqa: BLE001 - never block startup on warm-up
        print(f"[startup] model warm-up skipped ({exc!r}); lazy path will load on first request")


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

    # --- Method 1: heterogeneity forward-simulator (ADDITIVE ANNOTATION) ------
    # Computed, not predicted: arithmetic over published glycan masses and the
    # definition of resolving power. It annotates the score, it never changes
    # it, and no existing response field is touched. The field is omitted
    # entirely for proteins where heterogeneity does not apply (no sequons), so
    # their responses are byte-identical to before. Any failure here is
    # swallowed -- an annotation must never cost the caller their prediction.
    if heterogeneity_report is not None:
        try:
            het = heterogeneity_report(seq)
            if het.get("applies"):
                result["heterogeneity_risk"] = {
                    "level": het["level"],
                    "reason": het["reason"],
                    "computed": True,
                    "mode": het.get("mode"),
                    "n_sites": het.get("n_sites"),
                    "n_proteoforms": het.get("n_proteoforms"),
                    "envelope_width_da": het.get("predicted_envelope_width_da"),
                    "resolved_at_resolution": het.get("resolved_at_resolution"),
                    "resolution": het.get("resolution"),
                    "min_resolution_required": het.get("min_resolution_required"),
                    # Assumptions travel with the number, always. A computed
                    # claim is only honest if its premises ship alongside it.
                    "assumptions": het.get("assumptions", []),
                    "model_assumptions": het.get("model_assumptions", {}),
                }
        except Exception as exc:  # pragma: no cover - defensive
            print(f"[main] heterogeneity annotation skipped: {exc}")

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
    # Accept both the wire alias ("construct") and the Python field name
    # ("construct_desc"), so the JSON contract is unchanged while the attribute
    # name avoids shadowing BaseModel.construct.
    model_config = {"populate_by_name": True}

    sequence: str = Field(..., description="Sequence the user tested (will be hashed for privacy).")
    predicted_score: int = Field(..., ge=0, le=100, description="Score the model returned.")
    user_outcome: Literal["worked", "failed", "not_tested"] = Field(
        ..., description="What actually happened in the lab."
    )
    note: Optional[str] = Field(None, max_length=500, description="Optional context (e.g., conditions used).")
    # Structured experimental conditions (all optional). A native-MS outcome is
    # only meaningful with its conditions, so we capture them as their own fields
    # rather than relying on free-text `note`. None of these block a worked/failed
    # report; they enrich it when the user is willing to provide them.
    buffer: Optional[str] = Field(None, max_length=200, description="Buffer/pH used, e.g. '200 mM ammonium acetate, pH 7.5'.")
    construct_desc: Optional[str] = Field(None, alias="construct", max_length=200, description="Construct, e.g. 'full-length' or 'residues 1-256 ectodomain'.")
    expression_system: Optional[str] = Field(None, max_length=120, description="Expression system, e.g. 'E. coli', 'HEK293'.")
    instrument: Optional[str] = Field(None, max_length=200, description="Instrument, e.g. 'Waters Synapt G2-S'.")
    resolution: Optional[str] = Field(None, max_length=40, description="Resolving power, e.g. '30000'.")
    failure_mode: Optional[str] = Field(None, max_length=200, description="Observed failure mode (only when outcome=failed).")
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

    # Data-quality gate. A real outcome (worked or failed) that becomes a
    # training label has to carry the conditions that make it interpretable:
    # without buffer/instrument the result cannot be tied to the run, and a
    # failure especially must be told apart from a bad-conditions failure. Both
    # worked and failed therefore require at least a buffer or an instrument, and
    # a failure additionally requires a failure mode. "not_tested" is exempt (the
    # experiment has not run yet). The minimum is deliberately small so anyone who
    # actually ran the experiment can meet it; the dropdowns make it a few taps.
    _buf = (req.buffer or "").strip()
    _instr = (req.instrument or "").strip()
    _fmode = (req.failure_mode or "").strip()
    if req.user_outcome in ("worked", "failed") and not (_buf or _instr):
        raise HTTPException(
            status_code=422,
            detail=(
                f"A '{req.user_outcome}' outcome needs at least the buffer or the "
                "instrument used, so the result can be tied to the run (a failure "
                "especially must be told apart from a conditions problem). Add one "
                "and resubmit, or choose 'not tested' if you have not run it yet."
            ),
        )
    if req.user_outcome == "failed" and not _fmode:
        raise HTTPException(
            status_code=422,
            detail=(
                "A 'failed' outcome needs a failure mode (what went wrong, e.g. "
                "'no ionization', 'unresolved heterogeneity', 'salt adducts'). "
                "That is the part that actually trains the model. Add it and resubmit."
            ),
        )

    # Hash sequence for privacy (don't store raw user proteins)
    sequence_hash = hashlib.sha256(seq.encode("utf-8")).hexdigest()[:16]

    record = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sequence_hash": sequence_hash,
        "sequence_length": len(seq),
        "predicted_score": req.predicted_score,
        "user_outcome": req.user_outcome,
        "note": (req.note or "").strip()[:500] or None,
        "buffer": (req.buffer or "").strip()[:200] or None,
        "construct": (req.construct_desc or "").strip()[:200] or None,
        "expression_system": (req.expression_system or "").strip()[:120] or None,
        "instrument": (req.instrument or "").strip()[:200] or None,
        "resolution": (req.resolution or "").strip()[:40] or None,
        "failure_mode": (req.failure_mode or "").strip()[:200] or None,
        "model_version": req.model_version,
        "email_for_followup": (req.email_for_followup or "").strip()[:200] or None,
    }
    # Dedup guard: a double submission (double-click, SDK retry) must not create
    # a duplicate record. The server timestamp differs per request, so the guard
    # keys on the record's *content* — every field except the timestamp and the
    # optional follow-up email. An identical prior record is treated as a resend
    # and skipped; a genuinely different report (a new outcome or new conditions
    # for the same sequence) is still recorded. Safe by construction: it only
    # suppresses exact duplicates, never a distinct experimental outcome.
    def _dedup_signature(rec: Dict[str, Any]) -> tuple:
        return tuple(
            rec.get(k)
            for k in (
                "sequence_hash", "user_outcome", "predicted_score", "note",
                "buffer", "construct", "expression_system", "instrument",
                "resolution", "failure_mode", "model_version",
            )
        )

    incoming_sig = _dedup_signature(record)
    for _prior in _read_jsonl(FEEDBACK_LOG):
        if _dedup_signature(_prior) == incoming_sig:
            return {
                "status": "duplicate_ignored",
                "message": "This outcome was already recorded; duplicate skipped.",
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


class PruneRequest(BaseModel):
    """Criteria for removing feedback records. A record is a match only if every
    provided field equals the record's stored value (pass null to require null).
    The prune refuses unless the number of matches equals `expect`, so a broad or
    mistaken filter can never quietly delete the wrong rows."""
    match: Dict[str, Any] = Field(..., description="Field->value criteria a record must match exactly.")
    expect: int = Field(1, ge=1, le=50, description="Exact number of matches required, or the prune refuses.")


@app.post("/admin/feedback/prune")
def admin_feedback_prune(req: PruneRequest, token: str = "") -> Dict[str, Any]:
    """Remove specific feedback records by exact-match criteria (admin only).

    Safe by construction: it snapshots the current log to a timestamped backup on
    the same volume before writing, matches only records where *every* provided
    field equals the stored value, and refuses (writing nothing) unless the match
    count equals `expect`. Requires ?token=... matching NATIVEREADY_ADMIN_TOKEN.
    """
    if not ADMIN_TOKEN:
        raise HTTPException(status_code=503, detail="admin_token_not_configured_on_server")
    if token.strip() != ADMIN_TOKEN:
        raise HTTPException(status_code=403, detail="forbidden")

    records = list(_read_jsonl(FEEDBACK_LOG))
    crit = req.match

    def _matches(rec: Dict[str, Any]) -> bool:
        return all(rec.get(k) == v for k, v in crit.items())

    matched = [r for r in records if _matches(r)]
    if len(matched) != req.expect:
        raise HTTPException(
            status_code=409,
            detail=(
                f"refused: matched {len(matched)} record(s) but expect={req.expect}. "
                "Nothing was deleted. Refine the match criteria."
            ),
        )

    # Snapshot before rewriting, on the same durable volume.
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    backup_path = FEEDBACK_LOG.parent / f"feedback-prune-backup-{stamp}.jsonl"
    try:
        FEEDBACK_LOG.parent.mkdir(parents=True, exist_ok=True)
        backup_path.write_text(
            "".join(json.dumps(r) + "\n" for r in records), encoding="utf-8"
        )
        kept = [r for r in records if not _matches(r)]
        with open(FEEDBACK_LOG, "w", encoding="utf-8") as f:
            for r in kept:
                f.write(json.dumps(r) + "\n")
    except OSError as exc:
        raise HTTPException(status_code=500, detail=f"prune write failed: {exc}")

    return {
        "status": "pruned",
        "removed": len(matched),
        "remaining": len(records) - len(matched),
        "backup": str(backup_path),
    }


@app.get("/admin/export")
def admin_export(token: str = "") -> PlainTextResponse:
    """Full feedback-log export (JSONL) for off-Railway backup.

    Requires ?token=... matching NATIVEREADY_ADMIN_TOKEN (same auth as
    /admin/stats). Returns the raw feedback.jsonl so the entire flywheel dataset
    can be downloaded over HTTP with no SSH or volume access. Pair with a cron/
    scheduled `curl` to keep an off-Railway copy of the moat.
    """
    if not ADMIN_TOKEN:
        raise HTTPException(status_code=503, detail="admin_token_not_configured_on_server")
    if token.strip() != ADMIN_TOKEN:
        raise HTTPException(status_code=403, detail="forbidden")
    try:
        content = FEEDBACK_LOG.read_text(encoding="utf-8") if FEEDBACK_LOG.exists() else ""
    except OSError as exc:
        raise HTTPException(status_code=500, detail=f"feedback log read failed: {exc}")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    return PlainTextResponse(
        content,
        media_type="application/x-ndjson",
        headers={"Content-Disposition": f'attachment; filename="feedback-{stamp}.jsonl"'},
    )
