"""HTTP client for the NativeReady prediction API."""
from __future__ import annotations

import re
import time
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import requests

DEFAULT_API_URL = "https://nativeready-production.up.railway.app"
DEFAULT_TIMEOUT = 120  # seconds; first request is slow due to ESM-2 cold start
DEFAULT_RETRIES = 2

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{}.fasta"


@dataclass
class PredictionResult:
    """Structured result of a single prediction."""

    score: int  # 0-100 calibrated suitability
    label: str  # "Excellent", "Good", "Fair", "Poor", "Unsuitable"
    confidence_lower: int
    confidence_upper: int
    model_version: str
    risk_factors: List[Dict[str, str]]
    recommendations: List[str]
    is_ood: bool = False
    ood_message: Optional[str] = None
    uniprot_id: Optional[str] = None  # if predicted from UniProt accession or FASTA header
    sequence_length: Optional[int] = None
    raw_response: Dict[str, Any] = field(default_factory=dict, repr=False)

    @classmethod
    def from_api(
        cls,
        response: Dict[str, Any],
        uniprot_id: Optional[str] = None,
        sequence_length: Optional[int] = None,
    ) -> "PredictionResult":
        ci = response.get("confidence_interval", {})
        warning = response.get("warning")
        return cls(
            score=int(response.get("suitability_score", 0)),
            label=str(response.get("suitability_label", "")),
            confidence_lower=int(ci.get("lower", 0)),
            confidence_upper=int(ci.get("upper", 0)),
            model_version=str(response.get("model_version", "")),
            risk_factors=list(response.get("risk_factors", [])),
            recommendations=list(response.get("recommendations", [])),
            is_ood=warning == "out_of_distribution",
            ood_message=response.get("warning_message"),
            uniprot_id=uniprot_id,
            sequence_length=sequence_length,
            raw_response=dict(response),
        )

    def as_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d.pop("raw_response", None)
        return d

    def __repr__(self) -> str:
        ood = " (OOD)" if self.is_ood else ""
        return (
            f"PredictionResult(score={self.score}, label='{self.label}', "
            f"ci=[{self.confidence_lower}-{self.confidence_upper}], "
            f"model='{self.model_version}'{ood})"
        )


class NativeReadyError(Exception):
    """Base exception for the SDK."""


class APIError(NativeReadyError):
    """Raised when the NativeReady API returns an error."""


class SequenceError(NativeReadyError):
    """Raised when the input sequence is invalid."""


def _clean_sequence(raw: str) -> str:
    """Strip optional FASTA header and whitespace; return uppercase sequence."""
    if not isinstance(raw, str):
        raise SequenceError("Sequence must be a string")
    text = raw.strip()
    if text.startswith(">"):
        parts = text.split("\n", 1)
        text = parts[1] if len(parts) > 1 else ""
    cleaned = re.sub(r"\s+", "", text).upper()
    return cleaned


def _extract_uniprot_from_header(raw: str) -> Optional[str]:
    """Try to extract a UniProt accession from a FASTA header line."""
    if not raw or not raw.lstrip().startswith(">"):
        return None
    header = raw.lstrip().split("\n", 1)[0]
    m = re.search(r">(?:sp|tr)\|([A-Z0-9]+)\|", header)
    if m:
        return m.group(1)
    m = re.search(r">([A-Z0-9]+)\b", header)
    if m and re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}", m.group(1)):
        return m.group(1)
    return None


class Client:
    """HTTP client for the NativeReady API.

    Parameters
    ----------
    base_url : str, optional
        Override the default production URL. Useful for local development
        or self-hosted deployments.
    timeout : int, optional
        Request timeout in seconds. The first request triggers an ESM-2
        cold-start on the server (~30 s); subsequent requests are ~3-5 s.
        Default is 120 to comfortably cover the cold start.
    retries : int, optional
        Number of automatic retries on transient errors. Default 2.
    session : requests.Session, optional
        Provide your own session for connection pooling, custom headers,
        proxies, etc.
    """

    def __init__(
        self,
        base_url: str = DEFAULT_API_URL,
        timeout: int = DEFAULT_TIMEOUT,
        retries: int = DEFAULT_RETRIES,
        session: Optional[requests.Session] = None,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.retries = retries
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", f"nativeready-python/0.4.1")
        self.session.headers.setdefault("Content-Type", "application/json")

    # ------------------------------------------------------------------
    # Single-sequence prediction
    # ------------------------------------------------------------------
    def predict(
        self,
        sequence: str,
        uniprot_id: Optional[str] = None,
    ) -> PredictionResult:
        """Predict for one sequence (FASTA header optional)."""
        # If user passed a FASTA-formatted string with a header, try to
        # extract the UniProt accession before stripping the header.
        if uniprot_id is None:
            uniprot_id = _extract_uniprot_from_header(sequence)
        cleaned = _clean_sequence(sequence)
        if len(cleaned) < 10:
            raise SequenceError(
                f"Sequence too short ({len(cleaned)} aa, minimum 10)"
            )
        return PredictionResult.from_api(
            self._post_predict(cleaned),
            uniprot_id=uniprot_id,
            sequence_length=len(cleaned),
        )

    def predict_uniprot(self, accession: str) -> PredictionResult:
        """Fetch the FASTA from UniProt and predict.

        Parameters
        ----------
        accession : str
            A canonical UniProt accession (e.g., 'P00918' for human carbonic anhydrase 2).
        """
        accession = accession.strip().upper()
        if not re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}", accession):
            raise SequenceError(f"'{accession}' does not look like a UniProt accession")
        try:
            r = requests.get(UNIPROT_FASTA_URL.format(accession), timeout=30)
            r.raise_for_status()
        except requests.RequestException as e:
            raise APIError(f"Failed to fetch UniProt {accession}: {e}") from e
        return self.predict(r.text, uniprot_id=accession)

    # ------------------------------------------------------------------
    # Batch prediction
    # ------------------------------------------------------------------
    def predict_batch(
        self,
        sequences: Sequence[Union[str, Dict[str, str]]],
        progress: bool = True,
    ) -> List[PredictionResult]:
        """Predict for many sequences in a row.

        Parameters
        ----------
        sequences : iterable
            Each element may be:
              - a raw amino acid string,
              - a FASTA-formatted string (with or without header),
              - a dict {'id': ..., 'sequence': ...} for explicit labeling.
        progress : bool
            Show a progress indicator if `tqdm` is installed.

        Returns
        -------
        list of PredictionResult
        """
        items = list(sequences)
        iterator = enumerate(items)
        try:
            from tqdm import tqdm  # type: ignore
            if progress:
                iterator = tqdm(iterator, total=len(items), desc="NativeReady")
        except ImportError:
            pass

        results: List[PredictionResult] = []
        for i, item in iterator:
            try:
                if isinstance(item, dict):
                    seq = item.get("sequence") or item.get("seq") or ""
                    uid = item.get("id") or item.get("uniprot_id")
                    results.append(self.predict(seq, uniprot_id=uid))
                else:
                    results.append(self.predict(item))
            except (SequenceError, APIError) as e:
                # Insert a placeholder so the result list aligns with input
                results.append(PredictionResult(
                    score=-1, label="ERROR",
                    confidence_lower=-1, confidence_upper=-1,
                    model_version="", risk_factors=[],
                    recommendations=[f"Prediction failed: {e}"],
                    raw_response={"error": str(e)},
                ))
        return results

    def predict_fasta(
        self,
        path_or_text: str,
        progress: bool = True,
    ) -> List[PredictionResult]:
        """Read a FASTA file (or raw FASTA text) and predict for every record."""
        from .fasta import parse_fasta
        records = parse_fasta(path_or_text)
        return self.predict_batch(records, progress=progress)

    # ------------------------------------------------------------------
    # Health and metadata
    # ------------------------------------------------------------------
    def health(self) -> Dict[str, str]:
        """Quick health check against the API."""
        r = self._request("GET", "/")
        return r

    def feedback_stats(self) -> Dict[str, Any]:
        """Public aggregate counts of community-submitted experimental outcomes."""
        return self._request("GET", "/feedback/stats")

    def report_outcome(
        self,
        sequence: str,
        predicted_score: int,
        outcome: str,
        note: Optional[str] = None,
        buffer: Optional[str] = None,
        construct: Optional[str] = None,
        expression_system: Optional[str] = None,
        instrument: Optional[str] = None,
        resolution: Optional[str] = None,
        failure_mode: Optional[str] = None,
        model_version: Optional[str] = None,
        email_for_followup: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Report what actually happened when you ran a protein on the instrument.

        Real experimental outcomes are the data the model cannot get any other
        way, since the literature almost never reports native-MS failures. Every
        outcome you send back (success or failure) becomes a training label and
        makes the model better for everyone. This is the one call that closes the
        loop from prediction to reality.

        Parameters
        ----------
        sequence : str
            The sequence you tested. Hashed before storage; the raw protein is
            never stored server-side.
        predicted_score : int
            The 0-100 score the model returned for this sequence (use
            ``result.score`` from a PredictionResult).
        outcome : str
            What happened in the lab: ``"worked"``, ``"failed"``, or
            ``"not_tested"`` (planned but not yet run).
        note : str, optional
            Free-text context, e.g. buffer, instrument, observed failure mode.
        model_version : str, optional
            The model version the prediction came from (use
            ``result.model_version``) so the outcome is tied to the right model.
        email_for_followup : str, optional
            If given, opts you in to a single follow-up message in 2-4 weeks.

        Returns
        -------
        dict
            The API acknowledgement (``{"status": "received", ...}``).

        Examples
        --------
        >>> client = NativeReady()
        >>> r = client.predict("MKT...")
        >>> # ...you run it on the instrument...
        >>> client.report_outcome(
        ...     "MKT...", r.score, "failed",
        ...     buffer="200 mM ammonium acetate, pH 7.5",
        ...     construct="residues 1-256 ectodomain",
        ...     instrument="Thermo Q Exactive UHMR",
        ...     failure_mode="no resolvable signal",
        ...     model_version=r.model_version,
        ... )

        The condition fields (buffer, construct, expression_system, instrument,
        resolution, failure_mode) are all optional. Providing them turns each
        outcome into a conditioned observation, which is what makes a "failed"
        result meaningful: the same protein can fail in one buffer and work in
        another.
        """
        valid = ("worked", "failed", "not_tested")
        if outcome not in valid:
            raise ValueError(f"outcome must be one of {valid}, got {outcome!r}")
        cleaned = _clean_sequence(sequence)
        if len(cleaned) < 10:
            raise SequenceError(
                f"Sequence too short ({len(cleaned)} aa, minimum 10)"
            )
        payload: Dict[str, Any] = {
            "sequence": cleaned,
            "predicted_score": int(predicted_score),
            "user_outcome": outcome,
        }
        if note is not None:
            payload["note"] = note
        # Structured conditions: send only what the caller provided.
        for _k, _v in (
            ("buffer", buffer),
            ("construct", construct),
            ("expression_system", expression_system),
            ("instrument", instrument),
            ("resolution", resolution),
            ("failure_mode", failure_mode),
        ):
            if _v is not None:
                payload[_k] = _v
        if model_version is not None:
            payload["model_version"] = model_version
        if email_for_followup is not None:
            payload["email_for_followup"] = email_for_followup
        return self._request("POST", "/feedback", json=payload)

    # ------------------------------------------------------------------
    # Internal HTTP plumbing
    # ------------------------------------------------------------------
    def _post_predict(self, sequence: str) -> Dict[str, Any]:
        return self._request("POST", "/predict", json={"sequence": sequence})

    def _request(
        self,
        method: str,
        path: str,
        json: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        url = self.base_url + path
        last_exc: Optional[Exception] = None
        for attempt in range(self.retries + 1):
            try:
                r = self.session.request(
                    method, url, json=json, timeout=self.timeout
                )
                if r.status_code == 200:
                    return r.json()
                if r.status_code == 400:
                    raise APIError(f"400 Bad Request: {r.text[:200]}")
                if r.status_code in (502, 503, 504) and attempt < self.retries:
                    time.sleep(2 ** attempt)
                    continue
                raise APIError(f"HTTP {r.status_code}: {r.text[:200]}")
            except requests.Timeout as e:
                last_exc = e
                if attempt < self.retries:
                    time.sleep(2 ** attempt)
                    continue
                raise APIError(f"Timeout after {self.timeout}s: {e}") from e
            except requests.RequestException as e:
                last_exc = e
                if attempt < self.retries:
                    time.sleep(2 ** attempt)
                    continue
                raise APIError(f"Request failed: {e}") from e
        raise APIError(f"All {self.retries + 1} attempts failed: {last_exc}")
