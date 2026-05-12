"""NativeReady Python SDK.

Predict native mass spectrometry suitability for a protein sequence.

Quick start
-----------
>>> from nativeready import predict
>>> result = predict("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
>>> print(result.score, result.label)
97 Excellent

For batch predictions, FASTA files, and UniProt accession lookup, use the Client:

>>> from nativeready import Client
>>> client = Client()
>>> results = client.predict_fasta("my_proteins.fasta")
>>> for r in results:
...     print(r.uniprot_id, r.score, r.label)

See https://github.com/brhanufen/nativeready for the full documentation,
the open dataset, and the trained model.
"""
from .client import Client, PredictionResult
from .fasta import parse_fasta

__version__ = "0.4.0"

# Module-level convenience function for the simplest possible API.
_default_client = None


def predict(sequence, **kwargs):
    """Predict native MS suitability for a single sequence.

    Convenience wrapper around Client().predict(). Caches a module-level
    Client so repeated calls reuse the HTTP connection.

    Parameters
    ----------
    sequence : str
        Protein sequence (FASTA header optional, will be stripped).
    **kwargs
        Forwarded to Client().predict().

    Returns
    -------
    PredictionResult

    Examples
    --------
    >>> from nativeready import predict
    >>> r = predict("MQIFVKTLTGKTITLEVEPSDTI...")
    >>> r.score, r.label, r.is_ood
    (97, 'Excellent', False)
    """
    global _default_client
    if _default_client is None:
        _default_client = Client()
    return _default_client.predict(sequence, **kwargs)


__all__ = [
    "Client",
    "PredictionResult",
    "parse_fasta",
    "predict",
    "__version__",
]
