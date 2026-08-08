# nativeready

Python SDK for the [NativeReady](https://nativeready.bio) API. Predict whether a protein sequence is likely to give usable native mass spectrometry data, in seconds.

```bash
pip install nativeready
```

## Quick start

```python
from nativeready import predict

result = predict("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
print(result)
# PredictionResult(score=97, label='Excellent', ci=[88-100], model='0.4-esm2-glyco-tm')

print(result.score)         # 97
print(result.label)         # 'Excellent'
print(result.is_ood)        # False
print(result.recommendations)
```

## Common use cases

### Predict from a FASTA file

```python
from nativeready import Client

client = Client()
results = client.predict_fasta("my_proteins.fasta")
for r in results:
    print(r.uniprot_id, r.score, r.label)
```

### Predict from a UniProt accession

```python
from nativeready import Client

client = Client()
result = client.predict_uniprot("P00918")  # Carbonic anhydrase 2
print(result.score, result.label)
```

### Batch prediction with progress bar

```python
from nativeready import Client

client = Client()
sequences = [
    {"id": "ubiquitin", "sequence": "MQIFVKTLTGKTITLEV..."},
    {"id": "lysozyme",  "sequence": "KVFGRCELAAAMKR..."},
    "MSHHWGYGKHNGPEHWHKDF...",  # raw string also works
]
results = client.predict_batch(sequences)  # tqdm progress bar if installed
```

### Report a real outcome (close the loop)

Predictions are only half the story. When you actually run a protein on the
instrument, report back what happened. Real outcomes, especially failures, are
the data the model cannot get from the literature, and every one you send makes
the model better for the whole community.

```python
from nativeready import Client

client = Client()
r = client.predict("MQIFVKTLTGKTITLEV...")

# ...you run it on the instrument, and it fails to give resolvable signal...
client.report_outcome(
    "MQIFVKTLTGKTITLEV...",
    r.score,
    "failed",                       # "worked" | "failed" | "not_tested"
    note="200 mM AmAc pH 7.5, Q Exactive UHMR, no resolvable charge series",
    model_version=r.model_version,  # ties the outcome to the model that predicted it
)
```

The sequence is hashed before storage; the raw protein is never stored
server-side. Aggregate community outcomes are visible via
`client.feedback_stats()`.

### Pandas DataFrame output (with `pandas` extra)

```bash
pip install nativeready[pandas]
```

```python
import pandas as pd
from nativeready import Client

client = Client()
results = client.predict_fasta("my_proteins.fasta")
df = pd.DataFrame([r.as_dict() for r in results])
df.to_csv("results.csv", index=False)
```

## Command line interface

The package installs a `nativeready` CLI:

```bash
# Single sequence
nativeready predict --sequence "MQIFVKTLTGKTITLEV..."

# UniProt accession
nativeready predict --uniprot P00918

# Multiple UniProt accessions (comma-separated)
nativeready predict --uniprot P00918,P0CG48,P00698

# FASTA file with CSV output
nativeready predict --fasta proteins.fasta --output results.csv

# FASTA file with JSON output
nativeready predict --fasta proteins.fasta --output results.json

# Health check
nativeready health

# Public feedback statistics
nativeready stats
```

## What the predictions mean

Each `PredictionResult` contains:

- `score` — calibrated suitability score, 0 to 100
- `label` — `Excellent` (>= 80), `Good` (>= 65), `Fair` (>= 50), `Poor` (>= 35), `Unsuitable` (< 35)
- `confidence_lower`, `confidence_upper` — 95 percent confidence interval (wider when out-of-distribution)
- `is_ood` — `True` if the sequence is unusual relative to training data; trust the score with extra caution
- `risk_factors` — per-feature risk levels (length, MW, hydrophobicity, pI, instability, cysteine content)
- `recommendations` — buffer choice, sample-prep guidance, and instrument-mode notes
- `model_version` — server-side model identifier (e.g., `0.4-esm2-glyco-tm`)

## API server

By default the SDK calls `https://nativeready-production.up.railway.app`. To use a self-hosted or local deployment:

```python
from nativeready import Client
client = Client(base_url="http://localhost:8000")
```

## Honest scope

NativeReady is currently most reliable as a **positive-suitability triage tool**, not a validated failure detector. With only two evidence-based real-failure records in the training set, the negative-class performance is not yet statistically meaningful. High-confidence positive predictions can be trusted; low-confidence predictions should be treated as a flag for manual review, not a verdict. See the [bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.05.03.722506v1) for full methodology and limitations.

## Citing

```
Znabu BFZ, Atif Z. NativeReady: an open benchmark and sequence-based triage
model for native mass spectrometry suitability. bioRxiv, 2026.
https://doi.org/10.64898/2026.05.03.722506
```

## License

MIT. See `LICENSE`.

## Links

- API server source: https://github.com/brhanufen/nativeready
- Web tool: https://nativeready.bio
- Open dataset (634 proteins, CC-BY 4.0): in the `data/` folder of the main repository
- Issues: https://github.com/brhanufen/nativeready/issues
