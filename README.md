# NativeReady

**Predict whether a protein will give usable native mass spec data, before you run the experiment.**

A sequence-based triage model and open benchmark for native MS suitability. Sequence in, calibrated probability out, in seconds.

Live tool: https://nativeready.netlify.app
Python SDK: `pip install nativeready` ([PyPI](https://pypi.org/project/nativeready/))
Preprint: [bioRxiv 2026.05.03.722506](https://doi.org/10.64898/2026.05.03.722506) (CC-BY 4.0)

## Install (Python SDK)

```bash
pip install nativeready
```

```python
from nativeready import predict

result = predict("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
print(result.score, result.label)  # 97 Excellent
```

For batch predictions, FASTA file support, and CLI usage, see [`python-sdk/README.md`](python-sdk/README.md).

---

## What it does

Paste a protein sequence (FASTA or raw amino acids). The model returns:
- A **suitability score** (0-100) for native MS, calibrated by isotonic regression
- A **confidence interval**
- **Risk factors** (length, MW, hydrophobicity, pI, instability, Cys content)
- **Recommendations** (buffer choice, sample prep, when to denature first)
- An **out-of-distribution flag** when the input is unusual relative to training data

## What's new in v0.3

- **634 unique proteins** in the training set (up from 232 in v0.2), drawn from a hand-curated literature base, RCSB PDB full-text searches for native MS terms, and EuropePMC supplementary mining
- **ESM-2 protein-language embeddings** (facebook/esm2_t33_650M_UR50D) added as a feature representation alongside the original 36 BioPython physicochemical features
- **Cluster-aware cross-validation** (GroupKFold over ESM-2 embedding-similarity clusters) to defend against homology leakage between train and test folds
- **OOD detector retrained** on the v0.3 feature distribution
- **Schema released** (`data/LABEL_SCHEMA.md`) with a 5-level ordinal outcome label and a 7-term failure_mode controlled vocabulary

## Performance (cluster-aware 5-fold CV, n = 634)

| Variant | Features | ROC-AUC |
|---|---:|---|
| V1 BioPython only (baseline) | 36 | 0.852 +/- 0.074 |
| V2 ESM-2 linear probe | 1,280 | 0.842 +/- 0.040 |
| V3 ESM-2 + PCA + RF | up to 256 | 0.821 +/- 0.057 |
| **V4 combined (production)** | **up to 292** | **0.869 +/- 0.036** |

Positive recall under V4 cluster-aware: **99.4%**.
Negative recall under V4 cluster-aware: **9.4%** (limited by scarcity of evidence-based real-failure data; see preprint Section 4.2).

## Honest scope

NativeReady is currently most reliable as a **positive-suitability triage tool**, not as a validated failure detector. With only 2 of 96 negatives being evidence-based real experimental failures (the remaining 94 are 64 random Swiss-Prot proxies and 30 curated property-targeted records), the negative-class evaluation is not yet statistically meaningful. The full preprint discusses this limitation in detail and proposes a user-contribution mechanism to accumulate real failure outcomes over time.

Use a high-confidence positive prediction with reasonable trust. Treat a low-confidence prediction as a flag for manual review, not a verdict.

## Project layout

```
nativeready/
├── README.md                                 ← this file
├── NativeReady_preprint_v1.docx              ← bioRxiv-ready manuscript
├── LICENSE                                    ← MIT
├── Dockerfile                                 ← Railway deployment with pre-baked ESM-2
├── railway.toml
│
├── data/
│   ├── LABEL_SCHEMA.md                       ← v0.1 schema, ordinal label + failure_mode CV
│   ├── dataset_combined_v4_2026-05-02.json   ← canonical 634-protein release
│   ├── positives_with_sequences.json         ← original 69-protein curated base
│   ├── negatives_with_sequences.json         ← 64 Swiss-Prot proxies
│   ├── expansion_with_sequences.json         ← 99 expansion records
│   ├── extracted_new_2026-05-01.json         ← 260 PDB-extracted positives (May 1)
│   ├── extracted_pdb_expanded_2026-05-02.json ← 98 PDB-extracted positives (May 2)
│   ├── extracted_pmc_pilot_2026-05-02.json   ← 8 EuropePMC records (May 2)
│   ├── extracted_pmc_batch2b_2026-05-02.json ← 36 EuropePMC records (May 2)
│   └── extract_*.py, fetch_*.py              ← reproducible extraction scripts
│
├── model/
│   ├── model_v3.joblib                       ← v0.3 production model (combined ESM-2 + BioPython)
│   ├── ood_detector_v3.joblib                ← v0.3 OOD detector
│   ├── esm2_embeddings_634.npy               ← cached ESM-2 embeddings for all 634 proteins
│   ├── esm2_embeddings_metadata.json
│   ├── training_report_v3.json               ← stratified CV results
│   ├── training_report_v3_robust.json        ← cluster-aware CV results
│   ├── compute_esm2_embeddings.py            ← embedding pipeline
│   ├── train_v3_compare.py                   ← stratified CV training
│   ├── train_v3_cluster_aware.py             ← cluster-aware CV training
│   ├── model_v2.joblib                       ← v0.2 fallback (BioPython-only)
│   └── train.py, train_v2.py                 ← legacy training scripts
│
├── backend/
│   ├── main.py                               ← FastAPI service
│   ├── predictor_v3.py                       ← v0.3 inference (ESM-2 + BioPython)
│   ├── predictor_v2.py                       ← v0.2 fallback inference
│   ├── features.py                           ← BioPython feature extraction
│   ├── requirements.txt                      ← includes torch + transformers
│   └── README.md
│
├── frontend/
│   ├── index.html, style.css, script.js
│   └── README.md
│
├── figures/
│   ├── Figure1_pipeline.png
│   ├── Figure2_ROC_calibration.png
│   ├── Figure3_per_class_accuracy.png
│   └── make_figures.py
│
├── tests/
└── docs/
```

## API endpoints

| Method | Path | Purpose |
|---|---|---|
| GET  | `/` | Health check |
| POST | `/predict` | Predict native MS suitability for a sequence |
| POST | `/feedback` | Submit experimental outcome (worked / failed / not_tested) |
| GET  | `/feedback/stats` | Public aggregate counts |
| GET  | `/docs` | Auto-generated Swagger UI |

Example request:
```bash
curl -X POST https://nativeready.netlify.app/api/predict \
  -H "Content-Type: application/json" \
  -d '{"sequence": "MSHHWGYG..."}'
```

## Running locally

```bash
cd backend
python3 -m pip install -r requirements.txt
python3 -m uvicorn main:app --reload --port 8000
```

In another terminal:
```bash
cd frontend
python3 -m http.server 8080
```

Open `http://localhost:8080` in a browser.

## Reproducing the model

```bash
# Compute ESM-2 embeddings (first run downloads ~2.5 GB; cached after)
cd model
python3 compute_esm2_embeddings.py

# Train and compare all 4 variants under stratified CV
python3 train_v3_compare.py

# Train and compare under cluster-aware CV (the headline numbers)
python3 train_v3_cluster_aware.py
```

Reports are saved to `training_report_v3.json` and `training_report_v3_robust.json`.

## Citing NativeReady

```
Znabu BFZ, Atif Z. NativeReady: an open benchmark and sequence-based triage
model for native mass spectrometry suitability. bioRxiv, 2026.
https://doi.org/10.64898/2026.05.03.722506
```

The dataset (`data/dataset_combined_v4_2026-05-02.json`) is released under CC-BY 4.0. The code is released under MIT (see `LICENSE`).

## Honest caveats (read this)

- **Small dataset (n = 634).** Performance estimates have meaningful variance.
- **Negatives are mostly proxy.** 94 of 96 negatives are random or property-targeted Swiss-Prot proteins; only 2 are evidence-based real failures. The model's failure-detection capability is not yet statistically validated.
- **OOD detection** flags sequences very different from training. Trust scores less when this warning appears.
- **Sequences over 1,022 aa are truncated** for ESM-2 inference. 68 proteins (10.7 percent) in the training set were truncated; the model's representation for these reflects only the N-terminal 1,022 residues.
- **No FDA path**: this is a research tool, not a clinical diagnostic.

## License and disclaimer

MIT license (see `LICENSE`). Research prototype, not a medical device, not a diagnostic. Use at your own risk for research planning. All scores are calibrated probability predictions, not guarantees.

---

*v0.3 released May 2026. Solo build with methodology guidance from Zohaib Atif. Public data only, scientific integrity preserved.*
