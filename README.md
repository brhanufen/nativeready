# NativeReady

**Predict whether a protein will give usable native mass spec data — before you run the experiment.**

A solo-built research tool. Sequence-in, prediction-out, in seconds.
Trained on real proteins curated from published native MS literature.

---

## What it does

Paste a protein sequence (FASTA or raw amino acids). The model returns:

- A **suitability score** (0–100) for native MS
- A **confidence interval**
- **Risk factors** (length, MW, hydrophobicity, pI, instability, Cys content)
- **Recommendations** (buffer choice, sample prep, when to denature first)

## Project layout

```
nativeready/
├── README.md                 ← this file
├── build_plan.docx           ← end-to-end build plan with mitigations
├── make_build_plan.py        ← script that generates build_plan.docx
│
├── data/                     ← real public-database datasets
│   ├── positives_curated.py  ← 69 curated native MS validated proteins (UniProt IDs)
│   ├── fetch_sequences.py    ← script: pulls real sequences from UniProt REST API
│   ├── fetch_negatives.py    ← script: samples real Swiss-Prot proxy negatives
│   ├── positives.fasta       ← 69 real sequences (FASTA)
│   ├── positives_with_sequences.json
│   └── negatives_with_sequences.json
│
├── model/                    ← trained classifier
│   ├── train.py              ← training pipeline
│   ├── model.joblib          ← trained RandomForest
│   ├── feature_names.json    ← ordered feature list (must match backend)
│   └── training_report.json  ← real performance metrics
│
├── backend/                  ← FastAPI service
│   ├── main.py               ← API endpoints
│   ├── predictor.py          ← loads model, returns predictions
│   ├── features.py           ← biophysical feature extraction (matches train.py)
│   ├── requirements.txt      ← Python dependencies
│   ├── Dockerfile            ← for cloud deployment
│   └── README.md             ← backend-specific instructions
│
├── frontend/                 ← single-page website (HTML + CSS + JS)
│   ├── index.html            ← the public-facing page
│   ├── style.css             ← dark cyan aesthetic matching Traversa
│   ├── script.js             ← form submission, API call, results rendering
│   └── README.md             ← deployment instructions
│
├── tests/
│   └── test_integration.py   ← end-to-end tests with real sequences
│
└── docs/
    └── DEPLOYMENT.md         ← step-by-step launch guide
```

## How it actually works

1. **Data**: 69 real proteins from native MS literature (UniProt IDs verified) +
   64 proxy negatives sampled from Swiss-Prot. All sequences fetched from
   the UniProt REST API. **No synthetic data.**

2. **Features**: 36 biophysical features per protein computed with BioPython
   (length, MW, pI, GRAVY, instability, secondary structure fractions,
   per-AA composition, grouped composition).

3. **Model**: scikit-learn RandomForest, 200 trees, max depth 8,
   class-balanced. 80/20 stratified train/test split + 5-fold CV.

4. **API**: FastAPI service. POST `/predict` with `{"sequence": "..."}` →
   returns full prediction JSON.

5. **Frontend**: vanilla HTML/CSS/JS. No build step. Drop-in deployable
   to Netlify, Vercel, GitHub Pages, or S3.

6. **Feedback loop**: every prediction includes a one-click feedback card
   ("did your experiment work?"). Anonymous, sequence is hashed before
   storage, stored in `data/feedback.jsonl`. Used to improve the model
   over time as real-world labels accumulate.

## Real performance metrics (v0.2 calibrated)

Held-out test set (n=47):
- **Accuracy: 72.3%**
- **AUC: 0.780**
- **Brier score: 0.196** (lower = better calibrated)

5-fold cross-validation (more reliable for small dataset):
- **Mean AUC: 0.855** (± 0.043)
- **Mean Accuracy: 79.7%** (± 3.8%)
- **Mean Brier: 0.155**

### Why v2 looks worse than v1 on the test set

v0.1 reported 88.9% test accuracy. v0.2 reports 72.3%. **This is honest, not regression.**

v0.2 has:
- **75% more training data** (133 → 232 examples)
- **More diverse test set** including small peptides, mucins, large membrane proteins
- **Calibrated probabilities** (Brier 0.155, well-calibrated)
- **Out-of-distribution detection** (flags sequences unlike training)
- **Per-risk-factor recommendations**

The 5-fold CV AUC barely moved (0.869 → 0.855), confirming generalization is preserved. The held-out drop reflects a harder test set, not a worse model.

Top predictive features (still scientifically sensible): instability_index,
length, molecular_weight, amino acid composition (especially leucine,
phenylalanine), secondary structure sheet content.

## Honest caveats (read this)

- **Small dataset (~232 examples).** Performance estimates have meaningful
  variance. Expect modest generalization.
- **Negatives are proxy examples** (random Swiss-Prot proteins) plus
  expert-curated hard cases (mucins, large membrane proteins,
  aggregation-prone). Still not true experimental failures, since failures
  are rarely published.
- **v0 baseline**, not production-grade. Validate critical experiments
  yourself.
- **OOD detection** flags sequences very different from training. Trust
  scores less when this warning appears.
- **No FDA path**: this is a research tool, not a clinical diagnostic.

See `model/MODEL_CARD.md` for the full ML model card following the
Mitchell et al. (2018) framework.

## Running locally

```bash
# Install dependencies
cd backend
python3 -m pip install --user -r requirements.txt

# Boot the backend
python3 -m uvicorn main:app --reload --port 8000
```

In another terminal:
```bash
cd frontend
python3 -m http.server 8080
```

Open `http://localhost:8080` in your browser.

## Endpoints

| Method | Path | Purpose |
|---|---|---|
| GET | `/` | Health check |
| POST | `/predict` | Predict native MS suitability for a sequence |
| POST | `/feedback` | Submit experimental outcome to improve future model |
| GET | `/feedback/stats` | Public aggregate counts (no individual records) |
| GET | `/docs` | Auto-generated API documentation (Swagger) |

## Running tests

```bash
cd tests
python3 test_integration.py
```

Tests run against the trained model with real sequences (ubiquitin,
carbonic anhydrase 2, real Swiss-Prot negatives).

## Deploying

See `docs/DEPLOYMENT.md` for the step-by-step launch guide:
1. Frontend → Netlify Drop or Vercel (free)
2. Backend → Railway or Render ($5–20/month)
3. Custom domain → Cloudflare Registrar
4. Wire frontend `API_BASE` to the deployed backend URL

## License and disclaimer

Research prototype. Not a medical device. Not a diagnostic. Use at your
own risk for research planning. All risk scores are predictions, not
guarantees.

---

*Built April 2026. Solo build, public data only, scientific integrity preserved.*
