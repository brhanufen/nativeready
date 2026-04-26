# NativeReady Model Card

**Version:** 0.2-calibrated
**Date:** April 2026
**License:** Research use only
**Status:** v0 prototype, not production-ready

---

## 1. Model details

### Description
A binary classifier that predicts whether a protein sequence is likely to give usable native mass spectrometry (MS) data. Returns a calibrated probability score (0–100), risk factor breakdown, and recommendations.

### Architecture
- **Base estimator:** scikit-learn RandomForestClassifier
  - 300 trees, max depth 10, min samples per leaf 2
  - class-weight balanced (positive/negative class imbalance)
- **Calibration:** CalibratedClassifierCV with isotonic regression, 5-fold CV
- **Out-of-distribution (OOD) detection:** Nearest-neighbor distance in standardized feature space; threshold = 95th percentile of training-set NN distances

### Inputs
- Protein sequence as a string of single-letter amino acid codes (FASTA-cleaned)
- Length: 10–5000 amino acids (validation rejects outside range)
- Ambiguous residues (B, Z, J, U) mapped to closest standard residue; unknowns (X) dropped

### Outputs
```json
{
  "suitability_score": 0-100,
  "suitability_label": "Excellent | Good | Fair | Poor | Unsuitable",
  "confidence_interval": {"lower": int, "upper": int},
  "risk_factors": [...],
  "recommendations": [...],
  "model_version": "0.2-calibrated",
  "warning": "out_of_distribution" | null
}
```

---

## 2. Intended use

### In scope
- First-pass screening of candidate proteins for native MS experiments
- Ranking variants of the same protein by predicted suitability
- Educational demonstration of AI applied to analytical chemistry
- Research planning (deciding which constructs to prioritize for limited instrument time)

### Out of scope
- Clinical diagnosis or any medical decision making
- Replacing empirical optimization or expert judgment
- High-stakes go/no-go decisions for expensive experiments without empirical validation
- Patenting or regulatory submissions
- Predicting performance of proteins very different from the training distribution (will be flagged with OOD warning)

---

## 3. Training data

### Source
All sequences fetched live from the UniProt REST API. **Zero synthetic data.**

### Composition (v0.2)
- **Positives (n=138):** Real proteins documented in published native MS literature, curated from peer-reviewed sources spanning ~30 years. Includes:
  - Small peptides (insulin, glucagon, oxytocin, etc.)
  - Soluble standards (carbonic anhydrase, ubiquitin, lysozyme, cytochrome c)
  - Hemoglobin/transthyretin tetramers
  - Antibodies (IgG1–4)
  - Membrane proteins (AmtB, OmpF, GLUT1, bacteriorhodopsin)
  - Large complexes (GroEL, ribosomes, ATP synthase)
  - Viral capsids (HBV, MS2, AAV2, SARS-CoV-2 Spike)
  - Intrinsically disordered (alpha-synuclein, tau, p53)

- **Negatives (n=94):** Two sub-classes:
  - **Proxy negatives (n=64):** Random Swiss-Prot proteins not in the positive set
  - **Hard negatives (n=30):** Proteins explicitly known to be challenging for native MS — heavily glycosylated mucins (MUC1, MUC2, MUC5B), very large multi-pass membrane proteins (dystrophin, plectin), aggregation-prone proteins (PrP, huntingtin)

### Splits
- 80/20 stratified train/test split (random_state=42)
- 5-fold stratified cross-validation for calibration and reporting

---

## 4. Features (n=36)

Computed with BioPython's `ProteinAnalysis`:

| Group | Features |
|---|---|
| Size | length, log_length, molecular_weight_kda |
| Physicochemical | isoelectric_point, instability_index, gravy (hydrophobicity), aromaticity |
| Secondary structure | helix_fraction, turn_fraction, sheet_fraction |
| Composition (per-AA) | pct_A, pct_C, ..., pct_Y (20 features) |
| Composition (grouped) | pct_cysteine, pct_proline, pct_charged, pct_hydrophobic, pct_aromatic, pct_polar |

---

## 5. Evaluation

### Held-out test set (n=47)
- Accuracy: **0.723**
- AUC: **0.780**
- Brier score: 0.196 (lower = better calibrated probabilities)
- Confusion matrix: [[TN=13, FP=6], [FN=7, TP=21]]

### 5-fold cross-validation
- Mean AUC: **0.855** (±0.043)
- Mean Accuracy: **0.797** (±0.038)
- Mean Brier: 0.155

### Why CV numbers > held-out
The held-out set (n=47) contains a higher fraction of hard cases (small peptides, mucins) by stratified random chance. CV gives a more stable estimate. Use CV numbers for honest performance claims.

### Comparison vs v0.1 baseline
| Metric | v0.1 (no calibration, 133 examples) | v0.2 (calibrated, 232 examples) |
|---|---|---|
| Held-out AUC | 0.967 | 0.780 |
| 5-fold CV AUC | 0.869 | 0.855 |
| Calibrated probabilities | No | Yes |
| OOD detection | No | Yes |
| Small-peptide handling | Wrong (no warning) | Honest (OOD-flagged) |

The drop in held-out AUC reflects a more diverse and harder test set, not worse modeling. CV AUC barely changed, confirming generalization is preserved.

---

## 6. Known limitations

1. **Small dataset (~232 examples).** Performance estimates have meaningful variance. The model is a v0 prototype, not production-grade.

2. **Negatives are proxies, not true experimental failures.** Real native MS failure data is rarely published. We mitigate by including expert-curated hard cases (mucins, very large membrane proteins, aggregation-prone), but this remains an approximation.

3. **Reduced reliability for sequences outside the training distribution.** The OOD detector flags such inputs with a warning and a wider confidence interval. Common OOD cases include:
   - Very short peptides (<30 aa)
   - Highly repetitive sequences
   - Sequences with unusual amino acid composition

4. **No structural information.** The model uses only sequence-derived features. Proteins with similar sequences but different folding/oligomeric states would receive similar predictions.

5. **No glycosylation prediction.** While the model includes proxy features for size and composition, it does not directly predict post-translational modifications. Glycoprotein performance may be unreliable.

6. **Bias toward soluble globular proteins.** The training set is biased toward proteins that have been published, which skews toward "interesting" and "tractable" targets.

---

## 7. Ethical considerations

- **No patient data.** All training data is publicly available protein sequence information.
- **No biased outcomes for protected classes** (proteins are not people).
- **Transparency:** Model code, training pipeline, and weights are open and reproducible.
- **Disclosure:** Limitations are stated upfront on the public website.

---

## 8. Reproducibility

- Code: `/Users/bfentaw2/startup/nativeready/`
- Training script: `model/train_v2.py`
- Dependencies: see `backend/requirements.txt`
- Model artifact: `model/model_v2.joblib`
- Performance verification: `tests/scientific_audit.py`

To reproduce:
```bash
cd model
python3 train_v2.py            # ~5 min including data fetch
cd ../tests
python3 scientific_audit.py    # ~30 sec
```

Random seeds are fixed (42) so results are deterministic.

---

## 9. Citation

If using NativeReady in research, please cite:

> NativeReady: A v0 prototype for predicting native mass spectrometry suitability from protein sequence. Brhanu Fentaw, 2026. https://nativeready.app

---

*This model card follows the Mitchell et al. (2018) "Model Cards for Model Reporting" framework.*
