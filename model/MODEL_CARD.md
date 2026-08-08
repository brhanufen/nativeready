# NativeReady Model Card

**Version:** 0.4-esm2-glyco-tm (served model, `predictor_v3`)
**Date:** 2026-08-07 (rewritten from the stale v0.2 card)
**License:** MIT (code and model)
**Status:** Preprint-stage open benchmark and triage model. Not a validated failure-exclusion tool.

---

## 1. Model details

### Description
A sequence-based model that predicts whether a protein is likely to yield usable
native mass spectrometry (native MS) data. Returns a calibrated suitability score
(0-100), a label (Excellent/Good/Fair/Poor/Unsuitable), a confidence interval, an
out-of-distribution (OOD) flag, and (v0.4) a physics-computed heterogeneity annotation
for glycoproteins and antibody-drug conjugates.

### Architecture (v0.4 served model)
- **Served entry point:** `backend/predictor_v3.py`, `model_version = "0.4-esm2-glyco-tm"`.
- **Feature representation:** 47 sequence-derived features (BioPython physicochemical
  + composition + v0.4 glycosylation and transmembrane-topology features), optionally
  concatenated with ESM-2 PCA components in the combined variant.
- **Base estimator:** scikit-learn RandomForestClassifier with isotonic calibration
  (CalibratedClassifierCV).
- **OOD detection:** nearest-neighbor distance in standardized feature space; inputs
  beyond the 95th-percentile training distance are flagged with a warning and a wider
  interval.
- **ESM-2:** `facebook/esm2_t33_650M_UR50D` embeddings used for the combined variant and
  for the cluster-aware validation split. See the honest note in Section 6 on ESM-2's
  limited contribution.

### Inputs
- Protein sequence (single-letter amino acids, FASTA-cleaned). Length 10-5000 aa.
  Ambiguous residues mapped to nearest standard; X dropped.

### Outputs
```json
{
  "suitability_score": 0,
  "suitability_label": "Excellent | Good | Fair | Poor | Unsuitable",
  "confidence_interval": {"lower": 0, "upper": 0},
  "risk_factors": [],
  "recommendations": [],
  "heterogeneity_risk": {"level": "high", "reason": "...", "computed": true},
  "model_version": "0.4-esm2-glyco-tm",
  "warning": null
}
```
The `heterogeneity_risk` field appears only for glycoproteins/ADCs and is a physics
calculation (Method 1), not a learned prediction; see Section 5.

---

## 2. Intended use

### In scope
- First-pass positive triage: ranking and prioritizing candidate proteins likely to be
  native-MS-compatible, to allocate limited instrument time.
- Ranking variants of the same protein by predicted suitability.
- Physics-based heterogeneity screening for glycoproteins and ADCs (Method 1).
- An open benchmark and reproducible baseline for the community.

### Out of scope
- Definitive exclusion / go-no-go on expensive experiments. The model is a positive-
  triage tool, not a validated failure predictor (Section 6).
- Clinical, diagnostic, or regulatory use.
- Replacing empirical optimization or expert judgment.

---

## 3. Training data

### Source
All sequences fetched from the UniProt REST API. Zero synthetic data.

### Composition (v7 dataset, `data/dataset_combined_v7_2026-05-11.json`)
- **n = 635 unique proteins**: 538 positives, 97 negatives.
- **Positives (538):** proteins with documented native MS outcomes in the published
  literature (soluble standards, tetramers, antibodies, membrane proteins, complexes,
  viral capsids, intrinsically disordered proteins).
- **Negatives (97):** 64 randomly sampled Swiss-Prot proteins (proxy negatives), ~30
  property-targeted hard cases (mucins, very large membrane proteins, aggregation-prone),
  and only ~3 evidence-based documented experimental failures.

### Honest note on the negatives
64 of the 97 negatives were never experimentally tested by native MS; they are proxies,
not confirmed failures. Only ~3 are evidence-based real failures. This is the central
data limitation and it governs how the model may be used (Section 6).

---

## 4. Features (n = 47)

BioPython physicochemical + amino-acid composition (36, as in v0.2) plus 11 v0.4
features: N-glycosylation sequon count / strict count / density, mucin-like window
count and fraction, percent S/T/P, and transmembrane topology (TM-helix count, TM
residue total and fraction, longest TM helix, polytopic flag). Full list in
`model/feature_names.json`.

---

## 5. Evaluation (honest, cluster-aware primary)

### Primary: cluster-aware 5-fold cross-validation
GroupKFold over 80 ESM-2-similarity clusters (cosine distance 0.10 threshold), so
close homologs cannot split across train and test. This is the honest split; standard
stratified CV inflates results via homology leakage.

- **Combined model (V4) cluster-aware AUC = 0.835 +/- 0.029** (303-feature combined
  variant; the manuscript reports 0.869 for the primary combined configuration).
- Per-fold AUCs: 0.873, 0.801, 0.854, 0.802, 0.845.
- BioPython-only baseline (V1) cluster-aware AUC = 0.838; ESM-2 linear probe = 0.851.

### Read the AUC with this caveat (size confound)
Negatives are systematically longer than positives (median 777 vs 288 aa). On the same
pooled out-of-fold predictions, sequence length alone reaches AUC 0.769 and molecular
weight alone 0.773, versus the combined model's 0.847. Much of the apparent
discrimination reflects a size confound in the negative sampling, not learned native-MS
physics.

### Failure-detection performance: not yet established
We do not report a single negative-recall number. Because 64 of 97 negatives were never
tested, counting them as ground-truth failures assumes a 100 percent failure rate among
untested proteins. Under any plausible true failure rate at or below 30 percent, the
model already flags as many proxy negatives as a perfect classifier could. The true
native-MS failure rate on the deployment distribution is not identifiable from these
data and must be measured experimentally. Failure-detection performance is therefore
reported as not yet established, not as a percentage.

### Model-comparison note
A paired signed-rank test on five folds cannot reach significance at alpha = 0.05 (its
minimum attainable two-sided p is 0.0625), so it is not used to compare variants. The
observation that a length-only rule reaches AUC 0.769 against the model's 0.847 is the
basis for the conclusion that ESM-2 adds little discriminative signal here.

### Heterogeneity annotation (Method 1)
For glycoproteins and ADCs, the served model adds a physics-computed heterogeneity flag:
it simulates the proteoform mass envelope from N-glycosylation sequons and computes the
instrument resolution required to resolve it. This is arithmetic from published glycan
masses, validated against successes, not a learned classifier. It is a heterogeneity
detector, not a success predictor: heavily heterogeneous proteins that are nonetheless
studied by native MS (for example SARS-CoV-2 Spike) will be flagged. Assumptions
(75 percent sequon occupancy, 1 percent detection threshold, O-glycosylation not
modeled) are emitted with every report.

---

## 6. Known limitations

1. **Failure prediction is not validated.** The negative class is dominated by untested
   proxies; only ~3 real documented failures exist in the data. Use the model for
   positive triage, not exclusion.
2. **AUC is partly a length/size confound** (Section 5). Length alone recovers most of
   the model's skill-above-chance.
3. **ESM-2 adds little on this task.** Under the honest cluster-aware split the ESM-2
   signal is weak; a one-class experiment (`training_report_v5_pu.json`) shows ESM-2
   embeddings at near-chance for separating the class once leakage is removed. The
   biophysical features carry most of the usable signal.
4. **The true base rate is unknown.** The native-MS success rate on the deployment
   distribution is not identifiable from these data and requires a prospective
   measurement.
5. **Sequence only, no structure.** Proteins with similar sequence but different fold or
   oligomeric state receive similar predictions. OOD inputs are flagged.
6. **Publication bias in the labels.** The training label is closer to
   published-native-MS-success than to native-MS-suitability, which biases the model
   toward well-studied proteins. This is a property of the field's data, not the pipeline.

---

## 7. Ethical considerations
- No patient or human-subject data; all inputs are public protein sequences.
- Code, training pipeline, dataset, and weights are open (MIT) and reproducible.
- Limitations are disclosed on the public site and in the preprint.

---

## 8. Reproducibility
- Code: `/Users/bfentaw2/startup/nativeready/`
- Served predictor: `backend/predictor_v3.py`
- Training: `model/train_v3_cluster_aware.py`; report `model/training_report_v3_robust.json`
- Out-of-fold predictions: `model/v3_robust_oof_predictions.npz`
- Dataset: `data/dataset_combined_v7_2026-05-11.json`
- Dependencies: `backend/requirements.txt`

---

## 9. Citation
> Znabu, B. F. and Atif, Z. NativeReady: an open benchmark and sequence-based triage
> model for native mass spectrometry suitability. bioRxiv 2026.05.03.722506. https://nativeready.bio

---

*Follows the Mitchell et al. (2018) Model Cards framework. This card supersedes the
stale v0.2 card and describes the served v0.4 model honestly, including the corrections
in INTEGRITY_DISCLOSURE_DRAFT.md.*
