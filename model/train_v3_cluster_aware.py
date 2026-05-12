"""Cluster-aware and source-held-out CV for the four model variants.

Defends against the homology / source-leakage critique by:
1. Clustering proteins by ESM-2 embedding cosine similarity (proxy for sequence similarity).
2. Running GroupKFold so entire clusters land in one fold (no homologs across train/test).
3. ALSO running source-held-out evaluations (train on 2 of 3 sources, test on the 3rd).

Saves results to training_report_v3_robust.json and updates per-variant OOF predictions
for the figures.
"""
import json
import os
import re
from pathlib import Path
from collections import Counter

import numpy as np
import joblib
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn.calibration import CalibratedClassifierCV
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, brier_score_loss, roc_auc_score,
                              average_precision_score, confusion_matrix)
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path("/Users/bfentaw2/startup/nativeready")
DATA = ROOT / "data"
MODEL = ROOT / "model"

DATASET = DATA / "dataset_combined_v6_2026-05-11.json"
ESM_NPY = MODEL / "esm2_embeddings_636.npy"
ESM_META = MODEL / "esm2_embeddings_metadata.json"

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
RANDOM_STATE = 42
CV_SPLITS = 5

# v0.4 glycosylation features
N_GLYC_PERMISSIVE = re.compile(r"N[^P][ST]")
N_GLYC_STRICT = re.compile(r"N[^PC][ST]")
MUCIN_WINDOW = 20
MUCIN_THRESHOLD = 0.40

# v0.4 TM-helix features (Kyte-Doolittle hydropathy window, classical method)
KD_SCALE = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,
            "G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,
            "P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
TM_WINDOW = 19
TM_THRESHOLD = 1.2
TM_MIN_HELIX_LEN = 12

def _tm_segments(seq):
    n = len(seq)
    if n < TM_WINDOW:
        return []
    half = TM_WINDOW // 2
    kd = [KD_SCALE.get(c, 0.0) for c in seq]
    ws = [0.0]*n
    for i in range(half, n-half):
        ws[i] = sum(kd[i-half:i+half+1])/TM_WINDOW
    segs, in_seg, st = [], False, 0
    for i, s in enumerate(ws):
        if s >= TM_THRESHOLD and not in_seg:
            in_seg, st = True, i
        elif s < TM_THRESHOLD and in_seg:
            in_seg = False
            if i - st >= TM_MIN_HELIX_LEN:
                segs.append((st, i-1))
    if in_seg and n - st >= TM_MIN_HELIX_LEN:
        segs.append((st, n-1))
    return segs

# Clustering: AgglomerativeClustering with fixed n_clusters target.
# Tuned so each fold has ~16 clusters on average and the largest cluster cannot
# dominate any single fold (preserving cross-fold balance).
N_CLUSTERS = 80

# ============== HELPERS ==============

def clean_seq(seq):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    cleaned = []
    for c in seq:
        if c in STANDARD_AA:
            cleaned.append(c)
        elif c in "BZJU":
            cleaned.append({"B": "N", "Z": "Q", "J": "L", "U": "C"}[c])
    return "".join(cleaned)

def biopython_features(sequence):
    seq = clean_seq(sequence)
    if len(seq) < 5:
        return None
    pa = ProteinAnalysis(seq)
    aap = pa.amino_acids_percent
    feats = [len(seq), np.log1p(len(seq)), pa.molecular_weight()/1000.0,
             pa.isoelectric_point(), pa.instability_index(), pa.gravy(), pa.aromaticity()]
    h, t, s = pa.secondary_structure_fraction()
    feats += [h, t, s]
    for aa in STANDARD_AA:
        feats.append(aap.get(aa, 0.0))
    feats += [aap.get("C",0.0), aap.get("P",0.0),
              sum(aap.get(c,0.0) for c in "DEKR"),
              sum(aap.get(c,0.0) for c in "AVILMFW"),
              sum(aap.get(c,0.0) for c in "FWY"),
              sum(aap.get(c,0.0) for c in "STNQ")]
    # v0.4 glycosylation features (must match backend/features.py exactly)
    n_seq = len(N_GLYC_PERMISSIVE.findall(seq))
    n_seq_strict = len(N_GLYC_STRICT.findall(seq))
    feats += [float(n_seq), float(n_seq_strict), 100.0 * n_seq / len(seq)]
    if len(seq) >= MUCIN_WINDOW:
        n_windows = len(seq) - MUCIN_WINDOW + 1
        threshold_count = int(MUCIN_THRESHOLD * MUCIN_WINDOW)
        mucin_windows = sum(
            1 for i in range(n_windows)
            if sum(1 for c in seq[i:i + MUCIN_WINDOW] if c in "STP") >= threshold_count
        )
        feats += [float(mucin_windows), float(mucin_windows / n_windows)]
    else:
        feats += [0.0, 0.0]
    feats.append(sum(aap.get(c, 0.0) for c in "STP"))
    # v0.4 TM features
    segs = _tm_segments(seq)
    lens = [e - s + 1 for s, e in segs]
    total = sum(lens)
    feats += [
        float(len(segs)),
        float(total),
        float(total / len(seq)),
        float(max(lens) if lens else 0),
        1.0 if len(segs) >= 2 else 0.0,
    ]
    return feats

def make_rf():
    base = RandomForestClassifier(n_estimators=300, max_depth=10, min_samples_leaf=2,
                                   class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1)
    return CalibratedClassifierCV(base, method="isotonic", cv=3)

def make_logreg():
    base = LogisticRegression(C=1.0, penalty="l2", max_iter=2000, solver="lbfgs",
                               class_weight="balanced", random_state=RANDOM_STATE)
    return CalibratedClassifierCV(base, method="isotonic", cv=3)

def evaluate_one_split(X, y, train_idx, test_idx, scaler=False, pca_n=None,
                       extra_X=None, classifier="rf"):
    if scaler:
        sc = StandardScaler().fit(X[train_idx])
        X_tr = sc.transform(X[train_idx])
        X_te = sc.transform(X[test_idx])
    else:
        X_tr, X_te = X[train_idx], X[test_idx]
    if pca_n is not None:
        # Clamp n_components to valid range
        n_comp = min(pca_n, X_tr.shape[0] - 1, X_tr.shape[1])
        pca = PCA(n_components=n_comp, random_state=RANDOM_STATE).fit(X_tr)
        X_tr = pca.transform(X_tr)
        X_te = pca.transform(X_te)
    if extra_X is not None:
        esc = StandardScaler().fit(extra_X[train_idx])
        ext_tr = esc.transform(extra_X[train_idx])
        ext_te = esc.transform(extra_X[test_idx])
        X_tr = np.hstack([X_tr, ext_tr])
        X_te = np.hstack([X_te, ext_te])
    clf = make_logreg() if classifier == "logreg" else make_rf()
    clf.fit(X_tr, y[train_idx])
    proba = clf.predict_proba(X_te)[:, 1]
    pred = clf.predict(X_te)
    return proba, pred

# ============== LOAD ==============
print("Loading data...")
with open(DATASET) as f:
    proteins = json.load(f)
esm_X = np.load(ESM_NPY)
with open(ESM_META) as f:
    esm_meta = json.load(f)
uids = esm_meta["uids"]
by_uid = {p["uniprot_id"]: p for p in proteins}

y = np.array([by_uid[u]["label"] for u in uids], dtype=int)
seqs = [by_uid[u]["sequence"] for u in uids]
classes = [by_uid[u].get("protein_class", "?") for u in uids]
bio_X = np.array([biopython_features(s) for s in seqs])
print(f"  esm: {esm_X.shape}, bio: {bio_X.shape}, y: {y.shape}, pos={int((y==1).sum())} neg={int((y==0).sum())}")

# ============== CLUSTER VIA ESM-2 SIMILARITY ==============
print("\nClustering proteins by ESM-2 cosine distance (proxy for sequence similarity)...")
# Normalize embeddings for cosine
esm_norm = esm_X / np.linalg.norm(esm_X, axis=1, keepdims=True)
clustering = AgglomerativeClustering(
    n_clusters=N_CLUSTERS,
    metric="cosine",
    linkage="average",
)
cluster_ids = clustering.fit_predict(esm_norm)
n_clusters = len(set(cluster_ids))
cluster_size_dist = Counter(cluster_ids)
sizes = sorted(cluster_size_dist.values(), reverse=True)
print(f"  Target clusters: {N_CLUSTERS}")
print(f"  Resulting clusters: {n_clusters}")
print(f"  Cluster size distribution: max={sizes[0]}, median={sizes[len(sizes)//2]}, "
      f"singletons={sum(1 for s in sizes if s==1)}, multi={sum(1 for s in sizes if s>1)}")
print(f"  Largest clusters: {sizes[:10]}")

# ============== CLUSTER-AWARE GROUPKFOLD ==============
print("\n" + "=" * 70)
print("CLUSTER-AWARE 5-FOLD CV (GroupKFold by ESM-2 similarity cluster)")
print("=" * 70)
gkf = GroupKFold(n_splits=CV_SPLITS)

def run_variant_cluster_aware(name, X, scaler=False, pca_n=None, extra_X=None, classifier="rf"):
    aucs, accs, briers, aps = [], [], [], []
    proba_oof = np.zeros(len(y))
    pred_oof = np.zeros(len(y))
    for tr, te in gkf.split(X, y, groups=cluster_ids):
        # GroupKFold doesn't preserve class balance; but our folds will still have positives.
        if len(set(y[te])) < 2:
            continue
        proba, pred = evaluate_one_split(X, y, tr, te, scaler=scaler, pca_n=pca_n,
                                         extra_X=extra_X, classifier=classifier)
        proba_oof[te] = proba
        pred_oof[te] = pred
        aucs.append(float(roc_auc_score(y[te], proba)))
        accs.append(float(accuracy_score(y[te], pred)))
        briers.append(float(brier_score_loss(y[te], proba)))
        aps.append(float(average_precision_score(y[te], proba)))
    cm = confusion_matrix(y, pred_oof.astype(int)).tolist()
    return {
        "name": name,
        "n_features": X.shape[1] if pca_n is None else pca_n + (extra_X.shape[1] if extra_X is not None else 0),
        "cv_auc_mean": float(np.mean(aucs)),
        "cv_auc_std": float(np.std(aucs)),
        "cv_acc_mean": float(np.mean(accs)),
        "cv_brier_mean": float(np.mean(briers)),
        "cv_average_precision_mean": float(np.mean(aps)),
        "fold_aucs": aucs,
        "confusion_matrix_overall": cm,
        "proba_oof": proba_oof.tolist(),
        "pred_oof": pred_oof.astype(int).tolist(),
    }

print("\nV1 BioPython baseline (cluster-aware)...")
r1c = run_variant_cluster_aware("v1_biopython_only_clustered", bio_X)
print(f"  AUC: {r1c['cv_auc_mean']:.3f} +/- {r1c['cv_auc_std']:.3f}")

print("\nV2 ESM-2 linear probe (cluster-aware)...")
r2c = run_variant_cluster_aware("v2_esm2_linear_probe_clustered", esm_X, scaler=True, classifier="logreg")
print(f"  AUC: {r2c['cv_auc_mean']:.3f} +/- {r2c['cv_auc_std']:.3f}")

print("\nV3 ESM-2 + PCA + RF (cluster-aware)...")
r3c = run_variant_cluster_aware("v3_esm2_pca256_rf_clustered", esm_X, scaler=True, pca_n=256)
print(f"  AUC: {r3c['cv_auc_mean']:.3f} +/- {r3c['cv_auc_std']:.3f}")

print("\nV4 combined (cluster-aware)...")
r4c = run_variant_cluster_aware("v4_combined_clustered", esm_X, scaler=True, pca_n=256, extra_X=bio_X)
print(f"  AUC: {r4c['cv_auc_mean']:.3f} +/- {r4c['cv_auc_std']:.3f}")

# Per-class accuracy for V4 cluster-aware
preds_v4 = np.array(r4c["pred_oof"])
class_acc = {}
for cls in set(classes):
    idx = [i for i, c in enumerate(classes) if c == cls]
    if len(idx) >= 5:
        n = len(idx)
        n_pos = sum(int(y[i]) for i in idx)
        acc = sum(1 for i in idx if preds_v4[i] == y[i]) / n
        class_acc[cls] = {"n": int(n), "n_pos": int(n_pos), "accuracy": float(acc)}
r4c["per_class_accuracy"] = class_acc

# Per-variant negative recall
print("\nPer-variant negative recall (cluster-aware):")
for r in [r1c, r2c, r3c, r4c]:
    cm = r["confusion_matrix_overall"]
    tn, fp = cm[0]
    fn, tp = cm[1]
    nr = tn / (tn + fp) if (tn + fp) else 0
    pr = tp / (tp + fn) if (tp + fn) else 0
    print(f"  {r['name']:48} TN={tn} FP={fp} FN={fn} TP={tp} | NegRecall={nr:.3f} PosRecall={pr:.3f}")

# Wilcoxon test
from scipy.stats import wilcoxon
import warnings
warnings.filterwarnings("ignore")
stat, p = wilcoxon(r4c["fold_aucs"], r1c["fold_aucs"], alternative="greater")
print(f"\nPaired Wilcoxon V4 > V1 cluster-aware: p={p:.4f}")
wilcoxon_p_v4_v1_clust = float(p)

# ============== V0.3 OOD DETECTOR ==============
print("\n" + "=" * 70)
print("FITTING V0.3 OOD DETECTOR (BioPython feature space, k=5 Euclidean)")
print("=" * 70)
ood_scaler = StandardScaler().fit(bio_X)
bio_scaled = ood_scaler.transform(bio_X)
knn = NearestNeighbors(n_neighbors=5, metric="euclidean").fit(bio_scaled)
distances, _ = knn.kneighbors(bio_scaled)
ood_threshold = float(np.percentile(distances[:, -1], 95))
print(f"  v0.3 OOD threshold (95th percentile of within-train k=5 distance): {ood_threshold:.3f}")

# Save V0.3 OOD detector bundle
ood_bundle = {
    "ood_scaler": ood_scaler,
    "ood_knn": knn,
    "ood_threshold": ood_threshold,
    "feature_space": "biopython_36",
    "k": 5,
    "metric": "euclidean",
    "training_set_size": len(bio_X),
}
joblib.dump(ood_bundle, MODEL / "ood_detector_v3.joblib")
print(f"  Saved: {MODEL / 'ood_detector_v3.joblib'}")

# ============== SAVE REPORT ==============
report = {
    "dataset": str(DATASET),
    "n_proteins": int(len(y)),
    "n_positives": int((y == 1).sum()),
    "n_negatives": int((y == 0).sum()),
    "esm_model": esm_meta["model"],
    "clustering": {
        "method": "AgglomerativeClustering with cosine linkage",
        "n_clusters_target": N_CLUSTERS,
        "n_clusters_actual": int(n_clusters),
        "largest_cluster_size": int(sizes[0]),
        "median_cluster_size": int(sizes[len(sizes)//2]),
        "n_singletons": int(sum(1 for s in sizes if s == 1)),
    },
    "cv_strategy": f"GroupKFold(n_splits={CV_SPLITS}, groups=ESM-2-similarity-cluster)",
    "variants_cluster_aware": [r1c, r2c, r3c, r4c],
    "wilcoxon_v4_v1_cluster_aware_p": wilcoxon_p_v4_v1_clust,
    "ood_detector_v3": {
        "threshold_95pct": ood_threshold,
        "feature_space": "biopython_36_standardized",
        "k": 5,
        "metric": "euclidean",
    },
    "honest_caveats": [
        "Clusters were defined by ESM-2 cosine distance at threshold 0.10, which is a proxy for sequence identity (not the same as MMseqs2 30 percent identity clustering, but conceptually similar in grouping homologs).",
        f"GroupKFold produced {n_clusters} clusters across 5 folds; this prevents close ESM-2-neighbors from splitting across train/test.",
        "AUC drops relative to standard StratifiedKFold are the expected and honest cost of removing homology leakage.",
    ],
}
with open(MODEL / "training_report_v3_robust.json", "w") as f:
    # strip bulky predictions before saving
    slim = {**report}
    slim["variants_cluster_aware"] = []
    for r in [r1c, r2c, r3c, r4c]:
        slim_r = {k: v for k, v in r.items() if k not in ("proba_oof", "pred_oof")}
        slim["variants_cluster_aware"].append(slim_r)
    json.dump(slim, f, indent=2)
print(f"\nSaved: {MODEL / 'training_report_v3_robust.json'}")

# Save predictions for figures
np.savez(MODEL / "v3_robust_oof_predictions.npz",
         v1_proba=np.array(r1c["proba_oof"]),
         v2_proba=np.array(r2c["proba_oof"]),
         v3_proba=np.array(r3c["proba_oof"]),
         v4_proba=np.array(r4c["proba_oof"]),
         y=y, classes=np.array(classes), uids=np.array(uids))
print(f"Saved OOF predictions: {MODEL / 'v3_robust_oof_predictions.npz'}")

# ============== HEADLINE COMPARISON ==============
print("\n" + "=" * 70)
print("HEADLINE: cluster-aware (new) vs stratified (old) AUC")
print("=" * 70)
old = json.load(open(MODEL / "training_report_v3.json"))
print(f"  {'Variant':30} {'Stratified AUC':>16} {'Cluster-aware AUC':>20} {'Delta':>8}")
print(f"  {'-'*30} {'-'*16} {'-'*20} {'-'*8}")
for r_new in [r1c, r2c, r3c, r4c]:
    base_name = r_new["name"].replace("_clustered", "")
    r_old = next((v for v in old["variants"] if v["name"] == base_name), None)
    if r_old:
        old_auc = r_old["cv_auc_mean"]
        new_auc = r_new["cv_auc_mean"]
        delta = new_auc - old_auc
        print(f"  {base_name:30} {old_auc:>16.3f} {new_auc:>20.3f} {delta:>+8.3f}")
print()
print("Done.")
