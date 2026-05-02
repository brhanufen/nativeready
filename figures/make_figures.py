"""Generate Figures 1, 2, and 3 for the NativeReady preprint."""
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import joblib
import torch
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, brier_score_loss
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from Bio.SeqUtils.ProtParam import ProteinAnalysis

ROOT = Path("/Users/bfentaw2/startup/nativeready")
DATA = ROOT / "data"
MODEL = ROOT / "model"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

DATASET = DATA / "dataset_combined_v4_2026-05-02.json"
ESM_NPY = MODEL / "esm2_embeddings_634.npy"
ESM_META = MODEL / "esm2_embeddings_metadata.json"

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
RANDOM_STATE = 42
CV_SPLITS = 5

CYAN = "#0A6B68"
ACCENT = "#5CF2EA"
GRAY = "#666666"
RED = "#C0392B"
GOOD = "#2A8C4A"
WARN = "#D9892A"


def clean_seq(seq):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    cleaned = []
    for c in seq:
        if c in STANDARD_AA:
            cleaned.append(c)
        elif c in "BZJU":
            cleaned.append({"B":"N","Z":"Q","J":"L","U":"C"}[c])
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
    return feats


# ============== LOAD DATA ==============
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
print(f"  esm: {esm_X.shape}, bio: {bio_X.shape}, y: {y.shape}")


# ============== FIGURE 2: ROC + Calibration ==============
print("\nFigure 2: ROC and calibration curves")

skf = StratifiedKFold(n_splits=CV_SPLITS, shuffle=True, random_state=RANDOM_STATE)

variants = []  # (label, color, fold_y, fold_proba)

# V1 BioPython only
def make_rf():
    base = RandomForestClassifier(n_estimators=300, max_depth=10, min_samples_leaf=2,
                                  class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1)
    return CalibratedClassifierCV(base, method="isotonic", cv=3)

def make_logreg():
    base = LogisticRegression(C=1.0, penalty="l2", max_iter=2000, solver="lbfgs",
                               class_weight="balanced", random_state=RANDOM_STATE)
    return CalibratedClassifierCV(base, method="isotonic", cv=3)

def cv_predictions(X_full, y, make_clf, scaler=False, pca_n=None, extra_X=None):
    """Return (y_oof, proba_oof) stitched across folds."""
    proba_oof = np.zeros(len(y))
    for tr, te in skf.split(X_full, y):
        if scaler:
            sc = StandardScaler().fit(X_full[tr])
            Xs_tr = sc.transform(X_full[tr])
            Xs_te = sc.transform(X_full[te])
        else:
            Xs_tr = X_full[tr]
            Xs_te = X_full[te]
        if pca_n is not None:
            pca = PCA(n_components=pca_n, random_state=RANDOM_STATE).fit(Xs_tr)
            Xp_tr = pca.transform(Xs_tr)
            Xp_te = pca.transform(Xs_te)
        else:
            Xp_tr, Xp_te = Xs_tr, Xs_te
        if extra_X is not None:
            esc = StandardScaler().fit(extra_X[tr])
            extra_tr = esc.transform(extra_X[tr])
            extra_te = esc.transform(extra_X[te])
            Xp_tr = np.hstack([Xp_tr, extra_tr])
            Xp_te = np.hstack([Xp_te, extra_te])
        clf = make_clf()
        clf.fit(Xp_tr, y[tr])
        proba_oof[te] = clf.predict_proba(Xp_te)[:, 1]
    return proba_oof

print("  Loading cluster-aware OOF probabilities...")
preds = np.load(MODEL / "v3_robust_oof_predictions.npz", allow_pickle=True)
proba_v1 = preds["v1_proba"]
proba_v2 = preds["v2_proba"]
proba_v3 = preds["v3_proba"]
proba_v4 = preds["v4_proba"]
print("    Loaded all four variants from saved cluster-aware run")

variant_data = [
    ("V1 BioPython only (baseline)",  CYAN,   proba_v1),
    ("V2 ESM-2 linear probe",         WARN,   proba_v2),
    ("V3 ESM-2 PCA-256 + RF",         GOOD,   proba_v3),
    ("V4 combined (winner)",          RED,    proba_v4),
]

# 2-panel figure: ROC | calibration
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), dpi=150)

# Panel A: ROC
ax = axes[0]
for label, color, proba in variant_data:
    fpr, tpr, _ = roc_curve(y, proba)
    auc_v = auc(fpr, tpr)
    ax.plot(fpr, tpr, color=color, lw=2, label=f"{label} (AUC = {auc_v:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("(A) ROC curves on stitched out-of-fold predictions")
ax.legend(loc="lower right", fontsize=8, frameon=False)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.grid(True, alpha=0.2)

# Panel B: Calibration
ax = axes[1]
for label, color, proba in variant_data:
    bs = brier_score_loss(y, proba)
    frac_pos, mean_pred = calibration_curve(y, proba, n_bins=10, strategy="quantile")
    ax.plot(mean_pred, frac_pos, "o-", color=color, lw=1.5, ms=5,
            label=f"{label} (Brier = {bs:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5, label="Perfect calibration")
ax.set_xlabel("Mean predicted probability")
ax.set_ylabel("Fraction of positives")
ax.set_title("(B) Calibration curves (10 quantile bins)")
ax.legend(loc="upper left", fontsize=8, frameon=False)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig(FIGS / "Figure2_ROC_calibration.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGS / 'Figure2_ROC_calibration.png'}")


# ============== FIGURE 3: Per-class accuracy bar chart ==============
print("\nFigure 3: Per-class accuracy")

# Use the V4 predictions (winner)
preds_v4 = (proba_v4 >= 0.5).astype(int)

# Compute per-class accuracy
from collections import defaultdict
class_acc = {}
for cls in set(classes):
    idx = [i for i, c in enumerate(classes) if c == cls]
    if len(idx) >= 5:
        n = len(idx)
        n_pos = sum(y[i] for i in idx)
        acc = sum(1 for i in idx if preds_v4[i] == y[i]) / n
        class_acc[cls] = (n, int(n_pos), acc)

# Sort: positives first by n desc, then negatives by n desc
sorted_classes = sorted(class_acc.items(),
                        key=lambda x: (-(x[1][1] / x[1][0]), -x[1][0]))

fig, ax = plt.subplots(figsize=(11, 5), dpi=150)
labels = [c if len(c) <= 24 else c[:22]+"..." for c, _ in sorted_classes]
ns = [v[0] for _, v in sorted_classes]
accs = [v[2] for _, v in sorted_classes]
n_pos = [v[1] for _, v in sorted_classes]
pos_pct = [100*p/n for p, n in zip(n_pos, ns)]

# Color: green if ~all positive, red if ~all negative, gray mixed
colors = []
for p, n in zip(n_pos, ns):
    pct = p/n
    if pct >= 0.95:
        colors.append(GOOD)
    elif pct <= 0.05:
        colors.append(RED)
    else:
        colors.append(GRAY)

bars = ax.bar(range(len(labels)), accs, color=colors, edgecolor="black", linewidth=0.6)
# Annotate n on top
for i, (bar, n, pct) in enumerate(zip(bars, ns, pos_pct)):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 0.01,
            f"n={n}\n({pct:.0f}% pos)", ha='center', va='bottom', fontsize=7)

ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Accuracy (out-of-fold, V4 combined model)")
ax.set_ylim(0, 1.18)
ax.set_title("Per-class accuracy. Green: positive-dominated. Red: negative-dominated. Gray: mixed.")
ax.axhline(0.5, color="gray", lw=0.5, linestyle="--", alpha=0.5)
ax.grid(True, axis="y", alpha=0.2)
plt.tight_layout()
plt.savefig(FIGS / "Figure3_per_class_accuracy.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGS / 'Figure3_per_class_accuracy.png'}")


# ============== FIGURE 1: Pipeline schematic ==============
print("\nFigure 1: Pipeline schematic")

fig, ax = plt.subplots(figsize=(11, 6), dpi=150)
ax.set_xlim(0, 10)
ax.set_ylim(0, 7)
ax.axis("off")

def box(x, y, w, h, label, sub=None, color=ACCENT, edge=CYAN, text_color="#08272A"):
    rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.15",
                           facecolor=color, edgecolor=edge, linewidth=1.5)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2 + (0.18 if sub else 0), label,
            ha="center", va="center", fontsize=10, weight="bold", color=text_color)
    if sub:
        ax.text(x + w/2, y + h/2 - 0.22, sub,
                ha="center", va="center", fontsize=8, color=text_color, style="italic")

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                  arrowstyle="->", mutation_scale=15,
                                  color=GRAY, lw=1.2))

# STAGE 1: Data sources (left column)
ax.text(1.4, 6.55, "Data sources", ha="center", fontsize=11, weight="bold", color=CYAN)
box(0.2, 5.4, 2.4, 0.9, "Hand-curated base", "n = 232 (literature)", color="#E8F8F7")
box(0.2, 4.2, 2.4, 0.9, "RCSB PDB", "n = 358 (full-text query)", color="#E8F8F7")
box(0.2, 3.0, 2.4, 0.9, "EuropePMC supplements", "n = 44 (evidence-based)", color="#E8F8F7")
box(0.2, 1.5, 2.4, 1.1, "Dedup by UniProt", "634 unique proteins\n(538 pos, 96 neg)", color=ACCENT)
arrow(1.4, 5.4, 1.4, 2.65)
arrow(1.4, 4.2, 1.4, 2.65)
arrow(1.4, 3.0, 1.4, 2.65)

# STAGE 2: Features
ax.text(4.4, 6.55, "Features", ha="center", fontsize=11, weight="bold", color=CYAN)
box(3.2, 4.7, 2.4, 1.0, "BioPython 36-dim", "physicochemical\n(MW, pI, GRAVY, AA%)", color="#E8F8F7")
box(3.2, 3.2, 2.4, 1.0, "ESM-2 1280-dim", "facebook/esm2_t33_650M\nmean-pooled", color="#E8F8F7")
arrow(2.6, 2.0, 3.2, 2.0)
arrow(2.6, 2.0, 3.2, 5.2)
arrow(2.6, 2.0, 3.2, 3.7)

# STAGE 3: Models
ax.text(7.2, 6.55, "Models (5-fold CV)", ha="center", fontsize=11, weight="bold", color=CYAN)
box(6.0, 5.4, 2.5, 0.6, "V1 BioPython + RF", color="#FFF6E6")
box(6.0, 4.6, 2.5, 0.6, "V2 ESM-2 LinearProbe", color="#FFF6E6")
box(6.0, 3.8, 2.5, 0.6, "V3 ESM-2 PCA + RF", color="#FFF6E6")
box(6.0, 3.0, 2.5, 0.6, "V4 Combined (winner)", color=ACCENT)
arrow(5.6, 5.2, 6.0, 5.7)
arrow(5.6, 5.2, 6.0, 4.9)
arrow(5.6, 3.7, 6.0, 4.1)
arrow(5.6, 3.7, 6.0, 3.3)
arrow(5.6, 5.2, 6.0, 3.3)

# STAGE 4: Output
ax.text(9.0, 6.55, "Outputs", ha="center", fontsize=11, weight="bold", color=CYAN)
box(8.7, 4.7, 1.2, 1.0, "Open\ndataset", "CC-BY 4.0\nZenodo", color="#E8F8F7")
box(8.7, 3.2, 1.2, 1.0, "Web tool", "nativeready\n.netlify.app", color=ACCENT)
arrow(8.5, 5.2, 8.7, 5.2)
arrow(8.5, 3.7, 8.7, 3.7)

# Title
ax.text(5.0, 0.5, "NativeReady pipeline: dataset construction, features, four model variants under identical 5-fold CV, and public deployment",
        ha="center", fontsize=10, style="italic", color=GRAY)

plt.tight_layout()
plt.savefig(FIGS / "Figure1_pipeline.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {FIGS / 'Figure1_pipeline.png'}")

print("\nAll figures generated.")
print("Files:")
for f in sorted(FIGS.glob("*.png")):
    print(f"  {f.name}: {f.stat().st_size:,} bytes")
