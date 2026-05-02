#!/usr/bin/env python3
"""Compute ESM-2 embeddings for all proteins in dataset_combined_v4_2026-05-02.json.

Design decisions (with reasoning):
- Model: facebook/esm2_t33_650M_UR50D. The 650M variant is the standard for
  downstream tasks and fits comfortably in 16 GB RAM. Going larger (3B, 15B)
  is overkill for 634 samples and adds unjustified compute cost.
- Device: MPS (Apple Silicon GPU) when available, fall back to CPU. ~3-5x speedup.
- Sequence truncation: cap at 1022 residues (model max = 1024 incl. [CLS] and [EOS]).
  About 95% of proteins fit; the few that don't (mostly viral polyproteins) get
  the first 1022 residues. Documented as a known limitation in the report.
- Pooling: mean over per-residue embeddings (excluding special tokens). Standard,
  defensible, what most ESM-2 downstream tasks use.
- Caching: per-protein .npy files so we can resume on interruption.
- Precision: float32 storage (3.2 MB total for 634 x 1280). Inference uses default
  precision per the model card.
"""

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel

DATA = Path("/Users/bfentaw2/startup/nativeready/data")
MODEL_DIR = Path("/Users/bfentaw2/startup/nativeready/model")
EMB_DIR = MODEL_DIR / "esm2_embeddings"
EMB_DIR.mkdir(exist_ok=True)

DATASET = DATA / "dataset_combined_v4_2026-05-02.json"
ESM_MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
MAX_LEN = 1022  # 1024 model max minus 2 special tokens

print(f"Loading dataset: {DATASET}")
with open(DATASET) as f:
    proteins = json.load(f)
print(f"  Loaded {len(proteins)} proteins")

print(f"\nLoading ESM-2 model: {ESM_MODEL_NAME}")
print(f"  (first run downloads ~2.5 GB; cached after)")
tokenizer = AutoTokenizer.from_pretrained(ESM_MODEL_NAME)
model = AutoModel.from_pretrained(ESM_MODEL_NAME)
model.eval()

# Device selection
if torch.backends.mps.is_available():
    device = torch.device("mps")
    print(f"  Using MPS (Apple Silicon GPU)")
elif torch.cuda.is_available():
    device = torch.device("cuda")
    print(f"  Using CUDA")
else:
    device = torch.device("cpu")
    print(f"  Using CPU (slower)")
model = model.to(device)

# Quick check what the embedding dim is
with torch.no_grad():
    test_in = tokenizer("MKQH", return_tensors="pt").to(device)
    test_out = model(**test_in)
    emb_dim = test_out.last_hidden_state.shape[-1]
    print(f"  Embedding dim: {emb_dim}")

print(f"\nProcessing {len(proteins)} proteins")
print(f"  Truncation: sequences > {MAX_LEN} aa truncated to first {MAX_LEN}")
print(f"  Cache dir: {EMB_DIR}")

start = time.time()
n_done = 0
n_skipped_existing = 0
n_truncated = 0
n_errors = 0
errors = []
truncated_uids = []

for i, p in enumerate(proteins):
    uid = p["uniprot_id"]
    seq = p["sequence"].upper()
    cache_path = EMB_DIR / f"{uid}.npy"

    if cache_path.exists():
        n_skipped_existing += 1
        continue

    orig_len = len(seq)
    if orig_len > MAX_LEN:
        seq = seq[:MAX_LEN]
        n_truncated += 1
        truncated_uids.append((uid, orig_len))

    try:
        with torch.no_grad():
            inputs = tokenizer(seq, return_tensors="pt", truncation=True, max_length=MAX_LEN + 2).to(device)
            outputs = model(**inputs)
            # last_hidden_state shape: [1, seq_len_with_special_tokens, embed_dim]
            # Remove special tokens (positions 0 and -1) before mean-pooling
            hidden = outputs.last_hidden_state[0, 1:-1, :]  # [actual_len, embed_dim]
            embedding = hidden.mean(dim=0).cpu().numpy().astype(np.float32)
        np.save(cache_path, embedding)
        n_done += 1
    except Exception as e:
        n_errors += 1
        errors.append((uid, str(e)[:200]))
        continue

    if (i + 1) % 25 == 0:
        elapsed = time.time() - start
        rate = (n_done + 1) / max(elapsed, 0.1)
        remaining = (len(proteins) - (i + 1)) / max(rate, 0.001)
        print(f"  [{i+1}/{len(proteins)}] done={n_done} skipped={n_skipped_existing} truncated={n_truncated} errors={n_errors} | {rate:.1f}/s; ~{remaining:.0f}s remaining")

elapsed = time.time() - start
print("\n" + "=" * 70)
print("ESM-2 EMBEDDING SUMMARY")
print("=" * 70)
print(f"  Total proteins: {len(proteins)}")
print(f"  Newly embedded: {n_done}")
print(f"  Skipped (cached): {n_skipped_existing}")
print(f"  Truncated to {MAX_LEN} aa: {n_truncated}")
print(f"  Errors: {n_errors}")
print(f"  Elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")
print(f"  Embeddings saved to: {EMB_DIR}")

if truncated_uids:
    print(f"\n  Truncated proteins (UniProt: original_length):")
    for uid, ln in truncated_uids[:10]:
        print(f"    {uid}: {ln}")
    if len(truncated_uids) > 10:
        print(f"    ... and {len(truncated_uids) - 10} more")

if errors:
    print(f"\n  Errors:")
    for uid, err in errors[:10]:
        print(f"    {uid}: {err}")

# Build a single combined matrix for downstream training
print("\nBuilding combined embedding matrix...")
emb_matrix = []
labels = []
uids = []
extra = []  # store sequence_length, protein_class, label, mw_kda for later

for p in proteins:
    uid = p["uniprot_id"]
    cache_path = EMB_DIR / f"{uid}.npy"
    if not cache_path.exists():
        continue
    emb = np.load(cache_path)
    emb_matrix.append(emb)
    labels.append(p["label"])
    uids.append(uid)
    extra.append({
        "uniprot_id": uid,
        "label": p["label"],
        "sequence_length": p.get("sequence_length"),
        "protein_class": p.get("protein_class"),
        "mw_kda": p.get("mw_kda"),
    })

emb_matrix = np.stack(emb_matrix)
labels = np.array(labels)
print(f"  Combined matrix shape: {emb_matrix.shape}")
print(f"  Label counts: positives={int((labels==1).sum())}  negatives={int((labels==0).sum())}")

np.save(MODEL_DIR / "esm2_embeddings_634.npy", emb_matrix)
with open(MODEL_DIR / "esm2_embeddings_metadata.json", "w") as f:
    json.dump({
        "uids": uids,
        "extra": extra,
        "embedding_dim": int(emb_matrix.shape[1]),
        "n_proteins": int(emb_matrix.shape[0]),
        "model": ESM_MODEL_NAME,
        "max_seq_len": MAX_LEN,
        "n_truncated": n_truncated,
        "truncated_uids": [u for u, _ in truncated_uids],
    }, f, indent=2)

print(f"  Saved: {MODEL_DIR / 'esm2_embeddings_634.npy'}")
print(f"  Saved: {MODEL_DIR / 'esm2_embeddings_metadata.json'}")
