"""Build esm2_embeddings_635.npy from per-protein cache, in v7 dataset order."""
import json
from pathlib import Path
import numpy as np

DATA = Path("/Users/bfentaw2/startup/nativeready/data")
MODEL = Path("/Users/bfentaw2/startup/nativeready/model")
EMB_DIR = MODEL / "esm2_embeddings"

with open(DATA / "dataset_combined_v7_2026-05-11.json") as f:
    records = json.load(f)
print(f"v7 records: {len(records)}")

mat = []
uids = []
extra = []
n_truncated = 0
truncated_uids = []
for p in records:
    uid = p["uniprot_id"]
    cache = EMB_DIR / f"{uid}.npy"
    if not cache.exists():
        raise SystemExit(f"Missing embedding cache for {uid}")
    emb = np.load(cache)
    mat.append(emb)
    uids.append(uid)
    extra.append({
        "uniprot_id": uid,
        "label": p["label"],
        "sequence_length": p.get("sequence_length"),
        "protein_class": p.get("protein_class"),
        "mw_kda": p.get("mw_kda"),
    })
    if p.get("sequence_length", 0) > 1022:
        n_truncated += 1
        truncated_uids.append(uid)

mat = np.stack(mat)
print(f"matrix shape: {mat.shape}")
assert mat.shape == (635, 1280), f"unexpected shape: {mat.shape}"

np.save(MODEL / "esm2_embeddings_635.npy", mat)
with open(MODEL / "esm2_embeddings_metadata.json", "w") as f:
    json.dump({
        "uids": uids,
        "extra": extra,
        "embedding_dim": int(mat.shape[1]),
        "n_proteins": int(mat.shape[0]),
        "model": "facebook/esm2_t33_650M_UR50D",
        "max_seq_len": 1022,
        "n_truncated": n_truncated,
        "truncated_uids": truncated_uids,
    }, f, indent=2)
print(f"saved: {MODEL / 'esm2_embeddings_635.npy'}")
print(f"saved: {MODEL / 'esm2_embeddings_metadata.json'}")
