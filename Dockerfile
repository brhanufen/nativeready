# Root-level Dockerfile for Railway deployment of NativeReady v0.3.
# v0.3 includes ESM-2 inference; image is ~4-5 GB after weights are baked in.
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    HF_HOME=/app/.cache/huggingface \
    TRANSFORMERS_OFFLINE=0

WORKDIR /app

# System deps for numpy/scikit-learn/torch wheels
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Python deps first (better layer caching). CPU-only torch saves ~1 GB vs CUDA build.
COPY backend/requirements.txt /app/backend/requirements.txt
RUN pip install --upgrade pip && pip install -r /app/backend/requirements.txt

# Pre-download the ESM-2 weights into the image so first-request latency is acceptable.
# Without this, the first user request would trigger a ~2.5 GB download from HuggingFace.
RUN python -c "from transformers import AutoTokenizer, AutoModel; \
    AutoTokenizer.from_pretrained('facebook/esm2_t33_650M_UR50D'); \
    AutoModel.from_pretrained('facebook/esm2_t33_650M_UR50D'); \
    print('ESM-2 weights cached')"

# Copy application code and model artifacts
COPY backend/ /app/backend/
COPY model/   /app/model/

# Default port (Railway overrides via $PORT)
ENV PORT=8000
EXPOSE 8000

WORKDIR /app/backend

# Exec form via sh -c so $PORT expands at runtime AND signals are forwarded properly.
CMD ["sh", "-c", "uvicorn main:app --host 0.0.0.0 --port ${PORT:-8000}"]
