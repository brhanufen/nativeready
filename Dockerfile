# Root-level Dockerfile for Railway deployment.
# Builds the backend FastAPI service with the model files included.
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# System deps for numpy/scikit-learn wheels
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Python deps first (better layer caching)
COPY backend/requirements.txt /app/backend/requirements.txt
RUN pip install --upgrade pip && pip install -r /app/backend/requirements.txt

# Copy application code and model artifacts
COPY backend/ /app/backend/
COPY model/   /app/model/

# Default port (Railway overrides via $PORT)
ENV PORT=8000
EXPOSE 8000

WORKDIR /app/backend

# Exec form via sh -c so $PORT expands at runtime AND signals are forwarded properly.
CMD ["sh", "-c", "uvicorn main:app --host 0.0.0.0 --port ${PORT:-8000}"]
