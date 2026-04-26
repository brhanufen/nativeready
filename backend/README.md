# NativeReady backend

FastAPI service that predicts whether a protein sequence is likely to yield
usable native mass spectrometry (native MS) data.

## Run locally

```bash
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
```

Then open http://localhost:8000/docs for the auto-generated Swagger UI.

## Endpoints

- `GET /` — health check
- `POST /predict` — body: `{"sequence": "MGSSHHHHH..."}`
- `GET /docs` — Swagger UI

## Test the predict endpoint

```bash
curl -X POST http://localhost:8000/predict \
  -H "Content-Type: application/json" \
  -d '{"sequence": "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDPMVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"}'
```

## Model

The trained model is loaded from `/Users/bfentaw2/startup/nativeready/model/model.joblib`.
If absent, the API falls back to a transparent heuristic predictor and labels
its responses with `model_version: "0.1-baseline-heuristic"` so callers know
they are not seeing trained-model output.

## Docker

```bash
docker build -t nativeready-backend .
docker run -p 8000:8000 nativeready-backend
```
