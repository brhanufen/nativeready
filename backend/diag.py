"""Diagnostic — run inside Railway container to check model loading."""
import os
import sys
import traceback

print("=== System ===")
import sklearn, joblib
print(f"sklearn {sklearn.__version__}")
print(f"joblib  {joblib.__version__}")

print("\n=== Files ===")
for p in [
    "/app/model/model_v2.joblib",
    "/app/model/model.joblib",
    "/app/model/feature_names.json",
]:
    if os.path.exists(p):
        print(f"  EXISTS  {p}  ({os.path.getsize(p):,} bytes)")
    else:
        print(f"  MISSING {p}")

print("\n=== Predictor ===")
sys.path.insert(0, "/app/backend")
try:
    import predictor_v2
    print(f"  MODEL_PATH_V2 = {predictor_v2.MODEL_PATH_V2}")
    ok = predictor_v2._try_load()
    print(f"  _try_load() = {ok}")
    print(f"  _version    = {predictor_v2._version}")
    print(f"  _bundle is None?: {predictor_v2._bundle is None}")
    if predictor_v2._bundle is not None:
        print(f"  _bundle type: {type(predictor_v2._bundle).__name__}")
        if isinstance(predictor_v2._bundle, dict):
            print(f"  _bundle keys: {list(predictor_v2._bundle.keys())}")
            for k, v in predictor_v2._bundle.items():
                print(f"    {k}: {type(v).__name__}")
    print(f"  _feature_names is None?: {predictor_v2._feature_names is None}")

    print("\n=== Prediction test ===")
    seq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    result = predictor_v2.predict(seq)
    print(f"  Score: {result['suitability_score']}, Label: {result['suitability_label']}")
except Exception as e:
    print(f"\n!!! EXCEPTION: {e}")
    traceback.print_exc()
