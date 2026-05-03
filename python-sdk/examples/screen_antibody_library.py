"""Real-world example: screen a small antibody-fragment library and rank by suitability.

Useful pattern for a biotech analytical-development team that wants to triage
candidates before committing instrument time.
"""
import csv
from nativeready import Client


# Small example library of antibody/protein fragments to screen.
# In practice you would load this from a CSV or from your internal sequence database.
LIBRARY = [
    ("Fab_v1", "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCAR..."),
    ("Fab_v2", "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCAR..."),
    ("scFv_A", "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNNKNYLAWYQQKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQ..."),
    # Real implementations will have full sequences (>= 50 aa)
    ("ubiquitin_control",
     "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"),
]

# Filter to entries that meet the minimum sequence length for the model
LIBRARY = [(name, seq) for name, seq in LIBRARY if len(seq) >= 50 and "..." not in seq]

print(f"Screening {len(LIBRARY)} sequences...")
client = Client()
records = [{"id": name, "sequence": seq} for name, seq in LIBRARY]
results = client.predict_batch(records)

# Sort by score, highest first
ranked = sorted(zip(LIBRARY, results), key=lambda x: -x[1].score)

print(f"\n{'rank':<5}{'name':<22}{'score':<7}{'label':<14}{'OOD'}")
print("-" * 55)
for rank, ((name, _), r) in enumerate(ranked, 1):
    print(f"{rank:<5}{name:<22}{r.score:<7}{r.label:<14}{'Y' if r.is_ood else ''}")

# Save full results for archival / audit
with open("screening_results.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["rank", "name", "score", "label", "ci_lower", "ci_upper",
                "is_ood", "model_version", "recommendations"])
    for rank, ((name, _), r) in enumerate(ranked, 1):
        w.writerow([rank, name, r.score, r.label,
                    r.confidence_lower, r.confidence_upper,
                    "Y" if r.is_ood else "",
                    r.model_version,
                    " | ".join(r.recommendations)])
print("\nResults written to screening_results.csv")
