"""Tests for the Method 1 heterogeneity forward-simulator.

Run either way:
    python3 tests/test_heterogeneity.py
    python3 -m pytest tests/test_heterogeneity.py -v

All protein sequences come from data/dataset_combined_v7_2026-05-11.json -- real
UniProt records, no synthetic sequences. The mass constants are re-derived from
atomic monoisotopic masses inside the tests so the ladder cannot silently drift.
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "backend"))

from features import extract_features  # noqa: E402
from heterogeneity import (  # noqa: E402
    DEFAULT_DETECTION_THRESHOLD,
    DEFAULT_N_GLYCAN_DIST,
    DEFAULT_SEQUON_OCCUPANCY,
    GLYCAN_MASSES,
    N_GLYCAN_COMPOSITIONS,
    detectable_species,
    glycan_mass,
    heterogeneity_report,
    minimum_resolution_required,
    n_glyc_sites,
    resolvability,
    simulate_dar_envelope,
    simulate_glycoform_envelope,
)

DATASET = ROOT / "data" / "dataset_combined_v7_2026-05-11.json"
_RECORDS = {r["uniprot_id"]: r for r in json.loads(DATASET.read_text())}


def seq(uniprot_id):
    return _RECORDS[uniprot_id]["sequence"]


# Real records used below:
#   P01267  bovine thyroglobulin  - heavily glycosylated; v7 records this as a
#           documented real failure with failure_mode uninterpretable_heterogeneity
#   P02787  human transferrin     - glycoprotein routinely run by native MS
#   P01857  human IgG1 heavy chain constant - antibody scaffold, used for the ADC case
#   P0CG48  polyubiquitin-C       - non-glycosylated
#   P00918  carbonic anhydrase 2  - non-glycosylated native-MS calibrant
THYROGLOBULIN = "P01267"
TRANSFERRIN = "P02787"
IGG1_HC = "P01857"
UBIQUITIN = "P0CG48"
CARBONIC_ANHYDRASE = "P00918"

# Nominal MC-vc-PAB-MMAE linker-payload mass. Construct-specific and supplied by
# the caller by design -- the module deliberately does not default a drug mass.
VCMMAE_NOMINAL_DA = 1316.6


# ---------------------------------------------------------------------------
# Mass ladder
# ---------------------------------------------------------------------------

def test_glycan_residue_masses_match_atomic_derivation():
    """Guard the published ladder against silent edits."""
    atom = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.9949146196}

    def mass(**formula):
        return sum(atom[e] * n for e, n in formula.items())

    expected = {
        "Hex": mass(C=6, H=10, O=5),
        "HexNAc": mass(C=8, H=13, N=1, O=5),
        "Fuc": mass(C=6, H=10, O=4),
        "NeuAc": mass(C=11, H=17, N=1, O=8),
        "NeuGc": mass(C=11, H=17, N=1, O=9),
        "Pent": mass(C=5, H=8, O=4),
    }
    for name, value in expected.items():
        assert abs(GLYCAN_MASSES[name] - value) < 5e-4, name


def test_named_glycan_compositions():
    assert abs(glycan_mass(N_GLYCAN_COMPOSITIONS["G0F"]) - 1444.5339) < 1e-3
    assert abs(glycan_mass(N_GLYCAN_COMPOSITIONS["G1F"]) - 1606.5867) < 1e-3
    assert abs(glycan_mass(N_GLYCAN_COMPOSITIONS["G2F"]) - 1768.6395) < 1e-3
    assert abs(glycan_mass(N_GLYCAN_COMPOSITIONS["Man5"]) - 1216.4229) < 1e-3


def test_isobaric_glycoforms_are_exactly_degenerate():
    """G0F + G2F and G1F + G1F have identical mass, so merging them is physics."""
    g0f = glycan_mass(N_GLYCAN_COMPOSITIONS["G0F"])
    g1f = glycan_mass(N_GLYCAN_COMPOSITIONS["G1F"])
    g2f = glycan_mass(N_GLYCAN_COMPOSITIONS["G2F"])
    assert abs((g0f + g2f) - 2 * g1f) < 1e-9


def test_default_distribution_is_normalisable():
    assert abs(sum(DEFAULT_N_GLYCAN_DIST.values()) - 1.0) < 1e-9
    assert all(v > 0 for v in DEFAULT_N_GLYCAN_DIST.values())
    for name in DEFAULT_N_GLYCAN_DIST:
        assert name in N_GLYCAN_COMPOSITIONS


# ---------------------------------------------------------------------------
# Site counting reuses features.py
# ---------------------------------------------------------------------------

def test_n_glyc_sites_delegates_to_features_module():
    """The simulator must count sites exactly as the shipped model does."""
    for uid in (THYROGLOBULIN, TRANSFERRIN, IGG1_HC, UBIQUITIN):
        s = seq(uid)
        assert n_glyc_sites(s) == int(
            extract_features(s)["n_glyc_sequon_count"]
        ), uid
        assert n_glyc_sites(s, strict=True) == int(
            extract_features(s)["n_glyc_sequon_strict_count"]
        ), uid


def test_n_glyc_sites_short_sequence_does_not_raise():
    assert n_glyc_sites("MA") == 0


# ---------------------------------------------------------------------------
# Case 1: a known glycoprotein
# ---------------------------------------------------------------------------

def test_glycoprotein_thyroglobulin_is_unresolvable():
    """Heavily glycosylated: physics should call it unresolvable.

    v7 independently records this protein as a real documented failure with
    failure_mode 'uninterpretable_heterogeneity'. The simulator never sees that
    label -- it reaches the same conclusion from sequon count and glycan masses.
    """
    rep = heterogeneity_report(seq(THYROGLOBULIN))
    assert rep["applies"] is True
    assert rep["mode"] == "glycan"
    assert rep["n_sites"] == 14
    assert rep["computed"] is True
    assert rep["resolved_at_resolution"] is False
    assert rep["level"] == "high"
    assert rep["predicted_envelope_width_da"] > 1000.0
    assert rep["n_proteoforms"] > 100
    assert rep["min_resolution_required"] > 30000
    assert "unresolvable" in rep["reason"]


def test_glycoprotein_transferrin_resolves_at_orbitrap_resolution():
    """Two sequons on a 77 kDa protein: glycoforms separate at R = 30,000."""
    rep = heterogeneity_report(seq(TRANSFERRIN))
    assert rep["applies"] is True
    assert rep["n_sites"] == 2
    assert rep["resolved_at_resolution"] is True
    assert rep["level"] == "low"
    assert rep["min_resolution_required"] < 30000


def test_glycoform_envelope_is_normalised_and_ordered():
    env = simulate_glycoform_envelope(77000.0, 2)
    assert len(env) > 1
    assert abs(sum(a for _, a in env) - 1.0) < 1e-9
    assert env == sorted(env)
    assert all(m >= 77000.0 for m, _ in env)


def test_zero_sites_returns_single_species():
    env = simulate_glycoform_envelope(50000.0, 0)
    assert env == [(50000.0, 1.0)]


def test_more_sites_widens_the_envelope():
    widths = []
    for n in (1, 2, 4, 8):
        env = simulate_glycoform_envelope(100000.0, n)
        widths.append(env[-1][0] - env[0][0])
    assert widths == sorted(widths)
    assert widths[0] < widths[-1]


# ---------------------------------------------------------------------------
# Case 2: an ADC
# ---------------------------------------------------------------------------

def test_dar_envelope_is_binomial_with_correct_spacing():
    base, drug = 148000.0, VCMMAE_NOMINAL_DA
    env = simulate_dar_envelope(base, drug, mean_dar=4.0, max_sites=8)
    assert len(env) == 9  # DAR 0..8
    assert abs(sum(a for _, a in env) - 1.0) < 1e-9
    # Spacing between adjacent DAR species is exactly one drug-linker mass.
    for i in range(len(env) - 1):
        assert abs((env[i + 1][0] - env[i][0]) - drug) < 1e-6
    # p = 4/8 = 0.5 -> symmetric binomial, C(8,k)/256.
    expected = [math.comb(8, k) / 256.0 for k in range(9)]
    for (_, actual), want in zip(env, expected):
        assert abs(actual - want) < 1e-9
    # Most abundant species sits at the mean DAR.
    top_mass = max(env, key=lambda ma: ma[1])[0]
    assert abs(top_mass - (base + 4 * drug)) < 1e-6


def test_adc_report_uses_dar_mode_and_flags_cooccurring_glycans():
    rep = heterogeneity_report(
        seq(IGG1_HC),
        protein_class="adc",
        drug_mass_da=VCMMAE_NOMINAL_DA,
        mean_dar=4.0,
    )
    assert rep["applies"] is True
    assert rep["mode"] == "dar"
    assert rep["computed"] is True
    assert rep["n_proteoforms"] == 9
    # IgG1 HC carries one sequon; the report must say the DAR ladder alone
    # understates the true envelope.
    assert any("sequon" in note for note in rep.get("assumptions", []))


def test_adc_without_drug_mass_declines_rather_than_guessing():
    rep = heterogeneity_report(seq(IGG1_HC), protein_class="adc")
    assert rep["applies"] is False
    assert "drug_mass_da" in rep["reason"]


def test_dar_rejects_impossible_mean():
    for bad in (-1.0, 9.0):
        try:
            simulate_dar_envelope(148000.0, 1316.6, mean_dar=bad, max_sites=8)
        except ValueError:
            continue
        raise AssertionError(f"mean_dar={bad} should have raised")


# ---------------------------------------------------------------------------
# Case 3: a non-glycosylated protein -> applies: false
# ---------------------------------------------------------------------------

def test_non_glycosylated_protein_does_not_apply():
    for uid in (UBIQUITIN, CARBONIC_ANHYDRASE):
        s = seq(uid)
        assert n_glyc_sites(s) == 0, uid
        rep = heterogeneity_report(s)
        assert rep["applies"] is False, uid
        assert rep["n_sites"] == 0
        assert rep["predicted_envelope_width_da"] == 0.0
        assert rep["resolved_at_resolution"] is True
        assert rep["level"] == "none"
        assert rep["computed"] is True


def test_too_short_sequence_returns_not_applicable_without_raising():
    rep = heterogeneity_report("MA")
    assert rep["applies"] is False
    assert rep["computed"] is True


# ---------------------------------------------------------------------------
# Resolvability arithmetic
# ---------------------------------------------------------------------------

def test_resolvability_single_species_is_trivially_resolved():
    out = resolvability([(50000.0, 1.0)])
    assert out["resolved"] is True
    assert out["fraction_resolved"] == 1.0
    assert out["worst_overlap_da"] == 0.0


def test_resolvability_matches_hand_arithmetic():
    """At M = 30,000 and R = 30,000 the FWHM is ~1 Da, so a 0.5 Da gap fails."""
    resolved = resolvability([(30000.0, 0.5), (30002.0, 0.5)], resolution=30000)
    assert resolved["resolved"] is True
    unresolved = resolvability([(30000.0, 0.5), (30000.5, 0.5)], resolution=30000)
    assert unresolved["resolved"] is False
    assert unresolved["worst_overlap_da"] > 0.0
    # FWHM = (M + z*proton)/R, just above 1 Da at this mass.
    assert 1.0 < unresolved["required_spacing_da"] < 1.1


def test_higher_resolution_never_resolves_fewer_pairs():
    env = simulate_glycoform_envelope(150000.0, 6)
    fractions = [
        resolvability(env, resolution=r)["fraction_resolved"]
        for r in (1000, 10000, 30000, 100000, 1000000)
    ]
    assert fractions == sorted(fractions)


def test_minimum_resolution_required_is_consistent_with_verdict():
    env = simulate_glycoform_envelope(150000.0, 4)
    need = minimum_resolution_required(env)
    assert need is not None and need > 0
    assert resolvability(env, resolution=need * 1.01)["resolved"] is True
    assert resolvability(env, resolution=need * 0.99)["resolved"] is False


def test_verdict_is_charge_independent():
    """The mass-domain criterion must not hinge on a guessed charge state."""
    env = simulate_glycoform_envelope(150000.0, 6)
    verdicts = {
        resolvability(env, resolution=30000, mz_charge=z)["resolved"]
        for z in (10, 25, 40, 60)
    }
    assert len(verdicts) == 1


# ---------------------------------------------------------------------------
# Contract required by METHODS_1_2_BUILD_PLAN.md
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Sequon occupancy (macroheterogeneity)
# ---------------------------------------------------------------------------

def test_full_occupancy_produces_no_aglycosylated_species():
    """occupancy=1.0 is the fully-occupied limit: nothing sits at base mass."""
    env = simulate_glycoform_envelope(100000.0, 3, occupancy=1.0)
    assert all(m > 100000.0 for m, _ in env)


def test_partial_occupancy_adds_the_aglycosylated_form():
    """With occupancy < 1 the fully-unglycosylated species must exist."""
    env = simulate_glycoform_envelope(100000.0, 3, occupancy=0.75)
    masses = [m for m, _ in env]
    assert min(masses) == 100000.0
    # P(all three sites empty) = 0.25^3
    abundance = dict(env)[100000.0]
    assert abs(abundance - 0.25 ** 3) < 1e-6


def test_occupancy_default_is_75_percent():
    assert abs(DEFAULT_SEQUON_OCCUPANCY - 0.75) < 1e-9


def test_zero_occupancy_collapses_to_a_single_species():
    assert simulate_glycoform_envelope(90000.0, 5, occupancy=0.0) == [(90000.0, 1.0)]


def test_occupancy_out_of_range_raises():
    for bad in (-0.1, 1.1):
        try:
            simulate_glycoform_envelope(90000.0, 2, occupancy=bad)
        except ValueError:
            continue
        raise AssertionError(f"occupancy={bad} should have raised")


def test_partial_occupancy_widens_the_simulated_span():
    """Occupancy extends the envelope down to base mass, it does not narrow it."""
    full = simulate_glycoform_envelope(100000.0, 4, occupancy=1.0)
    part = simulate_glycoform_envelope(100000.0, 4, occupancy=0.75)
    span_full = full[-1][0] - full[0][0]
    span_part = part[-1][0] - part[0][0]
    assert span_part > span_full


# ---------------------------------------------------------------------------
# Detection threshold
# ---------------------------------------------------------------------------

def test_detection_threshold_default_is_one_percent():
    assert abs(DEFAULT_DETECTION_THRESHOLD - 0.01) < 1e-9


def test_detectable_species_filters_on_base_peak():
    env = [(1000.0, 1.0), (1001.0, 0.02), (1002.0, 0.001)]
    keep = detectable_species(env, 0.01)
    assert [m for m, _ in keep] == [1000.0, 1001.0]


def test_faint_species_do_not_drive_the_verdict():
    """A 0.1% species jammed against the base peak must not call it unresolved."""
    env = [(30000.0, 1.0), (30000.2, 0.001), (30010.0, 0.5)]
    strict = resolvability(env, resolution=30000, min_detectable_abundance=0.0)
    lenient = resolvability(env, resolution=30000, min_detectable_abundance=0.01)
    assert strict["resolved"] is False       # 0.2 Da gap vs ~1 Da FWHM
    assert lenient["resolved"] is True       # the faint species is undetectable
    assert lenient["n_species"] == 2
    assert lenient["n_species_simulated"] == 3


def test_threshold_zero_recovers_strict_all_species_behaviour():
    env = simulate_glycoform_envelope(150000.0, 6)
    strict = resolvability(env, resolution=30000, min_detectable_abundance=0.0)
    assert strict["n_species"] == strict["n_species_simulated"]


def test_min_resolution_required_honours_the_same_threshold():
    env = simulate_glycoform_envelope(150000.0, 5)
    for thr in (0.0, 0.01, 0.05):
        need = minimum_resolution_required(env, min_detectable_abundance=thr)
        if need is None:
            continue
        assert resolvability(
            env, resolution=need * 1.01, min_detectable_abundance=thr
        )["resolved"] is True
        assert resolvability(
            env, resolution=need * 0.99, min_detectable_abundance=thr
        )["resolved"] is False


# ---------------------------------------------------------------------------
# Mandatory assumption disclosure
# ---------------------------------------------------------------------------

def test_every_report_states_its_assumptions():
    """Including reports that conclude heterogeneity does not apply."""
    cases = [
        seq(THYROGLOBULIN),      # applies, unresolvable
        seq(TRANSFERRIN),        # applies, resolvable
        seq(UBIQUITIN),          # does not apply
        seq(CARBONIC_ANHYDRASE), # does not apply
        "MA",                    # too short to featurise
    ]
    for s in cases:
        rep = heterogeneity_report(s)
        notes = " ".join(rep["assumptions"]).lower()
        assert "o-glycosylation is not modelled" in notes
        assert "occupancy assumed at 75%" in notes
        assert "glycan profile assumed" in notes
        assert "base peak" in notes
        detail = rep["model_assumptions"]
        assert detail["sequon_occupancy"] == DEFAULT_SEQUON_OCCUPANCY
        assert detail["o_glycosylation_modelled"] is False
        assert detail["detection_threshold_rel_abundance"] == DEFAULT_DETECTION_THRESHOLD
        assert detail["glycan_profile"] == DEFAULT_N_GLYCAN_DIST


def test_report_records_a_caller_supplied_profile_as_such():
    custom = {"Man5": 0.5, "Man9": 0.5}
    rep = heterogeneity_report(seq(TRANSFERRIN), glycans_per_site_dist=custom)
    assert rep["model_assumptions"]["glycan_profile_source"] == "caller-supplied"
    assert rep["model_assumptions"]["glycan_profile"] == custom
    assert "caller-supplied" in " ".join(rep["assumptions"])


def test_predict_annotation_is_strictly_additive():
    """/predict must gain a field and change nothing else.

    The model is stubbed so this exercises the integration without loading
    ESM-2. Skipped if fastapi's test client is unavailable.
    """
    try:
        from fastapi.testclient import TestClient
    except Exception:  # pragma: no cover - optional dependency
        print("      (skipped: fastapi TestClient unavailable)")
        return

    import copy
    import tempfile
    import main

    baseline = {
        "suitability_score": 87,
        "suitability_label": "Excellent",
        "confidence_interval": {"lower": 78, "upper": 96},
        "risk_factors": [],
        "recommendations": [],
        "model_version": "0.4-esm2-glyco-tm",
    }
    saved_predict, saved_log = main.run_prediction, main.PREDICTION_LOG
    main.run_prediction = lambda s: copy.deepcopy(baseline)
    main.PREDICTION_LOG = Path(tempfile.gettempdir()) / "nr_test_predictions.jsonl"
    try:
        client = TestClient(main.app)

        # Non-glycosylated: response must be byte-identical to the model output.
        plain = client.post("/predict", json={"sequence": seq(CARBONIC_ANHYDRASE)})
        assert plain.status_code == 200
        assert plain.json() == baseline
        assert "heterogeneity_risk" not in plain.json()

        # Glycoprotein: annotated, but every pre-existing field untouched.
        glyco = client.post("/predict", json={"sequence": seq(THYROGLOBULIN)})
        assert glyco.status_code == 200
        body = glyco.json()
        assert "heterogeneity_risk" in body
        for key, value in baseline.items():
            assert body[key] == value, f"existing field '{key}' was modified"
        assert body["heterogeneity_risk"]["computed"] is True
        assert body["heterogeneity_risk"]["level"] == "high"

        # A simulator failure must not cost the caller their prediction.
        saved_report = main.heterogeneity_report
        main.heterogeneity_report = lambda s: 1 / 0
        try:
            broken = client.post("/predict", json={"sequence": seq(THYROGLOBULIN)})
            assert broken.status_code == 200
            assert broken.json() == baseline
        finally:
            main.heterogeneity_report = saved_report
    finally:
        main.run_prediction, main.PREDICTION_LOG = saved_predict, saved_log


def test_report_always_carries_the_specified_keys():
    required = {
        "applies",
        "n_sites",
        "predicted_envelope_width_da",
        "resolved_at_resolution",
        "reason",
    }
    for s in (seq(THYROGLOBULIN), seq(TRANSFERRIN), seq(UBIQUITIN), "MA"):
        rep = heterogeneity_report(s)
        assert required.issubset(rep.keys())
        assert rep["computed"] is True


# ---------------------------------------------------------------------------
# Standalone runner (no pytest required)
# ---------------------------------------------------------------------------

def _main():
    tests = [
        (name, fn)
        for name, fn in sorted(globals().items())
        if name.startswith("test_") and callable(fn)
    ]
    passed, failed = 0, []
    print(f"Running {len(tests)} tests from {Path(__file__).name}\n")
    for name, fn in tests:
        try:
            fn()
        except Exception as exc:  # noqa: BLE001
            failed.append((name, exc))
            print(f"  FAIL  {name}\n          {type(exc).__name__}: {exc}")
        else:
            passed += 1
            print(f"  ok    {name}")
    print(f"\n{passed} passed, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(_main())
