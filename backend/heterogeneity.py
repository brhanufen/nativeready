"""Heterogeneity forward-simulator (Method 1).

Given a sequence, SIMULATE the proteoform mass envelope that glycosylation (or
ADC drug conjugation) will produce, then COMPUTE whether those proteoforms can
be resolved at a given instrument resolution.

This module contains no machine learning and no fitted parameters. Every number
it returns is arithmetic over published monosaccharide masses plus the standard
definition of mass-spectrometric resolving power. A reviewer can redo any of it
by hand. It never sees a failure label and it never trains on anything.

Physics / mass conventions
--------------------------
Monoisotopic *residue* masses (the mass added to a glycan chain, i.e. the free
monosaccharide minus one H2O). Derived from CIAAW monoisotopic atomic masses
(C 12.000000, H 1.00782503, N 14.00307400, O 15.99491462) and identical to the
standard glycomics tables following Domon & Costello, Glycoconj. J. 5:397-409
(1988):

    Hex    (Man/Gal/Glc)  C6H10O5   162.0528
    HexNAc (GlcNAc/GalNAc) C8H13NO5 203.0794
    dHex   (Fuc)          C6H10O4   146.0579
    NeuAc  (Neu5Ac)       C11H17NO8 291.0954
    NeuGc  (Neu5Gc)       C11H17NO9 307.0903
    Pent   (Xyl)          C5H8O4    132.0423

Resolving power is defined as R = m / FWHM in the m/z domain. Two adjacent
proteoforms are called resolved when their spacing is at least one FWHM (the
conventional ~50%-valley criterion). Converting to the neutral-mass domain, the
required spacing is (M + z*1.00728) / R, which to first order is simply M / R
and is therefore independent of charge state -- so the resolvability verdict
does not depend on guessing a charge. Charge only affects the reported m/z-domain
view, and where a charge is needed and not supplied we use the charge-residue
model z ~ 0.0778*sqrt(M) (Fernandez de la Mora, Anal. Chim. Acta 406:93 (2000)),
clearly labelled as an estimate.

Modelling choices
-----------------
These are ASSUMPTIONS, not measurements. Each is exposed as an overridable
argument, each is stated in every report this module emits, and each was chosen
on physical grounds without reference to the v7 positive/negative labels:

  1. DEFAULT_N_GLYCAN_DIST -- a nominal per-site glycan profile typical of a
     CHO-expressed recombinant IgG. Real profiles are cell-line, clone and
     process dependent. Pass `glycans_per_site_dist` when you know the profile.

  2. DEFAULT_SEQUON_OCCUPANCY = 0.75. An N-X-S/T sequon is a site that *can* be
     glycosylated, not one that always is. Occupancy is incomplete because the
     oligosaccharyltransferase competes with co-translational folding for access
     to the nascent chain, so sequons that fold quickly or sit in constrained
     structural context are skipped -- the phenomenon conventionally called
     macroheterogeneity. Occupancy is also sequence dependent (N-X-T is
     transferred more efficiently than N-X-S). Treating every sequon as fully
     occupied therefore overstates how many distinct glycan-bearing species
     exist. 0.75 is a nominal working figure standing in for a broad and
     protein-specific distribution; it is a modelling choice, not a constant of
     nature. Set occupancy=1.0 to recover the fully-occupied limit.

  3. DEFAULT_DETECTION_THRESHOLD = 0.01. Resolvability is judged only over
     species at or above 1% of the base-peak abundance. Demanding that an
     instrument separate proteoforms present at 0.01% is not a physically
     meaningful test: such species sit at or below the noise floor of a native
     spectrum of a heterogeneous glycoprotein and are not observed at all,
     resolved or otherwise. The simulator's internal pruning floor
     (min_rel_abundance, 1e-4) is a separate and much lower *computational*
     cutoff -- it bounds the convolution, it does not define what is detectable.

  4. ADC drug-linker masses are NOT defaulted. `simulate_dar_envelope` requires
     an explicit drug_mass_da, because published linker-payload masses vary by
     construct and quoting one to 0.1 Da would be false precision.

Known scope limit: only N-linked glycosylation is modelled. O-glycosylation
contributes no sequon motif to count and is not simulated, so mucin-type and
hinge-region O-glycan heterogeneity is invisible here. Every report says so
explicitly, including reports that conclude heterogeneity does not apply --
that is precisely the case where an O-glycosylated protein would be wrongly
cleared.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

# Reuse the sequon logic that the shipped v0.4 model already uses. Do NOT
# reimplement it here -- the whole point is that the simulator counts sites the
# same way the feature vector does.
from features import extract_features

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Monoisotopic residue masses in daltons. See module docstring for derivation.
GLYCAN_MASSES: Dict[str, float] = {
    "Hex": 162.0528,
    "HexNAc": 203.0794,
    "Fuc": 146.0579,      # dHex
    "NeuAc": 291.0954,
    "NeuGc": 307.0903,
    "Pent": 132.0423,     # xylose, plant/insect N-glycans
}

#: Mass of water, for converting between free-glycan and glycan-residue masses.
WATER_MASS = 18.0106

#: Proton mass, for neutral-mass <-> m/z conversion.
PROTON_MASS = 1.007276

#: Charge-residue-model coefficient, <z> ~ CRM_COEFF * sqrt(M).
#: Fernandez de la Mora, Anal. Chim. Acta 406:93 (2000).
CRM_COEFF = 0.0778

#: Named N-glycan compositions, as counts of the residues in GLYCAN_MASSES.
N_GLYCAN_COMPOSITIONS: Dict[str, Dict[str, int]] = {
    "Man5":  {"HexNAc": 2, "Hex": 5},
    "Man6":  {"HexNAc": 2, "Hex": 6},
    "Man7":  {"HexNAc": 2, "Hex": 7},
    "Man8":  {"HexNAc": 2, "Hex": 8},
    "Man9":  {"HexNAc": 2, "Hex": 9},
    "G0":    {"HexNAc": 4, "Hex": 3},
    "G0F":   {"HexNAc": 4, "Hex": 3, "Fuc": 1},
    "G1":    {"HexNAc": 4, "Hex": 4},
    "G1F":   {"HexNAc": 4, "Hex": 4, "Fuc": 1},
    "G2":    {"HexNAc": 4, "Hex": 5},
    "G2F":   {"HexNAc": 4, "Hex": 5, "Fuc": 1},
    "G2FS1": {"HexNAc": 4, "Hex": 5, "Fuc": 1, "NeuAc": 1},
    "G2FS2": {"HexNAc": 4, "Hex": 5, "Fuc": 1, "NeuAc": 2},
}

#: ASSUMPTION (not a measurement): nominal per-site N-glycan occupancy for a
#: CHO-expressed recombinant IgG, normalised to 1.0. Override via the
#: `glycans_per_site_dist` argument whenever the real profile is known.
DEFAULT_N_GLYCAN_DIST: Dict[str, float] = {
    "G0F":   0.45,
    "G1F":   0.30,
    "G2F":   0.10,
    "Man5":  0.06,
    "G0":    0.05,
    "G2FS1": 0.04,
}

#: ASSUMPTION (see module docstring, item 2): fraction of N-X-S/T sequons that
#: actually carry a glycan. Incomplete occupancy is the well-established
#: macroheterogeneity phenomenon; 1.0 is the fully-occupied limit.
DEFAULT_SEQUON_OCCUPANCY = 0.75

#: ASSUMPTION (see module docstring, item 3): species below this fraction of the
#: base-peak abundance are treated as undetectable and excluded from the
#: resolvability verdict. Distinct from the simulator's computational pruning
#: floor, which is far lower.
DEFAULT_DETECTION_THRESHOLD = 0.01

#: Conventional peak-separation criterion: spacing >= SEPARATION_FACTOR * FWHM.
SEPARATION_FACTOR = 1.0

#: Protein classes in dataset_combined_v7 that flag an antibody-drug conjugate.
_ADC_CLASSES = {"adc", "antibody-drug conjugate"}


# ---------------------------------------------------------------------------
# Sites
# ---------------------------------------------------------------------------

def glycan_mass(composition: Dict[str, int]) -> float:
    """Monoisotopic mass of a glycan given its residue composition.

    >>> round(glycan_mass(N_GLYCAN_COMPOSITIONS["G0F"]), 4)
    1444.5339
    """
    return sum(GLYCAN_MASSES[res] * n for res, n in composition.items())


def n_glyc_sites(sequence: str, strict: bool = False) -> int:
    """Number of N-linked glycosylation sequons (N-X-S/T, X != P).

    Delegates to ``features.extract_features`` so the simulator counts sites
    exactly the way the shipped model's feature vector does. ``strict=True``
    additionally excludes X = C, matching ``n_glyc_sequon_strict_count``.

    Returns 0 for sequences too short to featurise, rather than raising.
    """
    try:
        feats = extract_features(sequence)
    except (ValueError, KeyError):
        return 0
    key = "n_glyc_sequon_strict_count" if strict else "n_glyc_sequon_count"
    return int(feats.get(key, 0.0))


# ---------------------------------------------------------------------------
# Envelope simulation
# ---------------------------------------------------------------------------

def _normalise(dist: Dict[str, float]) -> Dict[str, float]:
    total = float(sum(dist.values()))
    if total <= 0:
        raise ValueError("glycan distribution must have positive total weight")
    return {k: v / total for k, v in dist.items()}


def _prune(
    species: Dict[float, float],
    min_rel_abundance: float,
    max_species: int,
) -> Tuple[Dict[float, float], bool]:
    """Drop negligible species and cap the count. Returns (species, truncated)."""
    truncated = False
    if not species:
        return species, truncated
    peak = max(species.values())
    kept = {m: a for m, a in species.items() if a >= min_rel_abundance * peak}
    if len(kept) != len(species):
        truncated = True
    if len(kept) > max_species:
        top = sorted(kept.items(), key=lambda kv: kv[1], reverse=True)[:max_species]
        kept = dict(top)
        truncated = True
    return kept, truncated


def simulate_glycoform_envelope(
    base_mass_da: float,
    n_sites: int,
    glycans_per_site_dist: Optional[Dict[str, float]] = None,
    occupancy: float = DEFAULT_SEQUON_OCCUPANCY,
    min_rel_abundance: float = 1e-4,
    max_species: int = 2000,
    max_sites_simulated: int = 24,
) -> List[Tuple[float, float]]:
    """Convolve a per-site glycan distribution across ``n_sites``.

    Each site is independently occupied with probability ``occupancy`` and, when
    occupied, carries one glycan drawn from ``glycans_per_site_dist`` (names must
    be keys of ``N_GLYCAN_COMPOSITIONS``). An unoccupied site contributes 0 Da.
    The resulting per-site distribution is convolved ``n_sites`` times;
    combinations that are exactly isobaric are merged, which is required rather
    than optional -- G0F+G2F and G1F+G1F have identical mass by construction.

    Modelling occupancy < 1 adds the partially-glycosylated species that a real
    glycoprotein carries. It does not simply narrow the envelope: it extends it
    down towards ``base_mass_da`` (the fully-aglycosylated form) while
    redistributing abundance, so the observable consequence has to be measured
    rather than assumed.

    Returns ``[(mass_da, relative_abundance), ...]`` sorted by mass, with
    abundances normalised to sum to 1.0. An unglycosylated protein (n_sites=0)
    or fully-unoccupied protein (occupancy=0) returns a single species at
    ``base_mass_da``.

    ``max_sites_simulated`` bounds the convolution for extreme cases such as
    mucins with dozens of sequons; the envelope is already unresolvable long
    before the cap, and ``heterogeneity_report`` records when it was hit.
    """
    if not 0.0 <= occupancy <= 1.0:
        raise ValueError("occupancy must lie in [0, 1]")
    if n_sites <= 0 or occupancy == 0.0:
        return [(float(base_mass_da), 1.0)]

    dist = _normalise(glycans_per_site_dist or DEFAULT_N_GLYCAN_DIST)

    # Resolve names to masses once.
    site_masses: List[Tuple[float, float]] = []
    for name, weight in dist.items():
        if isinstance(name, str):
            comp = N_GLYCAN_COMPOSITIONS.get(name)
            if comp is None:
                raise KeyError(
                    f"unknown glycan '{name}'; add it to N_GLYCAN_COMPOSITIONS "
                    "or pass explicit masses"
                )
            site_masses.append((glycan_mass(comp), weight))
        else:  # pragma: no cover - defensive
            raise TypeError("glycan distribution keys must be glycan names")

    # Bernoulli occupancy: an unoccupied site adds no mass. Weights already sum
    # to 1 across (unoccupied + all glycans).
    site_masses = [(m, w * occupancy) for m, w in site_masses]
    if occupancy < 1.0:
        site_masses.append((0.0, 1.0 - occupancy))

    sites = min(int(n_sites), max_sites_simulated)

    # Iterative convolution. Keys are masses rounded to 4 dp so that genuinely
    # isobaric combinations merge and near-isobars stay distinct; the
    # resolvability step, not the binning, decides what an instrument can see.
    species: Dict[float, float] = {0.0: 1.0}
    for _ in range(sites):
        nxt: Dict[float, float] = {}
        for acc_mass, acc_ab in species.items():
            for gmass, gweight in site_masses:
                key = round(acc_mass + gmass, 4)
                nxt[key] = nxt.get(key, 0.0) + acc_ab * gweight
        species, _ = _prune(nxt, min_rel_abundance, max_species)

    total = sum(species.values())
    envelope = sorted(
        (round(float(base_mass_da) + m, 4), a / total) for m, a in species.items()
    )
    return envelope


def simulate_dar_envelope(
    base_mass_da: float,
    drug_mass_da: float,
    mean_dar: float = 4.0,
    max_sites: int = 8,
    min_rel_abundance: float = 1e-4,
) -> List[Tuple[float, float]]:
    """Binomial drug-to-antibody-ratio ladder for an ADC.

    Models ``max_sites`` independent conjugation sites each occupied with
    probability ``p = mean_dar / max_sites``; the resulting DAR distribution is
    Binomial(max_sites, p) and species k has mass
    ``base_mass_da + k * drug_mass_da``. The default of 8 sites corresponds to
    cysteine conjugation across the four interchain disulfides of an IgG1.

    ``drug_mass_da`` is deliberately required: published linker-payload masses
    are construct-specific and defaulting one would be false precision.
    """
    if max_sites <= 0:
        raise ValueError("max_sites must be positive")
    if not 0.0 <= mean_dar <= max_sites:
        raise ValueError(f"mean_dar must lie in [0, {max_sites}]")

    p = mean_dar / max_sites
    species: List[Tuple[float, float]] = []
    for k in range(max_sites + 1):
        prob = math.comb(max_sites, k) * (p ** k) * ((1.0 - p) ** (max_sites - k))
        species.append((round(float(base_mass_da) + k * float(drug_mass_da), 4), prob))

    peak = max(a for _, a in species) or 1.0
    kept = [(m, a) for m, a in species if a >= min_rel_abundance * peak]
    total = sum(a for _, a in kept)
    return [(m, a / total) for m, a in kept]


# ---------------------------------------------------------------------------
# Resolvability
# ---------------------------------------------------------------------------

def estimate_native_charge(mass_da: float) -> float:
    """Charge-residue-model estimate of the native charge state, 0.0778*sqrt(M).

    An analytic estimate for a compact globular protein, not a measurement. Used
    only to render the m/z-domain view; the resolvability verdict is computed in
    the neutral-mass domain and does not depend on it.
    """
    return CRM_COEFF * math.sqrt(max(float(mass_da), 1.0))


def detectable_species(
    envelope: Sequence[Tuple[float, float]],
    min_detectable_abundance: float = DEFAULT_DETECTION_THRESHOLD,
) -> List[Tuple[float, float]]:
    """Species at or above ``min_detectable_abundance`` x the base-peak abundance.

    Relative abundance is referenced to the most abundant species, the usual MS
    convention. Species below the threshold are not observed in a native
    spectrum at all, so asking whether they are *resolved* is not a meaningful
    physical question.
    """
    pts = sorted((float(m), float(a)) for m, a in envelope)
    if not pts:
        return []
    peak = max(a for _, a in pts)
    if peak <= 0:
        return pts
    return [(m, a) for m, a in pts if a >= min_detectable_abundance * peak]


def resolvability(
    envelope: Sequence[Tuple[float, float]],
    resolution: float = 30000.0,
    mz_charge: Optional[float] = None,
    separation_factor: float = SEPARATION_FACTOR,
    min_detectable_abundance: float = DEFAULT_DETECTION_THRESHOLD,
) -> Dict[str, Any]:
    """Can the *detectable* species in ``envelope`` be told apart at ``resolution``?

    Resolving power R = m / FWHM. Two adjacent proteoforms are resolved when
    their spacing is at least ``separation_factor`` FWHM. Working in the
    neutral-mass domain the required spacing is (M + z*proton)/R, which reduces
    to M/R and is independent of z to first order -- so the verdict holds
    whatever charge state is observed.

    The verdict is computed over species at or above
    ``min_detectable_abundance`` x the base peak. Set it to 0 to judge every
    simulated species, including ones below any realistic noise floor.

    Returns ``{resolved, fraction_resolved, worst_overlap_da, ...}``.
    """
    all_pts = sorted((float(m), float(a)) for m, a in envelope)
    n_total = len(all_pts)
    pts = detectable_species(all_pts, min_detectable_abundance)
    n = len(pts)
    if n <= 1:
        return {
            "resolved": True,
            "fraction_resolved": 1.0,
            "worst_overlap_da": 0.0,
            "n_species": n,
            "n_species_simulated": n_total,
            "detection_threshold": float(min_detectable_abundance),
            "envelope_width_da": 0.0,
            "simulated_width_da": round(
                (all_pts[-1][0] - all_pts[0][0]) if n_total > 1 else 0.0, 4
            ),
            "resolution": float(resolution),
            "criterion": f"spacing >= {separation_factor:g} x FWHM (FWHM = M/R)",
            "n_adjacent_pairs": 0,
        }

    masses = [m for m, _ in pts]
    width = masses[-1] - masses[0]
    simulated_width = all_pts[-1][0] - all_pts[0][0] if n_total > 1 else 0.0
    charge = float(mz_charge) if mz_charge else estimate_native_charge(
        sum(masses) / n
    )

    resolved_pairs = 0
    worst_overlap = 0.0
    tightest = None
    for i in range(n - 1):
        spacing = masses[i + 1] - masses[i]
        local_mass = 0.5 * (masses[i] + masses[i + 1])
        # FWHM in the neutral-mass domain at this mass and charge.
        required = separation_factor * (local_mass + charge * PROTON_MASS) / float(
            resolution
        )
        if spacing >= required:
            resolved_pairs += 1
        else:
            overlap = required - spacing
            if overlap > worst_overlap:
                worst_overlap = overlap
                tightest = (masses[i], masses[i + 1], spacing, required)

    pairs = n - 1
    frac = resolved_pairs / pairs
    out: Dict[str, Any] = {
        "resolved": resolved_pairs == pairs,
        "fraction_resolved": round(frac, 4),
        "worst_overlap_da": round(worst_overlap, 4),
        "n_species": n,
        "n_species_simulated": n_total,
        "detection_threshold": float(min_detectable_abundance),
        "envelope_width_da": round(width, 4),
        "simulated_width_da": round(simulated_width, 4),
        "resolution": float(resolution),
        "criterion": f"spacing >= {separation_factor:g} x FWHM (FWHM = M/R)",
        "n_adjacent_pairs": pairs,
        "charge_used": round(charge, 2),
        "charge_is_estimate": mz_charge is None,
    }
    if tightest is not None:
        lo, hi, spacing, required = tightest
        out["tightest_pair_da"] = [round(lo, 4), round(hi, 4)]
        out["tightest_spacing_da"] = round(spacing, 4)
        out["required_spacing_da"] = round(required, 4)
    # m/z-domain view (informational).
    if charge > 0:
        mean_mass = sum(masses) / n
        out["mz_at_charge"] = round((mean_mass + charge * PROTON_MASS) / charge, 3)
        out["fwhm_mz"] = round(out["mz_at_charge"] / float(resolution), 5)
    return out


def minimum_resolution_required(
    envelope: Sequence[Tuple[float, float]],
    mz_charge: Optional[float] = None,
    separation_factor: float = SEPARATION_FACTOR,
    min_detectable_abundance: float = DEFAULT_DETECTION_THRESHOLD,
) -> Optional[float]:
    """Smallest resolving power at which every *detectable* adjacent pair separates.

    Uses the same detection threshold as ``resolvability`` so the two always
    agree. Returns None when fewer than two species are detectable, and math.inf
    if two detectable species are exactly isobaric (no resolution separates
    them).
    """
    pts = sorted(
        float(m)
        for m, _ in detectable_species(envelope, min_detectable_abundance)
    )
    if len(pts) <= 1:
        return None
    charge = float(mz_charge) if mz_charge else estimate_native_charge(
        sum(pts) / len(pts)
    )
    worst = 0.0
    for i in range(len(pts) - 1):
        spacing = pts[i + 1] - pts[i]
        local_mass = 0.5 * (pts[i] + pts[i + 1])
        if spacing <= 0:
            return math.inf
        needed = separation_factor * (local_mass + charge * PROTON_MASS) / spacing
        worst = max(worst, needed)
    return worst


# ---------------------------------------------------------------------------
# Top-level report
# ---------------------------------------------------------------------------

def _level(applies: bool, resolved: bool, width: float, n_species: int) -> str:
    if not applies:
        return "none"
    if resolved:
        return "low"
    if width >= 1000.0 or n_species >= 20:
        return "high"
    return "medium"


def _assumption_block(
    glycans_per_site_dist: Optional[Dict[str, float]],
    occupancy: float,
    min_detectable_abundance: float,
) -> Tuple[List[str], Dict[str, Any]]:
    """Build the assumption disclosure carried by every report.

    Emitted on every return path, including ``applies: False``. A protein with
    no N-sequons is exactly the case where an unmodelled O-glycan load would be
    silently cleared, so the caller has to be told what was not simulated even
    when the answer is "no risk".
    """
    supplied = glycans_per_site_dist is not None
    dist = glycans_per_site_dist or DEFAULT_N_GLYCAN_DIST
    profile_txt = ", ".join(
        f"{name} {100 * frac:.0f}%"
        for name, frac in sorted(dist.items(), key=lambda kv: -kv[1])
    )
    source = "caller-supplied" if supplied else "nominal CHO-IgG default"

    notes = [
        f"Glycan profile assumed ({source}): {profile_txt}. "
        "Real profiles are cell-line, clone and process dependent.",
        f"Sequon occupancy assumed at {100 * occupancy:.0f}%: sequons are sites "
        "that can be glycosylated, not sites that always are (macroheterogeneity). "
        "Occupancy is protein-specific and this is a nominal figure.",
        "O-glycosylation is NOT modelled. O-glycans carry no sequon motif to "
        "count, so mucin-type and hinge-region O-glycan heterogeneity is invisible "
        "to this calculation and any such load is additional to what is reported.",
        f"Resolvability judged over species at or above "
        f"{100 * min_detectable_abundance:g}% of the base peak; fainter species "
        "are below the practical detection floor of a native spectrum.",
    ]
    detail = {
        "sequon_occupancy": occupancy,
        "glycan_profile_source": source,
        "glycan_profile": dict(dist),
        "detection_threshold_rel_abundance": min_detectable_abundance,
        "o_glycosylation_modelled": False,
    }
    return notes, detail


def heterogeneity_report(
    sequence: str,
    protein_class: Optional[str] = None,
    resolution: float = 30000.0,
    glycans_per_site_dist: Optional[Dict[str, float]] = None,
    drug_mass_da: Optional[float] = None,
    mean_dar: float = 4.0,
    mz_charge: Optional[float] = None,
    base_mass_da: Optional[float] = None,
    occupancy: float = DEFAULT_SEQUON_OCCUPANCY,
    min_detectable_abundance: float = DEFAULT_DETECTION_THRESHOLD,
) -> Dict[str, Any]:
    """Top-level entry point: does proteoform heterogeneity threaten this target?

    Returns a dict with at least ``applies``, ``n_sites``,
    ``predicted_envelope_width_da``, ``resolved_at_resolution``, ``reason``,
    ``assumptions`` and ``model_assumptions``. ``applies`` is False -- with the
    numeric fields zeroed but the assumptions still stated -- for a protein with
    no N-glycosylation sequons that is not flagged as an ADC.

    ``predicted_envelope_width_da`` is the span of the *detectable* species,
    i.e. what an instrument would actually show. The full simulated span is
    reported separately as ``simulated_width_da``.

    Never raises on ordinary bad input; a sequence too short to featurise simply
    returns ``applies: False``.
    """
    notes, detail = _assumption_block(
        glycans_per_site_dist, occupancy, min_detectable_abundance
    )

    def not_applicable(reason: str, n_sites: int = 0) -> Dict[str, Any]:
        return {
            "applies": False,
            "mode": None,
            "n_sites": n_sites,
            "predicted_envelope_width_da": 0.0,
            "resolved_at_resolution": True,
            "reason": reason,
            "level": "none",
            "assumptions": notes,
            "model_assumptions": detail,
            "computed": True,
        }

    cls = (protein_class or "").strip().lower()
    is_adc = cls in _ADC_CLASSES or drug_mass_da is not None

    # Base mass from the same BioPython featureisation the model uses.
    n_sites = 0
    if base_mass_da is None:
        try:
            feats = extract_features(sequence)
            base_mass_da = float(feats["molecular_weight_kda"]) * 1000.0
            n_sites = int(feats.get("n_glyc_sequon_count", 0.0))
        except (ValueError, KeyError):
            return not_applicable("Sequence too short to analyse for heterogeneity.")
    else:
        n_sites = n_glyc_sites(sequence)

    mode = None
    envelope: List[Tuple[float, float]] = [(float(base_mass_da), 1.0)]

    if is_adc:
        if drug_mass_da is None:
            return not_applicable(
                "Flagged as an ADC but no drug-linker mass was supplied; "
                "pass drug_mass_da to simulate the DAR ladder.",
                n_sites=n_sites,
            )
        mode = "dar"
        envelope = simulate_dar_envelope(
            base_mass_da, drug_mass_da, mean_dar=mean_dar
        )
        notes = notes + [
            f"Binomial DAR ladder, mean DAR {mean_dar:g} over 8 cysteine sites, "
            f"drug-linker {drug_mass_da:g} Da."
        ]
        if n_sites >= 1:
            notes = notes + [
                f"This construct also has {n_sites} N-glycosylation sequon(s). "
                "The DAR ladder is simulated on its own; the real proteoform "
                "envelope is the convolution of the DAR ladder with the glycan "
                "envelope and will be wider and denser than reported here."
            ]
    elif n_sites >= 1:
        mode = "glycan"
        envelope = simulate_glycoform_envelope(
            base_mass_da,
            n_sites,
            glycans_per_site_dist=glycans_per_site_dist,
            occupancy=occupancy,
        )
        if n_sites > 24:
            notes = notes + [
                f"{n_sites} sequons found; convolution capped at 24 sites "
                "(the envelope is already unresolvable well below the cap)."
            ]
    else:
        return not_applicable(
            "No N-glycosylation sequons and not an ADC; proteoform "
            "heterogeneity from N-glycans or conjugation does not apply."
        )

    res = resolvability(
        envelope,
        resolution=resolution,
        mz_charge=mz_charge,
        min_detectable_abundance=min_detectable_abundance,
    )
    width = res["envelope_width_da"]
    n_species = res["n_species"]
    min_r = minimum_resolution_required(
        envelope,
        mz_charge=mz_charge,
        min_detectable_abundance=min_detectable_abundance,
    )
    level = _level(True, res["resolved"], width, n_species)

    pct = 100 * min_detectable_abundance
    if res["resolved"]:
        reason = (
            f"~{n_species} detectable proteoforms (>={pct:g}% of base peak) span "
            f"{width:.0f} Da; resolvable at R = {resolution:,.0f}."
        )
    else:
        need = (
            f"above R ~ {min_r:,.0f}"
            if min_r not in (None, math.inf)
            else "at any resolution (isobaric species present)"
        )
        reason = (
            f"~{n_species} detectable proteoforms (>={pct:g}% of base peak) span "
            f"{width:.0f} Da, unresolvable below R = {resolution:,.0f} ({need})."
        )

    report: Dict[str, Any] = {
        "applies": True,
        "mode": mode,
        "n_sites": n_sites,
        "base_mass_da": round(float(base_mass_da), 2),
        "n_proteoforms": n_species,
        "n_proteoforms_simulated": res["n_species_simulated"],
        "predicted_envelope_width_da": width,
        "simulated_width_da": res["simulated_width_da"],
        "resolved_at_resolution": res["resolved"],
        "fraction_resolved": res["fraction_resolved"],
        "worst_overlap_da": res["worst_overlap_da"],
        "resolution": float(resolution),
        "min_resolution_required": (
            None if min_r is None else (None if min_r == math.inf else round(min_r))
        ),
        "level": level,
        "reason": reason,
        "assumptions": notes,
        "model_assumptions": detail,
        "computed": True,
    }
    return report
