# NativeReady Label Schema v0.1

**Status**: DRAFT, awaiting Brhanu review
**Date**: 2026-05-01
**Scope**: defines the per-measurement record format for the v0.3+ training dataset

This schema is the gatekeeper for the next 6 months of data work. Lock this before extracting more records, because every change forces re-extraction.

---

## 1. Storage format

- **One row per measurement, NOT per protein.** A protein attempted under three different buffer conditions = three rows.
- **File format**: Parquet for the table; JSON sidecar for the failure-attempts log on negative records.
- **Backup format**: JSON (current files are JSON; Parquet conversion is a v0.4 task).
- **Filename convention**: `dataset_v{version}_{YYYY-MM-DD}.parquet` (e.g., `dataset_v0.3_2026-05-15.parquet`).

## 2. Required columns

| Column | Type | Required | Description |
|---|---|---|---|
| `uniprot_id` | string | yes | UniProt accession (e.g., P00918) |
| `sequence_sha256` | string | yes | SHA-256 hash of sequence (allows synthetic constructs without UniProt entries) |
| `name` | string | yes | Protein full name with organism in parentheses |
| `protein_class` | string | yes | One of the controlled vocabulary terms below |
| `mw_kda` | float | yes | Theoretical molecular weight in kDa |
| `sequence` | string | yes | Full amino acid sequence (one-letter code, standard 20 + B,X,Z if present) |
| `sequence_length` | int | yes | Length of `sequence` |
| `oligomeric_state_expected` | string | optional | Expected oligomeric state (monomer, dimer, tetramer, hexamer, etc.) |
| `oligomeric_state_observed` | string | optional | Observed in native MS (may differ from expected) |
| `mass_theoretical_da` | float | optional | Theoretical mass of intact species (oligomer if applicable) |
| `mass_observed_da` | float | optional | Observed mass in native MS spectrum |
| `mass_error_ppm` | float | optional | Computed: `((mass_observed - mass_theoretical) / mass_theoretical) * 1e6` |
| `instrument_psims_id` | string | optional | HUPO PSI-MS CV term (e.g., `MS:1003029` for Orbitrap UHMR) |
| `buffer` | string | optional | Buffer used (e.g., "ammonium acetate") |
| `buffer_concentration_mM` | float | optional | Buffer concentration |
| `analyte_concentration_uM` | float | optional | Sample concentration |
| `activation_energy_V` | float | optional | Activation/HCD voltage |
| `is_idp` | bool | optional, default false | Flag if intrinsically disordered protein |
| `membrane_mimetic` | string | optional | Detergent class, nanodisc lipid, peptidisc, etc. (only for membrane proteins) |
| `outcome_label` | int | yes | See section 3 below |
| `ionization_score` | int (0-3) | optional | See section 3 |
| `complex_preservation_score` | int (0-3) | optional | See section 3 |
| `csd_quality_score` | int (0-3) | optional | See section 3 |
| `heterogeneity_score` | int (0-3) | optional | See section 3 |
| `notes` | string | yes | Free-text provenance: source, citation, instrument notes, anything operator-relevant |
| `evidence_url_or_doi` | string | yes for new records | URL or DOI of the source paper/dataset |
| `submitter_id` | string | optional | Lab/person identifier (use `traversa-pdb` for PDB-extracted, etc.) |
| `extraction_date` | string (ISO date) | yes | When this row was added to the dataset |

## 3. Outcome label and sub-scores

### Primary label: `outcome_label`

Five-level ordinal label. **This is what the v0.3 model trains on first.**

| Level | Name | Definition |
|---|---|---|
| 1 | `clean_native` | Resolved charge envelope, mass within 50 ppm of theoretical, no significant adducts beyond expected Na/K, complex stoichiometry intact if applicable. |
| 2 | `interpretable_native` | Charge states resolved, but adducts, modest heterogeneity, or partial sub-population dissociation present; mass assignable within 200 ppm. |
| 3 | `partial` | Some intact signal but dominated by sub-complex or dissociation products; or signal-to-noise so low that quantitation of variants is not possible. |
| 4 | `denatured_only` | No signal under native ESI from ammonium acetate; protein is detected only after acid/organic denaturation. |
| 5 | `failure` | No signal, or signal so contaminated by aggregation/heterogeneity that no useful mass can be extracted. |

**For backwards compatibility with the existing 232+260 dataset**: levels 1-3 collapse to `label=1` (positive), levels 4-5 collapse to `label=0` (negative). The legacy `label` column is still emitted alongside `outcome_label`.

### Sub-scores (0-3 each, optional during v0.3 training)

| Sub-score | What it measures |
|---|---|
| `ionization_score` | Signal intensity above instrument noise floor at typical 1-10 µM analyte |
| `complex_preservation_score` | For assemblies, fraction of total ion current in the intact complex vs sub-complexes |
| `csd_quality_score` | Width and gaussian-ness of the charge envelope (UniDec-style quality criterion) |
| `heterogeneity_score` | Controlled vocabulary entry: glycan-driven, oxidation, truncation, conjugation distribution (DAR), unknown |

`0 = unusable, 1 = poor, 2 = good, 3 = excellent`

## 4. Negative-example sub-schema (when `outcome_label >= 4`)

Negative records carry an additional structured field `failure_record` (JSON sidecar) with these fields:

| Field | Type | Description |
|---|---|---|
| `failure_mode` | string (CV) | One of the 7 failure modes below |
| `failure_confidence` | string | `high` / `medium` / `low` (based on number of independent attempts and operators) |
| `attempts_log` | list of dicts | Each attempt records buffer, concentration, instrument, source temperature, activation V, operator |
| `is_protein_intrinsic_failure` | bool | Curator's judgment: is the failure a property of the molecule or of the conditions tried |
| `orthogonal_evidence` | string | SEC, DLS, mass photometry, AUC results that support the failure interpretation. Free text + DOI. |
| `would_retry_with` | list of strings (CV) | Suggested alternative conditions |

### `failure_mode` controlled vocabulary

| Term | Definition |
|---|---|
| `no_ionization` | No analyte signal above noise after at least three buffer/concentration conditions tried |
| `denatured_signal_only` | Signal present only when protein is acid- or organic-denatured; the fold or complex is gone |
| `aggregation_dominant` | Broad high-m/z hump, no resolvable charge envelope; orthogonal evidence (DLS, SEC) of solution-phase aggregation |
| `gas_phase_dissociation` | Complex falls apart in source/transfer regardless of activation energy reduction |
| `uninterpretable_heterogeneity` | Signal present, charge states present, but mass cannot be assigned within reasonable error because of overlapping proteoforms |
| `csd_uninformative` | One or two charge states only, or signal at solvent peaks |
| `fragmentation_uncontrolled` | Survives ESI but cannot be activated/desolvated without complete fragmentation |

### Training inclusion rule for negatives

Only include in v0.3 loss function if **ALL** of:
- `is_protein_intrinsic_failure = true`
- `failure_confidence >= medium`
- `len(attempts_log) >= 3` (at least three distinct buffer/concentration conditions attempted)

Everything else goes into a held-out "soft negative" pool. Useful for semi-supervised training later, but not in the v0.3 loss function.

## 5. Controlled vocabulary: `protein_class`

| Term | When to use |
|---|---|
| `antibody` | Monoclonal antibody (mAb), full IgG or fragment |
| `bispecific` | Bispecific antibody |
| `adc` | Antibody-drug conjugate |
| `enzyme` | Globular enzyme (use specific subclass when known: kinase, phosphatase, etc.) |
| `kinase` | Protein kinase specifically |
| `phosphatase` | Protein phosphatase |
| `gpcr` | G-protein coupled receptor |
| `membrane protein` | Any other integral membrane protein |
| `ion channel` | Ion channel specifically |
| `transporter` | Membrane transporter |
| `viral capsid` | Viral capsid protein or assembly |
| `chaperone` | Protein chaperone |
| `idp` | Intrinsically disordered protein |
| `transcription factor` | Transcription factor |
| `complex` | Defined hetero-complex (use stoichiometry in `oligomeric_state_expected`) |
| `tetramer` / `dimer` / `hexamer` | Specific homo-oligomers |
| `small peptide` | Sequence length < 30 aa (conotoxins, bioactive peptides) |
| `protein` | Generic globular protein, use only when more specific category does not apply |

**Audit note**: 86 of the 260 newly extracted PDB records are currently classified as generic `protein`. These should be re-classified during v0.4 cleanup using GO terms from UniProt.

## 6. Controlled vocabulary: `instrument_psims_id`

Use HUPO PSI-MS CV terms: https://github.com/HUPO-PSI/psi-ms-CV

Common instruments in native MS:

| PSI-MS CV | Instrument |
|---|---|
| `MS:1003029` | Q Exactive UHMR (Thermo, Orbitrap) |
| `MS:1002523` | Q Exactive HF (Thermo) |
| `MS:1002416` | Synapt G2-Si (Waters) |
| `MS:1003123` | Cyclic IMS (Waters) |
| `MS:1003132` | timsTOF Pro (Bruker) |
| `MS:1003229` | timsTOF SCP (Bruker) |
| `MS:1003092` | solariX FTMS (Bruker) |

Add to this table as new instruments appear in the dataset.

## 7. Controlled vocabulary: `membrane_mimetic`

For membrane proteins only. Based on Ruotolo 2025 Anal. Chem. taxonomy.

| Term | Description |
|---|---|
| `none` | Membrane protein analyzed without mimetic (rare; usually fails) |
| `c8e4` | C8E4 detergent |
| `ddm` | n-Dodecyl-β-D-maltoside |
| `lda` | Lauryldimethylamine N-oxide |
| `ogp` | Octyl glucopyranoside |
| `polyamine_detergent` | Polyamine-functionalized detergent (Heck 2023) |
| `nanodisc_msp1d1` | Nanodisc with MSP1D1 scaffold |
| `nanodisc_msp1e3d1` | Nanodisc with MSP1E3D1 scaffold |
| `nanodisc_saposin` | Saposin-bound nanodisc |
| `peptidisc` | Peptidisc preparation |
| `smalp` | Styrene-maleic acid lipid particle |

## 8. Migration plan from existing 232+260 dataset

Existing files:
- `positives_with_sequences.json` (69 records, label=1)
- `negatives_with_sequences.json` (64 records, label=0, mostly proxy)
- `expansion_with_sequences.json` (99 records, mixed)
- `extracted_new_2026-05-01.json` (260 records, label=1, all PDB-derived)

Migration script (to be written): `migrate_to_v0.3_schema.py`

Mapping rules:
- All existing `label=1` → `outcome_label=1` (`clean_native`) by default; flag for human review
- All existing `label=0` from `negatives_with_sequences.json` → `outcome_label=NULL` and exclude from v0.3 training (these are proxy negatives, not real failures)
- All extracted PDB records → `outcome_label=2` (`interpretable_native`) by default (PDB inclusion implies the technique worked but spectra quality varies)
- Add `extraction_date`, `evidence_url_or_doi`, `submitter_id` columns
- Compute `sequence_sha256` for all rows

## 9. Data inclusion criteria for v0.3 training

A record is INCLUDED in v0.3 training if:

- `sequence_length >= 30` (excludes the 3 conotoxin peptides until peptide model variant exists)
- `sequence` contains only standard 20 amino acids
- `mw_kda` is present and > 0
- Either `outcome_label IN (1, 2, 3)` (positive/usable) OR `outcome_label IN (4, 5)` AND meets negative inclusion criteria from section 4
- `protein_class` is set
- `evidence_url_or_doi` is set OR record is in legacy 232 set

Targets for v0.3:
- ≥ 600 included positive records (currently 138 + 260 = 398; need ≥ 200 more)
- ≥ 50 strict negative records (currently 0 strict; need to extract from EuropePMC supplements)

## 10. Versioning

- Schema version is in the filename and in a top-level `schema_version` field of the dataset.
- Breaking changes increment the major version.
- New optional fields increment the minor version.
- Current: v0.1 (this document).

## 11. Open questions for Brhanu

1. **Are conotoxins (sequence < 30 aa) in or out of scope?** Currently 3 records. They are real native MS analytes but functionally peptides. Recommend out for v0.3, in scope for a future "small molecule / peptide" model variant.
2. **How to handle the 64 proxy negatives?** Recommend excluding from v0.3 training as flagged above. They are not real failures and may be hurting model calibration.
3. **Should `outcome_label` for PDB-derived records default to 1 or 2?** Default 2 is honest (PDB inclusion does not guarantee a perfect spectrum). Default 1 maximizes training data. Recommend 2 with the option to reclassify after spot-checking 30-50 records by reading their source papers.
4. **Should `protein_class` for the 86 generic 'protein' records be refined now or in v0.4?** Recommend v0.4 cleanup using UniProt GO terms (script-based, not manual).
5. **Does Trisam want to be the curator of failure_mode controlled vocabulary?** He has the domain expertise to maintain it.

## 12. Sources backing this design

- Heck et al., "Standard Proteoforms and Their Complexes for Native Mass Spectrometry," JASMS 2019. https://pmc.ncbi.nlm.nih.gov/articles/PMC6592724/
- Tamara et al., "High-Resolution Native Mass Spectrometry," Chem. Rev. 2022. https://pubs.acs.org/doi/10.1021/acs.chemrev.1c00212
- Marty et al., "Parsimonious Charge Deconvolution for Native MS," JPR 2018 (UniDec). https://pmc.ncbi.nlm.nih.gov/articles/PMC5838638/
- Ruotolo group, "Membrane mimetics evaluation for native IM-MS," Anal. Chem. 2025. https://pubs.acs.org/doi/10.1021/acs.analchem.4c06629
- Beveridge & Barran, "Are CSDs a Reliable Tool for IDPs?" JASMS 2016. https://link.springer.com/article/10.1007/s13361-016-1490-1
- HUPO PSI-MS Controlled Vocabulary. https://github.com/HUPO-PSI/psi-ms-CV
- CTDP Native Top-Down Project. https://www.topdownproteomics.org/initiatives/native-top-down-project/
