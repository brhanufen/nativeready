# NativeReady — Canonical Facts & Errata (single source of truth)

*Created 2026-09-10 to reconcile version/number drift across the file set so every document tells one consistent story. When any doc disagrees with the table below, the table wins. Re-run the checks noted before quoting a moving number externally.*

---

## 1. Canonical facts (conform every doc to these)

| Fact | Canonical value | Authority / how verified |
|---|---|---|
| **Served model version** | `0.4-esm2-glyco-tm` (entry point `backend/predictor_v3.py`) | `model/MODEL_CARD.md` |
| **Python SDK / PyPI package version** | **`nativeready` v0.4.2** | `python-sdk/pyproject.toml` (`version = "0.4.2"`) — source of truth |
| **Canonical dataset** | **v7** — `data/dataset_combined_v7_2026-05-11.json** | newest dataset file on disk; matches `MODEL_CARD.md` |
| **Dataset composition** | **n = 635 (538 positive, 97 negative; ~3 evidence-based real failures)** | stable across all docs; do not restate as 634 |
| **Cluster-aware combined AUC (current, v0.4)** | **0.835 ± 0.029** | `MODEL_CARD.md`, `README.md` |
| **Stratified combined AUC (current, v0.4)** | 0.870 ± 0.037 | `README.md` |
| **The 0.869 AUC figure** | is the **earlier v0.3** cluster-aware number, reported in the **bioRxiv preprint** (dated 2026-05-03, before v0.4 shipped 2026-05-12). Not a contradiction — version evolution. The drop 0.869→0.835 is the deliberate cost of removing the partially-memorized CB2 negative. | `README.md` v0.4.1 notes; preprint date vs v0.4 date |
| **Logged real outcomes** | **1 real outcome** (1 worked, 0 failed) + 2 not-tested = **3 total feedback records** | live `/feedback/stats` via `./check_outcomes.sh`, 2026-09-10 |
| **PyPI downloads** | **~836 all-time (non-mirror)** — a *moving* number; treat as "≈800+". Distribution signal, **NOT** users. | roadmap (most recent); reconcile the older "~820" to this |
| **Preprint** | bioRxiv 2026.05.03.722506 (CC-BY) | README |
| **Boltz** | an **MIT open-source model (Boltz-1/2), NOT a company.** Cite it as a design *model* alongside ESM/AlphaFold/RFdiffusion; never list it as a company/customer alongside Chai/Profluent. | `VISION_NUMBERS_VERIFICATION.md` §B |

### Note: three version numbers legitimately differ
The **served model** (0.4-esm2-glyco-tm), the **SDK package** (0.4.2), and the **dataset** (v7) version independently — that is normal and not an error. The errors are only the *stale references* below.

---

## 2. Canonical market frame (pick this; stop presenting two)

**Lead with the verified proteomics/mass-spec funnel as the addressable market *today*. Present "the decision layer for AI protein design" as the expansion narrative, NOT as a sized market.**

| Layer | Canonical figure | Source | Confidence |
|---|---|---|---|
| **TAM** | **~$27.6B** global proteomics | BCC Research 2024 (13.0% CAGR); Grand View $27.8B corroborates | **VERIFIED (primary)** |
| **SAM** | **~$6.33B** mass spectrometry (→ $9.62B by 2030, 7.2%) | MarketsandMarkets 2024 | **VERIFIED (primary)** |
| **SOM** | **~$5–15M** serviceable (~$0.15–7.8M obtainable) | bottom-up: ~500–2,000 cores ×$5K + ~50–400 biopharma teams ×$15–40K | **MODELED** (not a census; per-company price is an assumption) |

**Do NOT present as a sized market (all UNVERIFIED / low-confidence proxy):**
- "AI protein design market" ($1.5B-ish, aggregator-only) — **UNVERIFIED**; use only as narrative.
- Protein-characterization "$18B" (Fact.MR) and AI-drug-discovery-software "$2.1B" (GMInsights) from `VISION_MARKET_SIZE.md` — **LOW-confidence syndicated proxies**; if used at all, label as "adjacent proxy for where the vision grows into," never as the headline TAM/SAM.
- "Developability spend" — no standalone figure; only "biopharma analytical testing **services** ~$4.8–8.6B" as a labeled proxy.

**Stale sub-figures to drop wherever they appear:** proteomics "$33.6B (2024)" (M&M rebased to $36.48B/2025); mass-spec "ResearchAndMarkets ~$6.9B" (that page re-syndicates $6.33B). **Native MS** is a *niche/specialized* validation technique, not a "standard" one (SPR/BLI + cryo-EM are the routine stack) — state it that way.

**Acquisition comps** (VERIFIED from acquirers' SEC 10-Qs): Blue Reference $14M (Waters, 2012), Nonlinear Dynamics ~$23M (Waters, 2013), Core Informatics $94M (Thermo, 2017), Finesse $220M (Thermo, 2017 — weak comp, bioprocessing hardware). Recent: Dotmatics/Siemens $5.1B (2025), ACD-Labs/Revvity ~$72M+$8M earnout (2026), Olink/Thermo $3.1B (2023→24). Roadmap valuation bands are **reasoned scenario estimates, not forecasts.**

---

## 3. Per-file corrections

### Files already fixed in this pass (markdown, editable)
- **`README.md`** — footer version line and the "canonical" dataset pointer updated to current (v0.4 / v0.4.2 SDK / v7 dataset). ✔
- **`RISK_REGISTER.md`** — download figure aligned to ~836; "~1 logged outcome" confirmed correct (now exactly 1). ✔
- **`NATIVEREADY_CURRENT_AND_VISION.md`** — outcome count set to 1, downloads to ~836, SDK to v0.4.2, dataset to v7, AUC framing corrected to version-evolution, contradictions section marked reconciled. ✔

### Files NOT auto-editable (binary .docx / .pptx / .png) — apply these on next regeneration
- **`startup vision/NativeReady_Roadmap_to_Acquisition.docx`** — Stage 0 says "0 real worked/failed outcomes logged" → change to **"1 real outcome logged (1 worked / 0 failed); 3 total feedback records."** Also SDK "v0.4.2" is already correct here; downloads "~836" is fine (make it the canonical figure).
- **`startup vision/NativeReady_Vision.docx`** — text "design-first companies like Chai Discovery and **Boltz**" → Boltz is a model, not a company; reword to "design-first companies like Chai and Profluent, and open models like ESM/AlphaFold/RFdiffusion/Boltz."
- **`startup vision/NativeReady_Vision_Figure.docx`** and **`startup vision/flow_vision.png`** — the "AI DESIGNS" box lists "ESM · AlphaFold · RFdiffusion · Boltz + labs (Chai, Profluent)". This is acceptable (Boltz shown as a model, Chai/Profluent as labs) — no change needed, but keep Boltz on the *models* side, never the labs side.
- **`startup vision/flow_moat.png`** — "~1 outcome logged to date" is now exactly correct; optionally update to "1 outcome logged (2026-09)."
- **Pitch decks (`pitch_template/*.pptx`, `.pdf`)** — ensure the market slide uses the canonical funnel in §2 ($27.6B / $6.33B / $5–15M), not the $18B/$2.1B proxy frame; ensure SDK shown as v0.4.2; ensure any downloads figure matches ~836.

### External memory / notes
- The assistant's project memory previously listed "PyPI SDK v0.3.1 … dataset v5" — **stale**; corrected to SDK v0.4.2, dataset v7, 1 logged outcome.

---

## 4. Re-verify before quoting externally
- **Outcome count:** `./check_outcomes.sh` (live `/feedback/stats`) — moving number, this is the real scoreboard.
- **Downloads:** pypistats.org/packages/nativeready — all-time is a moving number.
- **Market comps/sizes:** re-ground against `startup vision/VISION_NUMBERS_VERIFICATION.md` before any slide or data room; it is the reconciler of record for market figures.
