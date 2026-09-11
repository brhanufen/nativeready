# NativeReady — Current State & Vision (5-minute briefing)

*Synthesized 2026-09-10 from the repo (README, model card, `backend/heterogeneity.py`), the strategy docs (`RISK_REGISTER.md`, `DESIGN_PARTNER_LOGGING_PLAYBOOK.md`), `startup vision/` (`NativeReady_Vision.docx`, `NativeReady_Roadmap_to_Acquisition.docx`, `AGENTIC_VISION_Phase3.md`, the three `VISION_*` market docs), the flow figures, and `comparables/`. Every market number carries its source + confidence as stated in those docs. Where documents disagree, the conflict is flagged (see the last section).*

---

## 1. What NativeReady is today

**In one line:** paste a protein sequence, and NativeReady predicts — in seconds, before any experiment — whether that protein will give usable **native mass spectrometry (native MS)** data, and flags why it might not.

It is a live, product-complete triage tool. It is **not** yet a company with users or revenue (see §5).

**Who it's for today:** bench scientists, native-MS / MS-core facilities, and biopharma characterization teams who must decide which proteins are worth committing scarce, expensive instrument time to. (A failed native-MS run costs ~$5–15K in instrument time + weeks of rework on an often irreplaceable sample — the sourced pain figure behind the product.)

**Live surfaces:** web tool (`nativeready.bio`), a REST API (Railway), and a Python SDK (`pip install nativeready`). Built-in `/feedback` endpoint records experimental outcomes (worked / failed + conditions). Preprint: bioRxiv 2026.05.03.722506 (CC-BY). Open source (MIT), public data only.

### The two-layer technology

**Layer 1 — biophysical + ML suitability scorer** (`predictor_v3`, model `0.4-esm2-glyco-tm`).
47 sequence-derived features — BioPython physicochemical + amino-acid composition, plus v0.4 N-glycosylation-sequon and transmembrane-topology features — optionally concatenated with ESM-2 (`esm2_t33_650M_UR50D`) embeddings. Base estimator: a scikit-learn RandomForest with isotonic calibration, plus a nearest-neighbor **out-of-distribution flag**. Output: a 0–100 calibrated score, a label (Excellent→Unsuitable), a confidence interval, risk factors, and buffer/prep recommendations.

**Layer 2 — physics heterogeneity / resolvability calculator** ("Method 1," `backend/heterogeneity.py`).
This is the genuinely differentiated piece: **no ML, no fitted parameters — pure arithmetic.** For glycoproteins and antibody-drug conjugates it *simulates the proteoform mass envelope* (glycoform convolution from published monosaccharide masses, or a binomial DAR ladder for ADCs) and then *computes whether those proteoforms can be resolved* at a given instrument resolving power (R = m/FWHM, ~50%-valley criterion, charge-independent to first order). It returns a heterogeneity-risk level + reason, and discloses its assumptions on every call (75% sequon occupancy, 1% detection threshold, O-glycosylation not modeled). It is a **heterogeneity detector, not a success predictor** — it will flag heavily glycosylated proteins (e.g. SARS-CoV-2 Spike) that are nonetheless studied by native MS. This layer is original technical work, not a wrapper around a foundation model.

### Honest performance & scope (from the model card / README)
- Trained on **n = 635 proteins (538 positive, 97 negative)** — but only **~3 of the 97 negatives are evidence-based real failures**; 64 are random Swiss-Prot proxies, ~30 are property-targeted hard cases.
- Cluster-aware ROC-AUC ≈ **0.835 ± 0.029** (current v0.4 combined model; stratified 0.870). The **0.869** figure sometimes seen is the earlier **v0.3** number the bioRxiv preprint (dated 2026-05-03, before v0.4 shipped 2026-05-12) was built on — version evolution, not a discrepancy; the drop is the deliberate cost of removing a partially-memorized negative.
- **Major caveat, stated in the repo:** a size confound. Negatives are systematically longer (median 777 vs 288 aa); sequence length alone reaches AUC 0.769, so much of the apparent skill is a length artifact, and ESM-2 adds little over the biophysical features under the honest split.
- **Therefore NativeReady is a positive-suitability triage tool, not a validated failure detector.** Failure-detection performance is honestly reported as "not yet established," because the true native-MS failure rate can't be identified from this data. This limitation *is* the reason the whole company hinges on collecting real outcomes (§3–5).

**Figure — `flow_today.png`:** a linear pipeline, `PROTEIN SEQUENCE → NativeReady → SCORE + FLAGS + ADVICE → YOU DECIDE what to run`, captioned "one protein · one technique (native MS) · a scientist checks it." This is the accurate picture of *today*.

---

## 2. The vision — the decision layer for protein design

**The thesis:** the field flipped which step is scarce. Generative models (ESM, AlphaFold, RFdiffusion; design-first labs like Chai, Profluent) made *designing* proteins cheap and effectively unlimited. *Testing* them — expression, purification, characterization — is still slow, costly, and sample-limited. Every design pipeline now hits the same wall: **which of these thousands of designs do we actually make and test?**

**NativeReady's vision is to be the layer that answers that** — the triage/decision node between "design" and "build," generalized beyond native MS. Native MS is explicitly framed as **the wedge, not the destination.** The same core capability (predict from sequence whether an expensive, sample-consuming experiment will succeed) generalizes into the position the whole field is about to need.

The path is staged, each stage a standalone product feeding the next:
1. **Wedge (now):** native-MS suitability triage.
2. **Technique-agnostic triage:** the same core extended to cryo-EM, crystallography, HDX-MS, SEC-MALS — the pre-experiment decision layer across structural biology.
3. **"Will it work?" → "What makes it work?":** predict the *recipe* (buffer, construct, prep, instrument settings), turning the warning light into an autopilot — this is where the accumulated conditions data becomes a product.
4. **The agentic decision layer:** the API an autonomous design-build-test (DBTL) loop calls before spending physical resources.

**Why it's bigger:** a single-technique tool is a feature; the generalized decision layer that owns a one-of-a-kind outcome dataset is a *category position* — "the throttle on the single most expensive step." North star is **acquisition** (instrument makers Thermo/Waters/Bruker/SCIEX; informatics platforms Dotmatics/Genedata/Benchling-class; or the agentic-lab/Adaptyv-class wave).

**Figure — `flow_vision.png` (NEAR-TERM):** `AI DESIGNS (1,000s) → NativeReady DECISION LAYER ("which are worth testing?") → TEST the few / SKIP the rest → WET LAB`, with a "data feeds back" arrow closing the loop. **Figure — `flow_agentic.png` (LONG-TERM, Phase 3):** `YOUR GOAL → AGENT PLANS → YOU APPROVE (the honest human-in-the-loop gate) → PARTNER LAB runs it → AGENT LEARNS → next batch`. Note the deliberate discipline in both: NativeReady stays software; it does **not** own the wet lab (`AGENTIC_VISION_Phase3.md` explicitly rejects owning/deep-wiring a lab as capital-heavy, and keeps a human-approval gate rather than auto-retraining).

---

## 3. How today connects to the vision — the bridge

There is exactly **one metric that runs through every stage: real, conditions-annotated experimental outcomes logged by real users** (worked / failed + buffer, construct, expression system, instrument, resolution, failure mode).

That single number is simultaneously:
- **today's traction** (the thing that moves the company off "idea with a tool"), and
- **the vision's moat** (the un-scrapeable dataset that powers the generalized decision layer).

So the roadmap is *not* "build the product, then build the vision." It is: **grow one dataset, and let it carry the company from wedge → decision layer → acquisition.** Every design-partner lab = the first proof point *and* the first source of the technique-agnostic data moat. Every logged outcome = near-term traction *and* a permanent compounding asset. The structured-conditions logging that exists today is both an honesty feature and the foundation of the future "what makes it work" product. This is why the plan doesn't compete with the vision — it *is* the first brick of it.

---

## 4. The moat (one paragraph)

The defensibility is **not the model** — foundation models are shared and improving for everyone. It is a **data + workflow flywheel** (`flow_moat.png`): labs use NativeReady as the triage step before they commit instrument time → they log real outcomes (worked/failed + conditions) → that dataset — which *cannot be scraped from the literature, because failures are never published* — trains predictions no one else has → better predictions and deeper workflow embedding raise switching cost → more usage → more outcomes. Two reinforcing moats: a **data moat** (un-scrapeable, failures-included, compounding) and a **workflow moat** (the tool becomes the load-bearing triage step, so leaving gets more expensive over time). The original physics-based failure-mode simulators (Layer 2) add a second, non-commodity technical asset on top.

---

## 5. The honest state

- **Stage:** late pre-seed — product-complete, traction-pending. Pre-entity, pre-revenue.
- **Traction, stated without inflation:**
  - **Downloads: ~836 all-time PyPI** (non-mirror; a moving number) — a *distribution/interest* signal, explicitly **NOT users.**
  - **Active usage: 1 real logged outcome** (1 worked / 0 failed; live `/feedback/stats`, 2026-09-10) out of 3 feedback records — effectively pre-usage. **0 active design-partner labs.**
  - **Revenue: $0.**
- **Honest value-if-sold-today:** founder + IP + a working product ≈ a **~$0–5M acqui-hire** — which is exactly why the plan is *not* to sell now.
- **Capital plan:** capital-light by design (no wet lab); ~$500K pre-seed sized to ~18 months, leaning on a Nebraska non-dilutive stack. Team: a technical founder (Brhanu Znabu) actively recruiting a commercial CEO — naming the most common technical-founder gap rather than hiding it.
- **The single gating milestone:** **land a design-partner lab that logs real, conditions-annotated outcomes.** Per `RISK_REGISTER.md`, the #1 risk (and the #2 startup-death cause, "no market need") collapses into this one test; per the playbook, *one committed design partner* is the make-or-break signal (the investor-facing framing stretches it to "two or three"). Everything — the moat, pricing, the next round, the whole wedge→vision path — follows from getting the outcome count off zero. It is cheaply testable now: a few good conversations with native-MS labs (UNL and one degree out), not $500K.

---

## Contradictions found — now reconciled (see `ERRATA_AND_CANONICAL_FACTS.md`)

*All items below were reconciled to canonical values on 2026-09-10; the editable markdown docs are fixed, and corrections for the binary `.docx`/`.pptx`/`.png` files are listed in `ERRATA_AND_CANONICAL_FACTS.md`.*

1. **Logged-outcome count — RESOLVED:** live `/feedback/stats` (2026-09-10) shows **1 real outcome** (1 worked / 0 failed), 3 feedback records total. The roadmap's "0" was stale; canonical = 1.
2. **Two different market-sizing frames** (this is the important one):
   - *Deck frame* (`pitch_template/Market_Sizing_Sources.md`, 2026-08-24): **TAM ~$27.6B proteomics** (BCC Research 2024, 13.0% CAGR — **VERIFIED, primary**) / **SAM ~$6.33B mass spectrometry** (MarketsandMarkets 2024 → $9.62B by 2030, 7.2% — **VERIFIED, primary**) / **SOM ~$5–15M** serviceable (bottom-up, **MODELED**, not a census: ~500–2,000 cores ×$5K + ~50–400 biopharma teams ×$15–40K).
   - *Vision frame* (`startup vision/VISION_MARKET_SIZE.md`, 2026-08-29): **TAM ~$18B protein characterization** (Fact.MR 2025 — **LOW-confidence syndicated proxy**) / **SAM ~$2.1B AI-drug-discovery software** (GMInsights 2025 — **LOW proxy**) / **SOM ~$0.4–1.25M/yr** (bottom-up; per-company price is an explicit **ASSUMPTION**).
   These size two different things (a mass-spec tool vs. a decision layer for AI protein design). `VISION_NUMBERS_VERIFICATION.md` (2026-08-29) is the reconciler: it marks the $27.6B/$6.33B figures VERIFIED and deck-safe, and rules the **"AI protein design" market size UNVERIFIED — do not present as an established market** (only aggregators size it). Developability spend is likewise **UNVERIFIED standalone**; use "biopharma analytical testing **services** ~$4.8–8.6B" only as a labeled proxy.
   - `VISION_NUMBERS_VERIFICATION.md` also flags the deck's own stale/unconfirmed sub-figures: proteomics "$33.6B (2024)" is superseded (M&M rebased to $36.48B, 2025); mass-spec "ResearchAndMarkets ~$6.9B" is not confirmed (that page re-syndicates $6.33B).
   - **Canonical choice (see ERRATA):** lead with the VERIFIED proteomics/MS funnel — TAM $27.6B / SAM $6.33B / SOM $5–15M; treat the AI-protein-design and characterization figures as narrative / adjacent-proxy only, never the headline TAM/SAM.
3. **Customer-count honesty:** `VISION_CUSTOMER_REALITY.md` states the true "external-triage" cohort (AI-design companies with no wet lab that must gate physical testing) is **~3–5 today**, not "dozens"; the larger ~15–20 vertically-integrated design companies are an *internal-tool* sale with a build-vs-buy objection; the concrete near-term channel is cloud-lab pipelines (Adaptyv/Tamarind/Phylo). The vision docs' looser "dozens of protein-design companies already exist" should be read against this.
4. **"Native MS is a standard validation instrument" → softened:** `VISION_NUMBERS_VERIFICATION.md` classifies native MS (and HDX-MS) as the **niche/specialized** members of the validation stack; SPR/BLI + cryo-EM are the routine ones. NativeReady's own technique is the niche one — present it as a specialized-but-growing entry point.
5. **"Boltz" as a company:** `NativeReady_Vision.docx` / `flow_vision.png` list Boltz alongside Chai and Profluent as design-first *companies*. Boltz is an **MIT open-source model, not a company** (per `VISION_NUMBERS_VERIFICATION.md`) — fine to cite as a design *model*, wrong to call it a company/customer.
6. **Version drift — RESOLVED:** canonical SDK = **v0.4.2** (`python-sdk/pyproject.toml`), dataset = **v7** (`dataset_combined_v7_2026-05-11.json`), served model = `0.4-esm2-glyco-tm`. Stale references corrected in README (footer had "v0.3"; layout pointed at the v4 dataset) and project memory ("v0.3.1 / v5"). The composition (n=635; 538/97; ~3 real failures) was always stable. These three version numbers legitimately differ (package vs served model vs dataset) — the only errors were the stale pointers.
7. **Roadmap valuation bands are scenario estimates, not forecasts** (the doc says so): ~$0–5M now → ~$10–40M (Stage 2, decision layer + $50–300K ARR) → ~$30–80M (Stage 3, $1–3M ARR + unique dataset) → $100–500M+ (full-platform ceiling). Acquisition comps *are* VERIFIED from SEC 10-Qs (Blue Reference $14M, Nonlinear Dynamics ~$23M, Core Informatics $94M, Finesse $220M [weak comp]; recent: Dotmatics/Siemens $5.1B, ACD-Labs/Revvity ~$72M+$8M, Olink/Thermo $3.1B) — but those were post-hoc data tools *without* a proprietary predictive dataset, which is the roadmap's argument for a premium.

---

### The company in three sentences
NativeReady is a live, capital-light, honestly-scoped software tool that predicts native-MS suitability from sequence (an ML scorer plus an original physics-based heterogeneity calculator), with real distribution (~820–836 downloads) but essentially no active usage (0–1 logged outcomes) and no revenue yet. Its vision is to generalize the same "predict before you spend the experiment" capability into **the decision layer every AI-protein-design pipeline calls before committing a physical test** — a category position whose moat is an un-scrapeable, failures-included outcome dataset that compounds with use. Today and the vision are the same bet: **land one design-partner lab that logs real outcomes**, and let that one dataset carry the company from wedge to decision layer to acquisition.
