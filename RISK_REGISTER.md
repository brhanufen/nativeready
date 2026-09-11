# NativeReady — Startup Risk Register

*Why startups fail (evidence-based), NativeReady's specific exposure to each, and the mitigation in progress. For internal use and for answering the investor question "what keeps you up at night?"*

Source base: CB Insights startup post-mortem analysis (founders' own write-ups of why their company died), corroborated by academic entrepreneurship research and failure-rate statistics (~90% of startups fail; failure concentrates in years 2 to 5, when initial money runs out before traction lands). Causes sum to more than 100% because most failures have several at once.

---

## The three root categories

1. **Demand risk** — do people actually want it, enough to pay? (Hardest to fix, most fatal.)
2. **Execution risk** — can this team build and sell it? (Team, product, go-to-market, pricing.)
3. **Financing / timing risk** — is there enough money and runway to reach the next proof point?

Most companies that die had a demand problem that a financing problem exposed: they ran out of money precisely because demand was never strong enough to justify the next check.

---

## The ranked causes (from the post-mortem data)

| # | Cause | Frequency | Category |
|---|-------|-----------|----------|
| 1 | Ran out of cash / could not raise next round | ~38% | Financing |
| 2 | No market need | ~35% | Demand |
| 3 | Got outcompeted | ~20% | Execution |
| 4 | Flawed business model | ~19% | Execution |
| 5 | Regulatory / legal challenges | ~18% | External |
| 6 | Pricing / cost problems | ~15% | Execution |
| 7 | Wrong team / missing key skill | ~14% | Execution |
| 8 | Product mistimed (too early or too late) | ~10% | Timing |
| 9 | Poor product / usability | ~8% | Execution |
| 10 | Founder-investor disharmony, bad pivot, burnout | ~5-7% each | Execution |

---

## NativeReady's specific exposure (honest)

### Highest exposure — No market need (cause #2)
This is the number one risk, and it is known. Real usage today is 1 logged outcome (1 worked, 0 failed; live `/feedback/stats`, 2026-09-10). There are downloads (~836 all-time, PyPI) but no demonstrated willingness to pay. The entire "land a design-partner lab that logs real outcomes" plan *is* the demand-validation experiment. Everything hinges on it.
- **Mitigation in progress:** outreach to native-MS labs; the goal metric is one active lab logging real outcomes, not download count. This single milestone is the gate on everything else.

### Managed — Running out of cash (cause #1)
Manageable because the company is capital-light by design: no wet lab, small pre-seed, a Nebraska non-dilutive stack (state grants, MOVE, Invest Nebraska Seed match). That deliberately extends runway cheaply and hedges the most common killer.
- **Mitigation in progress:** ~$500K pre-seed sized to ~18 months; capital-efficient software model; multiple non-dilutive sources in parallel.

### Partially exposed — Wrong team / missing skill (cause #7)
Currently a technical founder recruiting a CEO. The gap being filled openly (commercial / go-to-market) is exactly the skill the data says technical-founder startups most often lack. Naming it and recruiting for it is the correct response, not a weakness to hide.
- **Mitigation in progress:** actively recruiting a full-time CEO to lead commercial, partnerships, and fundraising.

### Moderate — Got outcompeted (cause #3)
Larger AI-bio players could build a triage feature. The defense is the data moat (un-scrapeable logged outcomes, because public data is success-biased) plus workflow lock-in. Both only become real once labs are actually logging, so this shares a dependency with the demand risk.
- **Mitigation in progress:** build the logged-outcome flywheel early; become the triage step in the workflow so switching cost grows.

### Timing bet — Product mistimed (cause #8)
The "Why Now" thesis is that two waves just crossed: protein foundation models maturing, and native MS expanding into biopharma. If correct, timing is a tailwind. If real demand is still a few years out, it is a risk.
- **Mitigation in progress:** stay capital-light so the company can survive being early; let real usage data confirm or correct the timing assumption.

### Low exposure today — Regulatory (cause #5)
Not a factor at the current stage (a research/triage software tool, not a regulated diagnostic). Would only become relevant far in the future if the product moved toward clinical or QC-of-record use.

---

## The one-sentence version

The evidence says most startups die from **no real demand, revealed when the money runs out.** For NativeReady, both collapse into the same single task: **get one real lab using the tool and logging outcomes.** That is not a coincidence; it is why that milestone gates everything.

---

## How to answer an investor who asks "what is your biggest risk?"

"Demand. We are capital-light, so cash is not the thing that kills us early. The real risk is proving that labs will use this in their workflow and log outcomes. That is exactly why our pre-seed milestone is landing two or three design-partner labs actively logging real data, not chasing downloads. Everything, the moat, the pricing, the next round, follows from that one proof point."
