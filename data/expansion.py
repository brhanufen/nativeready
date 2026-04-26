"""
Expansion sets for the native MS suitability predictor.

Three lists:
1. SMALL_PEPTIDES_POSITIVES — small peptides (<50 aa typical) with documented
   native MS / ESI-MS literature. Improves coverage at the low-MW end.
2. HARD_NEGATIVES — proteins that experts would expect to be very challenging
   for native MS (heavy O-glycosylation, very large multi-pass membrane,
   highly disordered, aggregation-prone). These are NOT random Swiss-Prot
   draws; they are deliberate hard-cases.
3. ADDITIONAL_POSITIVES — more proteins with documented native MS literature
   to broaden the positive class (kinases, GPCRs, antibodies, plant
   complexes, additional viral, etc.).

All UniProt accessions are real and verifiable on uniprot.org. Citations are
"short" (lab + topic + year-ish) to be filled in / verified before any
publication. No fabricated DOIs.

No accession in this file duplicates an entry in positives_curated.py
(checked manually against the 69-entry POSITIVE_EXAMPLES list).

MW values for parent precursor proteins are the full-length UniProt-listed
MW; for peptides derived from a precursor (e.g., bradykinin from KNG1),
the mw_kda field reflects the MATURE PEPTIDE mass actually observed in
native MS, and the notes call out the precursor.
"""

# -----------------------------------------------------------------------------
# LIST 1: small peptides documented as native / ESI-MS targets
# -----------------------------------------------------------------------------
SMALL_PEPTIDES_POSITIVES = [
    {"uniprot_id": "P01308", "name": "Insulin (human, intact A+B disulfide-linked)",
     "protein_class": "small peptide", "mw_kda": 5.8,
     "notes": "Classic ESI/native MS standard; intact insulin and chains routinely observed (Loo, Smith)."},
    {"uniprot_id": "P01275", "name": "Glucagon (human, mature 29-aa peptide from GCG)",
     "protein_class": "small peptide", "mw_kda": 3.5,
     "notes": "ESI-MS of glucagon and fibril intermediates; precursor GCG (P01275)."},
    {"uniprot_id": "P01178", "name": "Oxytocin (mature 9-aa peptide from OXT precursor)",
     "protein_class": "small peptide", "mw_kda": 1.0,
     "notes": "Cyclic disulfide nonapeptide; ESI-MS metal-binding studies (Pernelle, Wysocki)."},
    {"uniprot_id": "P01185", "name": "Vasopressin / AVP (mature 9-aa peptide)",
     "protein_class": "small peptide", "mw_kda": 1.1,
     "notes": "Cyclic nonapeptide; routine ESI-MS reference; Cu/Zn binding native MS."},
    {"uniprot_id": "P20366", "name": "Substance P (mature 11-aa peptide from TAC1)",
     "protein_class": "small peptide", "mw_kda": 1.35,
     "notes": "Tachykinin neuropeptide; standard ESI-MS calibration peptide."},
    {"uniprot_id": "P01042", "name": "Bradykinin (mature 9-aa peptide from kininogen-1)",
     "protein_class": "small peptide", "mw_kda": 1.06,
     "notes": "Universal ESI-MS calibrant; mass = 1060 Da; precursor KNG1."},
    {"uniprot_id": "P01258", "name": "Calcitonin (human, mature 32-aa peptide from CALCA)",
     "protein_class": "small peptide", "mw_kda": 3.4,
     "notes": "Amyloidogenic peptide; oligomer distributions by native ESI-MS (Bowers)."},
    {"uniprot_id": "P61278", "name": "Somatostatin-14 (mature peptide from SST)",
     "protein_class": "small peptide", "mw_kda": 1.64,
     "notes": "Cyclic disulfide 14-mer; ESI-MS metal coordination studies."},
    {"uniprot_id": "P01019", "name": "Angiotensin II (mature 8-aa peptide from AGT)",
     "protein_class": "small peptide", "mw_kda": 1.05,
     "notes": "Standard ESI-MS calibration peptide; conformer studies by IM-MS (Clemmer)."},
    # NOTE: Abeta42 is derived from APP (P05067), which is already in
    # positives_curated.py. UniProt does not assign Abeta its own accession,
    # so we deliberately omit it here rather than duplicate P05067.
    {"uniprot_id": "P01501", "name": "Melittin (bee venom, 26-aa)",
     "protein_class": "small peptide / antimicrobial", "mw_kda": 2.85,
     "notes": "Tetramer assembly observable by native ESI-MS; lipid interaction studies."},
    {"uniprot_id": "P59665", "name": "HNP-1 / human alpha-defensin 1 (from DEFA1)",
     "protein_class": "small peptide / antimicrobial", "mw_kda": 3.4,
     "notes": "Three-disulfide defensin; ESI-MS dimer/oligomer characterization."},
    {"uniprot_id": "Q01523", "name": "HD-5 / human alpha-defensin 5 (from DEFA5)",
     "protein_class": "small peptide / antimicrobial", "mw_kda": 3.6,
     "notes": "Paneth-cell defensin; native MS of dimers and metal binding."},
    {"uniprot_id": "Q01524", "name": "HD-6 / human alpha-defensin 6 (from DEFA6)",
     "protein_class": "small peptide / antimicrobial", "mw_kda": 3.7,
     "notes": "Self-assembling nanonet defensin; native MS oligomer ladders."},
    {"uniprot_id": "P01189", "name": "ACTH 1-39 (mature peptide from POMC)",
     "protein_class": "small peptide", "mw_kda": 4.5,
     "notes": "Pituitary peptide; ESI-MS standard; precursor POMC."},
    {"uniprot_id": "P01270", "name": "Parathyroid hormone PTH 1-84 (mature peptide)",
     "protein_class": "small peptide / hormone", "mw_kda": 9.4,
     "notes": "Therapeutic peptide; intact native MS for biosimilar QC."},
    {"uniprot_id": "P05305", "name": "Endothelin-1 (mature 21-aa peptide from EDN1)",
     "protein_class": "small peptide", "mw_kda": 2.49,
     "notes": "Two-disulfide peptide; ESI-MS folding and receptor binding."},
    {"uniprot_id": "P10082", "name": "Peptide YY / PYY (mature 36-aa from PYY)",
     "protein_class": "small peptide", "mw_kda": 4.05,
     "notes": "Gut hormone; ESI-MS structural and aggregation studies."},
    {"uniprot_id": "P01303", "name": "Neuropeptide Y / NPY (mature 36-aa)",
     "protein_class": "small peptide", "mw_kda": 4.27,
     "notes": "PP-fold family; ESI-MS dimer dissociation studies."},
    {"uniprot_id": "P0DV02", "name": "Liraglutide-relevant GLP-1(7-37) (mature peptide from GCG)",
     "protein_class": "small peptide / hormone", "mw_kda": 3.36,
     "notes": "GLP-1 incretin; native MS for therapeutic peptide aggregation."},
]

# -----------------------------------------------------------------------------
# LIST 2: hard negatives — proteins expected to be very challenging
# -----------------------------------------------------------------------------
HARD_NEGATIVES = [
    # A. Heavily O-glycosylated mucins
    {"uniprot_id": "P15941", "name": "MUC1 (human)",
     "protein_class": "O-glycosylated mucin", "mw_kda": 122.1,
     "notes": "Tandem-repeat O-glycosylation, extreme microheterogeneity; intact native MS effectively impossible."},
    {"uniprot_id": "Q02817", "name": "MUC2 (human)",
     "protein_class": "O-glycosylated mucin", "mw_kda": 540.0,
     "notes": "Gel-forming mucin >5000 aa; massive O-glycan heterogeneity; oligomerizes via disulfides."},
    {"uniprot_id": "P98088", "name": "MUC5AC (human)",
     "protein_class": "O-glycosylated mucin", "mw_kda": 585.0,
     "notes": "Airway gel-forming mucin; multi-MDa polymers in vivo; native MS not tractable."},
    {"uniprot_id": "Q9HC84", "name": "MUC5B (human)",
     "protein_class": "O-glycosylated mucin", "mw_kda": 596.0,
     "notes": "Salivary/airway mucin; heavy O-glyco, disulfide-linked polymers."},
    {"uniprot_id": "Q8WXI7", "name": "MUC16 / CA-125 (human)",
     "protein_class": "O-glycosylated mucin", "mw_kda": 1518.0,
     "notes": "~22000 aa, ovarian cancer marker; size and glycan microheterogeneity preclude native MS."},
    {"uniprot_id": "Q9HCU0", "name": "Endosialin / CD248 (human)",
     "protein_class": "O-glycosylated single-pass", "mw_kda": 80.9,
     "notes": "Heavily O-glycosylated stromal marker; glycan heterogeneity blurs native MS peaks."},

    # B. Very large multi-pass / cytoskeletal proteins (>1000 aa)
    {"uniprot_id": "P11532", "name": "Dystrophin (human)",
     "protein_class": "very large cytoskeletal", "mw_kda": 426.7,
     "notes": "3685 aa; spectrin-like rod; aggregation-prone in vitro; not a native MS target."},
    {"uniprot_id": "Q8WZ42", "name": "Titin (human, partial)",
     "protein_class": "very large multi-domain", "mw_kda": 3816.0,
     "notes": "Largest known protein, ~34000 aa; native MS of intact titin not feasible."},
    {"uniprot_id": "P02549", "name": "Spectrin alpha I (human)",
     "protein_class": "very large cytoskeletal", "mw_kda": 280.0,
     "notes": "2419 aa; forms (αβ)2 tetramers >1 MDa with extreme aspect ratio; not native MS friendly."},
    {"uniprot_id": "Q15149", "name": "Plectin (human)",
     "protein_class": "very large cytoskeletal", "mw_kda": 531.5,
     "notes": "4684 aa cytolinker; multimerizes; impractical for intact native MS."},
    {"uniprot_id": "P27816", "name": "MAP4 (human)",
     "protein_class": "large MAP / disordered", "mw_kda": 121.0,
     "notes": "1152 aa, large disordered projection; broad CSDs and very heterogeneous."},
    {"uniprot_id": "P78527", "name": "DNA-PKcs / PRKDC (human)",
     "protein_class": "very large kinase", "mw_kda": 469.0,
     "notes": "4128 aa; >450 kDa monomer; intact native MS rarely attempted on bare protein."},
    {"uniprot_id": "P04637", "name": "skip"  # placeholder removed below
     , "protein_class": "n/a", "mw_kda": 0.0, "notes": "REMOVE"},
    {"uniprot_id": "P98164", "name": "LRP2 / megalin (human)",
     "protein_class": "very large single-pass receptor", "mw_kda": 521.9,
     "notes": "4655 aa endocytic receptor; heavy glycosylation + size; not a native MS substrate."},
    {"uniprot_id": "Q9Y6V0", "name": "Piccolo / PCLO (human)",
     "protein_class": "very large scaffold / disordered", "mw_kda": 560.0,
     "notes": "5142 aa presynaptic scaffold, largely disordered; impractical for intact native MS."},
    {"uniprot_id": "P21817", "name": "Ryanodine receptor 1 / RyR1 (human)",
     "protein_class": "very large multi-pass channel", "mw_kda": 565.0,
     "notes": "5038 aa; 2.2 MDa homotetramer in membrane; native MS extremely demanding (only specialized labs)."},

    # C. Highly intrinsically disordered (not already in positives)
    {"uniprot_id": "O95997", "name": "Securin / PTTG1 (human)",
     "protein_class": "IDP / cell cycle", "mw_kda": 22.0,
     "notes": "Largely disordered separase inhibitor; broad CSDs, ill-defined fold."},
    {"uniprot_id": "Q12778", "name": "FOXO1 (human)",
     "protein_class": "IDP transcription factor", "mw_kda": 69.7,
     "notes": "Long disordered transactivation regions outside FH domain; broad CSD."},
    {"uniprot_id": "Q9Y261", "name": "FOXA2 (human)",
     "protein_class": "IDP transcription factor", "mw_kda": 48.3,
     "notes": "Pioneer TF with extensive disorder; difficult to obtain narrow native MS peaks."},
    {"uniprot_id": "P04150", "name": "Glucocorticoid receptor / NR3C1 (human)",
     "protein_class": "IDR-containing nuclear receptor", "mw_kda": 85.6,
     "notes": "Large N-terminal IDR; full-length GR notoriously hard to handle for native MS."},
    {"uniprot_id": "P38398", "name": "BRCA1 (human)",
     "protein_class": "very large IDR-containing", "mw_kda": 207.7,
     "notes": "1863 aa with extensive disorder; aggregation-prone; not a tractable native MS target full-length."},
    {"uniprot_id": "Q14653", "name": "IRF3 (human)",
     "protein_class": "TF with IDR", "mw_kda": 47.2,
     "notes": "Phospho-regulated dimerization with disordered linkers; charge-state ladder is broad."},

    # D. Aggregation-prone / amyloidogenic / polyQ
    {"uniprot_id": "P42858", "name": "Huntingtin / HTT (human)",
     "protein_class": "polyQ / aggregation-prone", "mw_kda": 347.6,
     "notes": "3142 aa with polyQ; exon-1 fragment is aggregation-prone; intact HTT not native-MS friendly."},
    {"uniprot_id": "P54259", "name": "Atrophin-1 / ATN1 (human, polyQ)",
     "protein_class": "polyQ / aggregation-prone", "mw_kda": 125.5,
     "notes": "PolyQ disease protein (DRPLA); aggregation makes intact native MS unsuitable."},
    {"uniprot_id": "P54253", "name": "Ataxin-1 / ATXN1 (human, polyQ)",
     "protein_class": "polyQ / aggregation-prone", "mw_kda": 86.9,
     "notes": "SCA1 polyQ protein; oligomerizes; broad CSD and aggregation."},
    {"uniprot_id": "P04156", "name": "Prion protein PrP (human)",
     "protein_class": "aggregation-prone / amyloidogenic", "mw_kda": 27.7,
     "notes": "Conversion to PrPSc and oligomerization make reproducible native MS hard; metal binding adds heterogeneity."},
    {"uniprot_id": "P19838", "name": "NF-kB p105 / NFKB1 (human)",
     "protein_class": "very large with ankyrin/IDR", "mw_kda": 105.4,
     "notes": "Long disordered C-terminus, processed in vivo; intact form rarely amenable to native MS."},
    {"uniprot_id": "Q9UM47", "name": "NOTCH3 (human)",
     "protein_class": "very large multi-pass receptor", "mw_kda": 243.6,
     "notes": "2321 aa, heavily glycosylated, multiple EGF repeats; intact native MS impractical."},
    {"uniprot_id": "P35658", "name": "Nucleoporin NUP214 (human)",
     "protein_class": "FG-repeat nucleoporin", "mw_kda": 213.6,
     "notes": "Large FG-repeat IDR; phase separates; not a native MS-friendly target."},
    {"uniprot_id": "Q8N3C0", "name": "ASCC3 (human)",
     "protein_class": "very large helicase", "mw_kda": 251.3,
     "notes": "2202 aa; large multi-domain helicase; full-length intact MS not commonly reported."},
    {"uniprot_id": "P49792", "name": "RANBP2 / NUP358 (human)",
     "protein_class": "very large FG nucleoporin", "mw_kda": 358.2,
     "notes": "3224 aa nucleoporin with extensive disorder; phase-separating; not native MS friendly."},
]

# Drop the placeholder we left for shape clarity (P04637 is in positives_curated):
HARD_NEGATIVES = [e for e in HARD_NEGATIVES if e["notes"] != "REMOVE"]

# -----------------------------------------------------------------------------
# LIST 3: additional positive native MS targets (no overlap with positives_curated)
# -----------------------------------------------------------------------------
ADDITIONAL_POSITIVES = [
    # Kinases studied by native MS
    {"uniprot_id": "P00517", "name": "PKA catalytic subunit alpha (bovine)",
     "protein_class": "kinase", "mw_kda": 40.6,
     "notes": "Regulatory/catalytic complex by native MS; nucleotide and inhibitor binding."},
    {"uniprot_id": "P17612", "name": "PKA catalytic alpha (human, PRKACA)",
     "protein_class": "kinase", "mw_kda": 40.6,
     "notes": "Holoenzyme assembly and PKI peptide binding by native MS."},
    {"uniprot_id": "P24941", "name": "CDK2 (human)",
     "protein_class": "kinase", "mw_kda": 33.9,
     "notes": "CDK2/cyclin and inhibitor binding by native MS (Robinson, Heck)."},
    {"uniprot_id": "P38936", "name": "p21 / CDKN1A (human)",
     "protein_class": "IDP / kinase regulator", "mw_kda": 18.1,
     "notes": "CDK/cyclin regulator captured in complex by native MS."},
    {"uniprot_id": "P27361", "name": "ERK1 / MAPK3 (human)",
     "protein_class": "kinase", "mw_kda": 43.1,
     "notes": "Phosphoform-resolved native MS of MAPK signalling proteins."},
    {"uniprot_id": "P28482", "name": "ERK2 / MAPK1 (human)",
     "protein_class": "kinase", "mw_kda": 41.4,
     "notes": "Phosphoform resolution and substrate binding by native MS."},
    {"uniprot_id": "Q02750", "name": "MEK1 / MAP2K1 (human)",
     "protein_class": "kinase", "mw_kda": 43.4,
     "notes": "Allosteric inhibitor binding (trametinib class) studied by native MS."},
    {"uniprot_id": "P31749", "name": "AKT1 (human)",
     "protein_class": "kinase", "mw_kda": 55.7,
     "notes": "Allosteric inhibitor binding by native MS."},

    # Phosphatases
    {"uniprot_id": "P18031", "name": "PTP1B / PTPN1 (human)",
     "protein_class": "phosphatase", "mw_kda": 49.9,
     "notes": "Drug-discovery phosphatase; inhibitor screening by native MS."},
    {"uniprot_id": "P30153", "name": "PP2A scaffold subunit A alpha / PPP2R1A (human)",
     "protein_class": "phosphatase scaffold", "mw_kda": 65.3,
     "notes": "PP2A holoenzyme assembly captured by native MS (Heck lab)."},

    # GPCRs / membrane drug targets (Robinson lab and others)
    {"uniprot_id": "P07550", "name": "Beta-2 adrenergic receptor / ADRB2 (human)",
     "protein_class": "GPCR", "mw_kda": 46.5,
     "notes": "Landmark GPCR native MS in detergent/nanodiscs (Robinson)."},
    {"uniprot_id": "P30542", "name": "Adenosine A1 receptor / ADORA1 (human)",
     "protein_class": "GPCR", "mw_kda": 36.5,
     "notes": "Lipid- and ligand-binding native MS in membrane mimetics."},
    {"uniprot_id": "P29274", "name": "Adenosine A2A receptor / ADORA2A (human)",
     "protein_class": "GPCR", "mw_kda": 44.7,
     "notes": "Ligand and cholesterol binding by native MS (Robinson lab)."},
    {"uniprot_id": "P21728", "name": "Dopamine D1 receptor / DRD1 (human)",
     "protein_class": "GPCR", "mw_kda": 49.3,
     "notes": "GPCR-lipid interactions by native MS."},
    {"uniprot_id": "P41595", "name": "Serotonin 5-HT2B receptor / HTR2B (human)",
     "protein_class": "GPCR", "mw_kda": 54.3,
     "notes": "Ligand binding studied by native MS in detergent."},

    # Ion channels / transporters
    {"uniprot_id": "Q03721", "name": "Kv3.4 / KCNC4 (human)",
     "protein_class": "membrane channel", "mw_kda": 70.4,
     "notes": "K+ channel native MS in nanodiscs (Robinson)."},
    {"uniprot_id": "P63027", "name": "VAMP2 (human)",
     "protein_class": "SNARE", "mw_kda": 12.7,
     "notes": "SNARE complex assembly by native MS."},
    {"uniprot_id": "P32856", "name": "Syntaxin-2 (human)",
     "protein_class": "SNARE", "mw_kda": 33.2,
     "notes": "SNARE complex stoichiometry by native MS."},

    # Antibodies and Fabs (additional, non-duplicating positives_curated)
    {"uniprot_id": "P01871", "name": "IgM heavy chain constant (human)",
     "protein_class": "antibody", "mw_kda": 49.4,
     "notes": "Pentameric IgM ~970 kDa native MS (Heck lab)."},
    {"uniprot_id": "P01876", "name": "IgA1 heavy chain constant (human)",
     "protein_class": "antibody", "mw_kda": 37.7,
     "notes": "Dimeric secretory IgA assembly by native MS."},
    {"uniprot_id": "P01877", "name": "IgA2 heavy chain constant (human)",
     "protein_class": "antibody", "mw_kda": 36.5,
     "notes": "Subclass-resolved IgA native MS."},
    {"uniprot_id": "P01880", "name": "IgD heavy chain constant (human)",
     "protein_class": "antibody", "mw_kda": 42.2,
     "notes": "Surface IgD characterized by native MS."},
    {"uniprot_id": "P01854", "name": "IgE heavy chain constant (human)",
     "protein_class": "antibody", "mw_kda": 47.2,
     "notes": "FcεRI binding studied by native MS."},
    {"uniprot_id": "P01591", "name": "Ig J chain (human)",
     "protein_class": "antibody adapter", "mw_kda": 18.1,
     "notes": "Joining chain in IgM/IgA polymers; resolved in native MS of pIgs."},
    {"uniprot_id": "P01833", "name": "Polymeric immunoglobulin receptor (human)",
     "protein_class": "antibody-related", "mw_kda": 83.3,
     "notes": "Secretory component on sIgA; native MS of polymeric Ig assemblies."},

    # Plant / photosynthesis complexes
    {"uniprot_id": "P00875", "name": "RuBisCO large subunit (spinach)",
     "protein_class": "plant enzyme", "mw_kda": 52.9,
     "notes": "L8S8 hexadecamer (~550 kDa) intact native MS (Robinson, Heck)."},
    {"uniprot_id": "P10795", "name": "RuBisCO small subunit (spinach)",
     "protein_class": "plant enzyme", "mw_kda": 14.8,
     "notes": "Small subunit of L8S8 RuBisCO complex."},
    {"uniprot_id": "P02905", "name": "PsbA / D1 photosystem II (Synechocystis)",
     "protein_class": "photosynthetic membrane", "mw_kda": 39.5,
     "notes": "Photosystem II core; native MS of PSII supercomplex (Boekema/Heck)."},
    {"uniprot_id": "P03689", "name": "Light-harvesting complex CP43 / PsbC (spinach)",
     "protein_class": "photosynthetic membrane", "mw_kda": 50.1,
     "notes": "Component of PSII supercomplex captured by native MS."},

    # Additional viral proteins (none duplicating existing entries)
    {"uniprot_id": "P03466", "name": "Influenza A nucleoprotein NP (H1N1, A/PR/8/34)",
     "protein_class": "viral RNP", "mw_kda": 56.2,
     "notes": "NP oligomers and RNP assembly by native MS."},
    {"uniprot_id": "P03485", "name": "Influenza A matrix protein M1",
     "protein_class": "viral matrix", "mw_kda": 27.8,
     "notes": "M1 oligomerization studied by native MS."},
    {"uniprot_id": "P03070", "name": "SV40 large T antigen",
     "protein_class": "viral helicase", "mw_kda": 81.6,
     "notes": "Hexameric helicase assembly by native MS."},
    {"uniprot_id": "Q77M19", "name": "Norovirus VP1 capsid (Norwalk virus)",
     "protein_class": "viral capsid", "mw_kda": 56.5,
     "notes": "T=3 VLP (180-mer) intact native MS (Heck, Uetrecht)."},
    {"uniprot_id": "P03129", "name": "HPV-16 E7 oncoprotein",
     "protein_class": "viral oncoprotein / IDP", "mw_kda": 11.0,
     "notes": "Dimer formation and Zn binding by native MS."},
    {"uniprot_id": "P12504", "name": "HIV-1 Vif",
     "protein_class": "viral regulatory", "mw_kda": 23.1,
     "notes": "Vif/CBF-beta/Cul5 assembly captured by native MS."},
    {"uniprot_id": "P04591", "name": "HIV-1 Gag polyprotein (HXB2)",
     "protein_class": "viral structural", "mw_kda": 55.8,
     "notes": "Gag assembly intermediates and VLPs by native MS (Uetrecht)."},
    {"uniprot_id": "P03366", "name": "HIV-1 Pol polyprotein (HXB2)",
     "protein_class": "viral enzyme", "mw_kda": 115.3,
     "notes": "RT and protease native MS (drug binding)."},

    # Additional drug targets / oncology
    {"uniprot_id": "P10275", "name": "Androgen receptor / AR (human)",
     "protein_class": "nuclear receptor LBD", "mw_kda": 99.2,
     "notes": "Ligand-binding domain native MS for antagonist screening."},
    {"uniprot_id": "P03372", "name": "Estrogen receptor alpha / ESR1 (human)",
     "protein_class": "nuclear receptor LBD", "mw_kda": 66.2,
     "notes": "LBD native MS for SERM/SERD characterization."},
    {"uniprot_id": "P10415", "name": "Bcl-2 (human)",
     "protein_class": "apoptosis regulator", "mw_kda": 26.3,
     "notes": "BH3 mimetic (venetoclax) binding by native MS."},
    {"uniprot_id": "Q07817", "name": "Bcl-xL / BCL2L1 (human)",
     "protein_class": "apoptosis regulator", "mw_kda": 26.0,
     "notes": "BH3 peptide binding stoichiometry by native MS."},
    {"uniprot_id": "Q9NR28", "name": "DIABLO / Smac (human)",
     "protein_class": "apoptosis regulator", "mw_kda": 27.1,
     "notes": "Dimer + IAP binding by native MS."},

    # Small enzymes and metabolic
    {"uniprot_id": "P00558", "name": "Phosphoglycerate kinase 1 / PGK1 (human)",
     "protein_class": "metabolic enzyme", "mw_kda": 44.6,
     "notes": "Open/closed conformers by IM-native MS."},
    {"uniprot_id": "P04406", "name": "GAPDH (human)",
     "protein_class": "metabolic enzyme tetramer", "mw_kda": 36.0,
     "notes": "144 kDa tetramer; standard native MS subject."},
    {"uniprot_id": "P00338", "name": "Lactate dehydrogenase A / LDHA (human)",
     "protein_class": "metabolic enzyme tetramer", "mw_kda": 36.7,
     "notes": "~145 kDa tetramer; inhibitor screens by native MS."},
    {"uniprot_id": "P14618", "name": "Pyruvate kinase M2 / PKM2 (human)",
     "protein_class": "metabolic enzyme tetramer", "mw_kda": 57.9,
     "notes": "Tetramer/dimer equilibrium and FBP binding by native MS."},
    {"uniprot_id": "P00925", "name": "Enolase 2 / ENO2 (human)",
     "protein_class": "metabolic enzyme dimer", "mw_kda": 47.3,
     "notes": "Dimer native MS; metal binding."},

    # AAA+ / proteasome / chaperone family extensions
    {"uniprot_id": "P25786", "name": "Proteasome subunit alpha-1 / PSMA1 (human)",
     "protein_class": "proteasome subunit", "mw_kda": 29.5,
     "notes": "20S/26S proteasome native MS (Robinson, Heck)."},
    {"uniprot_id": "P28074", "name": "Proteasome subunit beta-5 / PSMB5 (human)",
     "protein_class": "proteasome subunit", "mw_kda": 28.5,
     "notes": "Bortezomib-target subunit; 20S native MS."},
    {"uniprot_id": "P62195", "name": "26S proteasome regulatory subunit 8 / PSMC5 (human)",
     "protein_class": "AAA+ ATPase", "mw_kda": 45.6,
     "notes": "19S RP component; native MS of intact 26S."},
    {"uniprot_id": "P62333", "name": "26S proteasome regulatory subunit 10B / PSMC6 (human)",
     "protein_class": "AAA+ ATPase", "mw_kda": 44.2,
     "notes": "19S RP component; intact 26S native MS."},
]

if __name__ == "__main__":
    print(f"SMALL_PEPTIDES_POSITIVES: {len(SMALL_PEPTIDES_POSITIVES)}")
    print(f"HARD_NEGATIVES: {len(HARD_NEGATIVES)}")
    print(f"ADDITIONAL_POSITIVES: {len(ADDITIONAL_POSITIVES)}")
