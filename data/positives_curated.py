"""
Curated list of REAL proteins documented as successful native mass spec
targets in published literature. Each entry corresponds to a real protein
with a verifiable UniProt accession.

DOIs are best-effort from literature knowledge and should be verified
before publication, but the proteins themselves are real and the native MS
literature on each is well-documented.

Errors flagged by the research agent have been corrected (streptavidin
accession fixed, GDH accession fixed, duplicates removed).
"""

POSITIVE_EXAMPLES = [
    # Small soluble model proteins
    {"uniprot_id": "P00918", "name": "Carbonic anhydrase 2 (human)", "protein_class": "enzyme", "mw_kda": 29.2, "notes": "Canonical native MS calibrant"},
    {"uniprot_id": "P0CG48", "name": "Ubiquitin (human)", "protein_class": "small soluble", "mw_kda": 8.6, "notes": "Standard CSD reference"},
    {"uniprot_id": "P00698", "name": "Lysozyme C (chicken)", "protein_class": "enzyme", "mw_kda": 14.3, "notes": "Disulfide-stabilized standard"},
    {"uniprot_id": "P00004", "name": "Cytochrome c (equine)", "protein_class": "heme protein", "mw_kda": 12.4, "notes": "Foundational ESI-MS folding study"},
    {"uniprot_id": "P02185", "name": "Myoglobin (sperm whale)", "protein_class": "heme protein", "mw_kda": 17.2, "notes": "Heme retention textbook test"},
    {"uniprot_id": "P02144", "name": "Myoglobin (horse heart)", "protein_class": "heme protein", "mw_kda": 17.0, "notes": "Heme-loss benchmark"},
    {"uniprot_id": "P02666", "name": "Beta-casein (bovine)", "protein_class": "intrinsically disordered", "mw_kda": 24.0, "notes": "IDP test case, broad CSD"},

    # Hemoglobin and tetramers
    {"uniprot_id": "P69905", "name": "Hemoglobin alpha (human)", "protein_class": "tetramer", "mw_kda": 15.1, "notes": "Native HbA tetramer demo"},
    {"uniprot_id": "P68871", "name": "Hemoglobin beta (human)", "protein_class": "tetramer", "mw_kda": 15.9, "notes": "Forms HbA tetramer"},
    {"uniprot_id": "P02766", "name": "Transthyretin (human)", "protein_class": "amyloidogenic tetramer", "mw_kda": 13.8, "notes": "55 kDa tetramer, drug binding"},
    {"uniprot_id": "P02769", "name": "Serum albumin (bovine, BSA)", "protein_class": "transport protein", "mw_kda": 66.5, "notes": "17 disulfides; drug binding"},
    {"uniprot_id": "P02768", "name": "Serum albumin (human, HSA)", "protein_class": "transport protein", "mw_kda": 66.5, "notes": "Drug-binding screens"},
    {"uniprot_id": "P00760", "name": "Trypsin (bovine)", "protein_class": "protease", "mw_kda": 23.3, "notes": "Inhibitor binding native MS"},
    {"uniprot_id": "P00974", "name": "BPTI", "protein_class": "protease inhibitor", "mw_kda": 6.5, "notes": "Disulfide-rich benchmark"},

    # Antibodies
    {"uniprot_id": "P01857", "name": "IgG1 heavy chain constant", "protein_class": "antibody", "mw_kda": 50.0, "notes": "Generic IgG1 framework"},
    {"uniprot_id": "P01834", "name": "Ig kappa light chain constant", "protein_class": "antibody", "mw_kda": 25.0, "notes": "Light chain"},
    {"uniprot_id": "P01859", "name": "IgG2 heavy chain constant", "protein_class": "antibody", "mw_kda": 50.0, "notes": "Disulfide isoforms"},
    {"uniprot_id": "P01860", "name": "IgG3 heavy chain constant", "protein_class": "antibody", "mw_kda": 56.0, "notes": "Long hinge, harder to ionize"},
    {"uniprot_id": "P01861", "name": "IgG4 heavy chain constant", "protein_class": "antibody", "mw_kda": 50.0, "notes": "Fab-arm exchange observable"},

    # Membrane proteins (Robinson lab)
    {"uniprot_id": "P0AF00", "name": "AmtB (E. coli)", "protein_class": "membrane protein", "mw_kda": 44.5, "notes": "Landmark membrane native MS"},
    {"uniprot_id": "P02943", "name": "Maltose-binding protein", "protein_class": "periplasmic", "mw_kda": 40.7, "notes": "Sugar-binding monomer"},
    {"uniprot_id": "P02931", "name": "OmpF (E. coli)", "protein_class": "beta-barrel membrane", "mw_kda": 39.3, "notes": "Trimeric porin"},
    {"uniprot_id": "P02932", "name": "OmpA (E. coli)", "protein_class": "beta-barrel membrane", "mw_kda": 35.2, "notes": "Nanodisc native MS"},
    {"uniprot_id": "P0ABB0", "name": "F1-ATPase alpha (E. coli)", "protein_class": "membrane complex", "mw_kda": 55.2, "notes": "Part of intact ATP synthase"},
    {"uniprot_id": "P0ABB4", "name": "F1-ATPase beta (E. coli)", "protein_class": "membrane complex", "mw_kda": 50.3, "notes": "Catalytic subunit"},
    {"uniprot_id": "P11166", "name": "GLUT1 (human)", "protein_class": "membrane transporter", "mw_kda": 54.1, "notes": "MFS-fold, lipid effects"},
    {"uniprot_id": "P30536", "name": "TSPO (human)", "protein_class": "membrane protein", "mw_kda": 18.8, "notes": "Drug binding by native MS"},
    {"uniprot_id": "P0DTC2", "name": "SARS-CoV-2 Spike", "protein_class": "viral glycoprotein", "mw_kda": 141.0, "notes": "Heavily glycosylated trimer"},

    # Large complexes / chaperones
    {"uniprot_id": "P0A6F5", "name": "GroEL (E. coli)", "protein_class": "chaperonin", "mw_kda": 57.2, "notes": "800 kDa tetradecamer"},
    {"uniprot_id": "P0A6F9", "name": "GroES (E. coli)", "protein_class": "co-chaperonin", "mw_kda": 10.4, "notes": "73 kDa heptamer"},
    {"uniprot_id": "P11142", "name": "Hsc70 (human)", "protein_class": "chaperone", "mw_kda": 70.9, "notes": "ATP/ADP conformations"},
    {"uniprot_id": "P07900", "name": "Hsp90 alpha (human)", "protein_class": "chaperone dimer", "mw_kda": 84.6, "notes": "170 kDa dimer"},
    {"uniprot_id": "P08238", "name": "Hsp90 beta (human)", "protein_class": "chaperone dimer", "mw_kda": 83.2, "notes": "Co-chaperone complexes"},
    {"uniprot_id": "P31947", "name": "14-3-3 sigma (human)", "protein_class": "scaffold dimer", "mw_kda": 27.8, "notes": "PPI stabilizer screens"},
    {"uniprot_id": "P63104", "name": "14-3-3 zeta (human)", "protein_class": "scaffold dimer", "mw_kda": 27.7, "notes": "PPI native MS"},
    {"uniprot_id": "P0A6Y8", "name": "DnaK (E. coli Hsp70)", "protein_class": "chaperone", "mw_kda": 69.1, "notes": "Nucleotide-state conformations"},

    # Viral capsids
    {"uniprot_id": "P0C6A0", "name": "HBV core antigen", "protein_class": "viral capsid", "mw_kda": 21.0, "notes": "Multi-MDa capsids"},
    {"uniprot_id": "P03612", "name": "MS2 phage coat", "protein_class": "viral capsid", "mw_kda": 13.7, "notes": "180-mer T=3 capsid"},
    {"uniprot_id": "P03611", "name": "Qbeta phage coat", "protein_class": "viral capsid", "mw_kda": 14.1, "notes": "Cargo-loading studies"},
    {"uniprot_id": "P69470", "name": "AAV2 capsid VP1", "protein_class": "viral capsid", "mw_kda": 87.0, "notes": "Gene therapy QC"},
    {"uniprot_id": "P03435", "name": "Influenza A HA (H1)", "protein_class": "viral glycoprotein", "mw_kda": 63.0, "notes": "Glycoform analysis"},

    # Intrinsically disordered proteins
    {"uniprot_id": "P37840", "name": "Alpha-synuclein (human)", "protein_class": "IDP", "mw_kda": 14.5, "notes": "Oligomer tracking"},
    {"uniprot_id": "P05067", "name": "APP / Abeta (human)", "protein_class": "IDP/amyloid", "mw_kda": 86.9, "notes": "Abeta oligomer distributions"},
    {"uniprot_id": "P10636", "name": "Tau (human)", "protein_class": "IDP", "mw_kda": 45.9, "notes": "Highly disordered, broad CSD"},
    {"uniprot_id": "P04637", "name": "p53 (human)", "protein_class": "transcription factor", "mw_kda": 43.7, "notes": "Tetramer + DNA"},

    # Glycoproteins / clinical biotherapeutics
    {"uniprot_id": "P01133", "name": "EGF (human)", "protein_class": "growth factor", "mw_kda": 6.2, "notes": "Receptor binding"},
    {"uniprot_id": "P00533", "name": "EGFR (human)", "protein_class": "RTK ectodomain", "mw_kda": 134.0, "notes": "Glycosylated, difficult"},
    {"uniprot_id": "P04626", "name": "HER2 / ERBB2 (human)", "protein_class": "RTK", "mw_kda": 138.0, "notes": "Trastuzumab/T-DM1 target"},
    {"uniprot_id": "P02787", "name": "Transferrin (human)", "protein_class": "glycoprotein", "mw_kda": 77.0, "notes": "Glycoform analysis"},
    {"uniprot_id": "P02788", "name": "Lactoferrin (human)", "protein_class": "glycoprotein", "mw_kda": 78.2, "notes": "Glycoform heterogeneity"},
    {"uniprot_id": "P01008", "name": "Antithrombin III", "protein_class": "serpin glycoprotein", "mw_kda": 52.6, "notes": "Heparin binding"},
    {"uniprot_id": "P00734", "name": "Prothrombin", "protein_class": "serine protease", "mw_kda": 70.0, "notes": "Drug binding"},
    {"uniprot_id": "P01112", "name": "HRAS (human)", "protein_class": "small GTPase", "mw_kda": 21.3, "notes": "Nucleotide state"},
    {"uniprot_id": "P01116", "name": "KRAS (human)", "protein_class": "small GTPase", "mw_kda": 21.7, "notes": "G12C inhibitor target"},
    {"uniprot_id": "P0DP23", "name": "Calmodulin (human)", "protein_class": "Ca-binding regulator", "mw_kda": 16.8, "notes": "Ca2+ stoichiometry"},
    {"uniprot_id": "P02753", "name": "RBP4 (human)", "protein_class": "lipocalin", "mw_kda": 23.0, "notes": "Retinoid binding"},
    {"uniprot_id": "P19429", "name": "Cardiac troponin I", "protein_class": "regulatory", "mw_kda": 24.0, "notes": "Troponin complex"},
    {"uniprot_id": "P63316", "name": "Cardiac troponin C", "protein_class": "Ca-binding subunit", "mw_kda": 18.4, "notes": "Ca-loading states"},
    {"uniprot_id": "P45379", "name": "Cardiac troponin T", "protein_class": "structural subunit", "mw_kda": 35.9, "notes": "Largest cTn subunit"},

    # Ribosome / RNP
    {"uniprot_id": "P0A7V0", "name": "30S ribosomal protein S2 (E. coli)", "protein_class": "ribosomal protein", "mw_kda": 26.7, "notes": "30S/70S RNP marker"},
    {"uniprot_id": "P0A7L0", "name": "50S ribosomal protein L1 (E. coli)", "protein_class": "ribosomal protein", "mw_kda": 24.7, "notes": "50S subunit native MS"},

    # SID benchmarks
    {"uniprot_id": "P22629", "name": "Streptavidin", "protein_class": "tetramer", "mw_kda": 13.2, "notes": "53 kDa biotin-binding tetramer"},
    {"uniprot_id": "P02701", "name": "Avidin (chicken)", "protein_class": "glycoprotein tetramer", "mw_kda": 16.4, "notes": "Biotin stoichiometry"},
    {"uniprot_id": "P00722", "name": "Beta-galactosidase (E. coli)", "protein_class": "enzyme tetramer", "mw_kda": 116.5, "notes": "465 kDa tetramer"},
    {"uniprot_id": "P00946", "name": "Concanavalin A", "protein_class": "lectin tetramer", "mw_kda": 25.6, "notes": "pH-dependent assembly"},
    {"uniprot_id": "P00366", "name": "Glutamate dehydrogenase (bovine)", "protein_class": "enzyme hexamer", "mw_kda": 61.4, "notes": "336 kDa hexamer calibrant"},
    {"uniprot_id": "P00563", "name": "Creatine kinase M (rabbit)", "protein_class": "enzyme dimer", "mw_kda": 43.0, "notes": "Nucleotide binding"},
    {"uniprot_id": "P02662", "name": "Alpha-S1-casein (bovine)", "protein_class": "IDP/phosphoprotein", "mw_kda": 24.5, "notes": "Phosphorylated IDP"},

    # Phosphorylase / additional standards
    {"uniprot_id": "P00489", "name": "Glycogen phosphorylase b (rabbit)", "protein_class": "enzyme dimer", "mw_kda": 97.4, "notes": "195 kDa dimer mass standard"},
]

if __name__ == "__main__":
    print(f"Total positive examples: {len(POSITIVE_EXAMPLES)}")
    classes = {}
    for p in POSITIVE_EXAMPLES:
        classes[p["protein_class"]] = classes.get(p["protein_class"], 0) + 1
    print("\nProtein class breakdown:")
    for c, n in sorted(classes.items(), key=lambda x: -x[1]):
        print(f"  {c}: {n}")
