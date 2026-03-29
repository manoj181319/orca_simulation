"""
extract_targets.py  —  Phase 0, Step 2
Reads raw JSON files from pubmed_fetch.py and extracts gene names, protein names,
and ligand mentions relevant to P. aeruginosa redox biology and AMR.

Architecture — Biochemistry-First Ligand Assignment
----------------------------------------------------
Previous versions assigned top_ligand by co-mention frequency from the corpus.
This is fundamentally unreliable because ligands like pyocyanin and iron appear
in almost every P. aeruginosa redox paper regardless of which protein is being
discussed. Co-mention frequency measures biological context overlap, not binding.

Root cause confirmed: pyocyanin is a redox-active virulence factor, QS signal,
and electron shuttle all at once — it appears in papers about catalases, SODs,
efflux pumps, thioredoxins, and QS regulators not because it binds them, but
because it is always present in the same biological context. The same applies
to iron and phenazine.

The correct architecture has three layers:
  1. GENE_BIOCHEMISTRY — authoritative source for each gene's simulation ligand,
     derived from PDB structures and primary literature. This is the ONLY source
     of top_ligand. Adding a new gene requires adding it here.
  2. LIGAND_NOISE_FILTER — globally prevalent ligands (pyocyanin, iron, phenazine,
     pyoverdine, pyochelin) that carry no gene-specific binding signal. Still
     recorded in co-mention counts for transparency but excluded from ligand
     diversity scoring in score_targets.py.
  3. Co-mention counts — collected for informational purposes and the diversity
     bonus score only. Never used to assign top_ligand.

Input:
    literature_mining/raw_papers/pubmed_raw.json
    literature_mining/raw_papers/epmc_raw.json

Output:
    literature_mining/extracted/extracted_targets.json

Usage:
    python extract_targets.py
"""

import json
import re
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_DIR  = Path(__file__).parent / "raw_papers"
OUTPUT_DIR = Path(__file__).parent / "extracted"
PUBMED_RAW  = INPUT_DIR  / "pubmed_raw.json"
EPMC_RAW    = INPUT_DIR  / "epmc_raw.json"
OUTPUT_FILE = OUTPUT_DIR / "extracted_targets.json"

# ---------------------------------------------------------------------------
# Known-gene list — identity and classification only.
# Ligand assignment is in GENE_BIOCHEMISTRY below, not here.
# ---------------------------------------------------------------------------
KNOWN_GENES = {
    # Catalases and peroxidases
    "katA":  {"protein": "Catalase KatA",                        "function": "H2O2 decomposition (primary catalase)",            "category": "oxidative_stress"},
    "katB":  {"protein": "Catalase KatB",                        "function": "H2O2 decomposition (inducible catalase)",          "category": "oxidative_stress"},
    "katE":  {"protein": "Catalase KatE",                        "function": "H2O2 decomposition (stationary phase)",            "category": "oxidative_stress"},
    "ahpC":  {"protein": "Alkyl hydroperoxide reductase C",       "function": "Organic peroxide reduction (2-Cys peroxiredoxin)", "category": "oxidative_stress"},
    "ahpF":  {"protein": "Alkyl hydroperoxide reductase F",       "function": "NADH/FAD-dependent AhpC reductase subunit",       "category": "oxidative_stress"},
    "oxyR":  {"protein": "OxyR",                                  "function": "H2O2-sensing transcription factor",               "category": "regulation"},
    "soxR":  {"protein": "SoxR",                                  "function": "Superoxide/redox-sensing transcription factor",    "category": "regulation"},
    "soxS":  {"protein": "SoxS",                                  "function": "SoxR-regulated AMR gene activator",               "category": "regulation"},
    # Superoxide dismutases
    "sodA":  {"protein": "Superoxide dismutase A (MnSOD)",        "function": "O2- dismutation (manganese-cofactored)",          "category": "oxidative_stress"},
    "sodB":  {"protein": "Superoxide dismutase B (FeSOD)",        "function": "O2- dismutation (iron-cofactored)",               "category": "oxidative_stress"},
    "sodM":  {"protein": "Superoxide dismutase M (cambialistic)", "function": "O2- dismutation",                                 "category": "oxidative_stress"},
    # Thioredoxin system
    "trxA":  {"protein": "Thioredoxin 1",                         "function": "Thiol-disulfide oxidoreductase",                  "category": "redox_enzyme"},
    "trxB":  {"protein": "Thioredoxin reductase",                 "function": "NADPH-dependent Trx regeneration",                "category": "redox_enzyme"},
    "trxC":  {"protein": "Thioredoxin 2 (periplasmic)",           "function": "Thiol-disulfide oxidoreductase",                  "category": "redox_enzyme"},
    # Glutaredoxin / glutathione system
    "grxA":  {"protein": "Glutaredoxin 1",                        "function": "Glutathione-dependent oxidoreductase",            "category": "redox_enzyme"},
    "grxB":  {"protein": "Glutaredoxin 2",                        "function": "Glutathione-dependent oxidoreductase",            "category": "redox_enzyme"},
    "gshA":  {"protein": "Glutamate-cysteine ligase",              "function": "Glutathione biosynthesis step 1",                "category": "redox_enzyme"},
    "gshB":  {"protein": "Glutathione synthetase",                 "function": "Glutathione biosynthesis step 2",                "category": "redox_enzyme"},
    "gor":   {"protein": "Glutathione reductase",                  "function": "NADPH-dependent GSH regeneration",               "category": "redox_enzyme"},
    # NADH/NADPH metabolism and electron transport
    "ndh":   {"protein": "NADH dehydrogenase II",                  "function": "Non-proton-pumping NADH oxidase",                "category": "electron_transport"},
    "nuoA":  {"protein": "NADH:ubiquinone oxidoreductase subunit A","function": "Complex I — NADH oxidase (proton-pumping)",      "category": "electron_transport"},
    "nuoB":  {"protein": "NADH:ubiquinone oxidoreductase subunit B","function": "Complex I — electron input module",              "category": "electron_transport"},
    "nuoN":  {"protein": "NADH:ubiquinone oxidoreductase subunit N","function": "Complex I — proton translocation module",        "category": "electron_transport"},
    "sdhA":  {"protein": "Succinate dehydrogenase subunit A",       "function": "Complex II — succinate oxidation, FAD-binding", "category": "electron_transport"},
    "sdhB":  {"protein": "Succinate dehydrogenase subunit B",       "function": "Complex II — Fe-S cluster electron transfer",   "category": "electron_transport"},
    # Cytochromes and terminal oxidases
    "cioA":  {"protein": "Cyanide-insensitive oxidase subunit A",  "function": "Terminal oxidase (antibiotic tolerance)",         "category": "electron_transport"},
    "cioB":  {"protein": "Cyanide-insensitive oxidase subunit B",  "function": "Terminal oxidase",                                "category": "electron_transport"},
    "cyoA":  {"protein": "Cytochrome bo3 oxidase subunit I",        "function": "Terminal oxidase (aerobic)",                    "category": "electron_transport"},
    "coxA":  {"protein": "Cytochrome c oxidase subunit I",          "function": "Terminal oxidase (microaerobic)",               "category": "electron_transport"},
    # Ferredoxins
    "fdxA":  {"protein": "Ferredoxin I",                            "function": "Electron carrier ([2Fe-2S])",                   "category": "redox_enzyme"},
    "fdxB":  {"protein": "Ferredoxin II",                           "function": "Electron carrier ([2Fe-2S])",                   "category": "redox_enzyme"},
    # Iron metabolism
    "fur":   {"protein": "Ferric uptake regulator Fur",             "function": "Fe2+-sensing transcription factor",             "category": "iron_metabolism"},
    "pvdA":  {"protein": "Pyoverdine biosynthesis protein A",       "function": "L-ornithine hydroxylase (FAD+NADPH-dependent)", "category": "iron_metabolism"},
    "pvdD":  {"protein": "Pyoverdine synthetase D",                 "function": "Non-ribosomal peptide synthetase",              "category": "iron_metabolism"},
    "pchA":  {"protein": "Salicylate biosynthesis protein",         "function": "Pyochelin siderophore biosynthesis",            "category": "iron_metabolism"},
    "bfrA":  {"protein": "Bacterioferritin A",                      "function": "Iron storage (heme b-containing)",              "category": "iron_metabolism"},
    "bfrB":  {"protein": "Bacterioferritin B",                      "function": "Iron storage (heme b-containing)",              "category": "iron_metabolism"},
    # Phenazine biosynthesis
    "phzA1": {"protein": "Phenazine biosynthesis protein A1",       "function": "Phenazine core biosynthesis",                   "category": "phenazine"},
    "phzB1": {"protein": "Phenazine biosynthesis protein B1",       "function": "Phenazine core biosynthesis",                   "category": "phenazine"},
    "phzM":  {"protein": "Phenazine methyltransferase",             "function": "SAM-dependent methylation to 5-Me-PCA",        "category": "phenazine"},
    "phzS":  {"protein": "Phenazine hydroxylase",                   "function": "FAD+NADH-dependent hydroxylation to pyocyanin", "category": "phenazine"},
    "phzH":  {"protein": "Phenazine amide synthetase",              "function": "Phenazine-1-carboxamide biosynthesis",          "category": "phenazine"},
    # Efflux pumps
    "mexR":  {"protein": "MexR",                                    "function": "Redox-sensing repressor of MexAB-OprM efflux",  "category": "regulation"},
    "mexA":  {"protein": "MexA (MexAB-OprM efflux)",               "function": "Drug efflux membrane fusion protein",           "category": "efflux"},
    "mexB":  {"protein": "MexB (MexAB-OprM efflux)",               "function": "RND pump — multidrug efflux",                   "category": "efflux"},
    "oprM":  {"protein": "OprM (MexAB-OprM efflux)",               "function": "Outer membrane channel (passive)",              "category": "efflux"},
    "mexC":  {"protein": "MexC (MexCD-OprJ efflux)",               "function": "RND pump membrane fusion",                     "category": "efflux"},
    "mexD":  {"protein": "MexD (MexCD-OprJ efflux)",               "function": "RND pump — fluoroquinolone efflux",             "category": "efflux"},
    "mexE":  {"protein": "MexE (MexEF-OprN efflux)",               "function": "RND pump membrane fusion",                     "category": "efflux"},
    "mexF":  {"protein": "MexF (MexEF-OprN efflux)",               "function": "RND pump — redox-regulated efflux",             "category": "efflux"},
    "mexX":  {"protein": "MexX (MexXY-OprM efflux)",               "function": "RND pump — aminoglycoside efflux",              "category": "efflux"},
    "mexY":  {"protein": "MexY (MexXY-OprM efflux)",               "function": "RND pump",                                     "category": "efflux"},
    # Quorum sensing
    "lasR":  {"protein": "LasR",                                    "function": "3OC12-HSL-sensing QS receptor/regulator",      "category": "quorum_sensing"},
    "lasI":  {"protein": "LasI",                                    "function": "3OC12-HSL autoinducer synthase",               "category": "quorum_sensing"},
    "rhlR":  {"protein": "RhlR",                                    "function": "C4-HSL-sensing QS receptor/regulator",         "category": "quorum_sensing"},
    "rhlI":  {"protein": "RhlI",                                    "function": "C4-HSL autoinducer synthase",                  "category": "quorum_sensing"},
    "mvfR":  {"protein": "MvfR (PqsR)",                            "function": "PQS/HHQ-binding QS receptor",                  "category": "quorum_sensing"},
    "pqsA":  {"protein": "PqsA",                                    "function": "PQS biosynthesis — anthraniloyl-CoA ligase",   "category": "quorum_sensing"},
}

# ---------------------------------------------------------------------------
# GENE_BIOCHEMISTRY — the single authoritative source for simulation ligands.
#
# top_ligand in all downstream outputs comes ONLY from here.
# Co-mention statistics from the corpus are NEVER used for ligand assignment.
#
# Fields:
#   simulation_ligand : small molecule for ORCA simulation
#   ligand_role       : "cofactor"|"substrate"|"effector"|"product"|"no_ligand"
#   basis             : evidence (PDB ID, UniProt, primary literature reference)
#   simulatable       : True = include in Recommended Pairs; False = exclude
#   note              : important qualifications
#
# "no_ligand" = protein has no small-molecule binding partner suitable for ORCA
# (passive channels, structural subunits, multi-domain enzymes). These are
# automatically excluded from the Recommended Pairs sheet.
# ---------------------------------------------------------------------------
GENE_BIOCHEMISTRY = {
    # Catalases — all heme b cofactor (universal class property)
    "katA":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "PDB 4E37 — Pa PAO1 KatA with HEM+NADPH confirmed",                    "simulatable": True,  "note": "NADPH also present at peripheral site in 4E37"},
    "katB":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "PMC177506 — KatB confirmed as heme catalase in Pa",                    "simulatable": True,  "note": "No Pa-specific PDB; use KatA as structural template"},
    "katE":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "Monofunctional catalase class — universally heme b",                   "simulatable": True,  "note": "Stationary-phase isoform"},
    # Peroxiredoxin system
    "ahpC":  {"simulation_ligand": "organic hydroperoxide","ligand_role": "substrate","basis": "ScienceDirect S2213231721002342 — AhpC neutralizes organic hydroperoxides","simulatable": True,"note": "2-Cys peroxiredoxin; Cys46/Cys165 thiol mechanism; not iron"},
    "ahpF":  {"simulation_ligand": "FAD",                "ligand_role": "cofactor",  "basis": "AhpF is a flavoprotein; FAD+NADH electron donor for AhpC recharging",  "simulatable": True,  "note": "AhpF directly binds FAD; it does NOT bind the hydroperoxide substrate"},
    # Redox-sensing regulators
    "oxyR":  {"simulation_ligand": "H2O2",               "ligand_role": "effector",  "basis": "ScienceDirect S0923250811001756 — OxyR = primary H2O2 sensor in Pa",   "simulatable": True,  "note": "H2O2 oxidises Cys199 to activate OxyR; no stable co-crystal in PDB"},
    "soxR":  {"simulation_ligand": "[2Fe-2S]",           "ligand_role": "cofactor",  "basis": "SoxR contains redox-active [2Fe-2S] cluster as sensing cofactor",      "simulatable": True,  "note": "No Pa-specific PDB; E. coli SoxR structures only (not valid for Pa)"},
    "soxS":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "SoxS is a transcription factor with no small-molecule binding pocket",  "simulatable": False, "note": "Exclude from recommended pairs"},
    "mexR":  {"simulation_ligand": "H2O2",               "ligand_role": "effector",  "basis": "OUP 8222440 — MexR inactivated by H2O2 via Cys30-Cys62 disulphide",   "simulatable": True,  "note": "PDB 1LNW (oxidised), 4LFH — Pa MexR confirmed; pyocyanin binds MexR indirectly"},
    # Superoxide dismutases — metal ion at active site IS the redox cofactor
    "sodA":  {"simulation_ligand": "Mn2+",               "ligand_role": "cofactor",  "basis": "MnSOD class — manganese ion is the catalytic metal at active site",     "simulatable": True,  "note": "NOT pyocyanin; SODs dismutate superoxide using metal cofactor"},
    "sodB":  {"simulation_ligand": "Fe2+",               "ligand_role": "cofactor",  "basis": "PMC177481 — sodB encodes FeSOD in Pa; iron is catalytic cofactor",     "simulatable": True,  "note": "NOT pyocyanin; SodB needed for PYO biosynthesis but does not bind it"},
    "sodM":  {"simulation_ligand": "Mn2+",               "ligand_role": "cofactor",  "basis": "Cambialistic SOD; preferentially uses Mn2+ in Pa",                     "simulatable": True,  "note": "NOT pyocyanin"},
    # Thioredoxin system
    "trxA":  {"simulation_ligand": "NADPH",              "ligand_role": "cofactor",  "basis": "TrxA reduced by TrxB/NADPH system; NADPH is electron source",          "simulatable": True,  "note": "TrxA has no small-molecule site; NADPH acts through TrxB; NOT glutathione"},
    "trxB":  {"simulation_ligand": "FAD",                "ligand_role": "cofactor",  "basis": "Wikipedia Thioredoxin_reductase — each monomer contains FAD cofactor", "simulatable": True,  "note": "FAD is the redox prosthetic group; NADPH is electron donor"},
    "trxC":  {"simulation_ligand": "NADPH",              "ligand_role": "cofactor",  "basis": "Periplasmic thioredoxin; same NADPH-dependent system as TrxA",         "simulatable": True,  "note": ""},
    # Glutaredoxin / glutathione system
    "grxA":  {"simulation_ligand": "glutathione",        "ligand_role": "cofactor",  "basis": "Glutaredoxin class — GSH is the defining reductant cofactor",          "simulatable": True,  "note": ""},
    "grxB":  {"simulation_ligand": "glutathione",        "ligand_role": "cofactor",  "basis": "Glutaredoxin class — GSH is the defining reductant cofactor",          "simulatable": True,  "note": ""},
    "gshA":  {"simulation_ligand": "glutathione",        "ligand_role": "product",   "basis": "GshA catalyses step 1 of GSH biosynthesis",                           "simulatable": True,  "note": "GSH is the biologically relevant molecule for simulation"},
    "gshB":  {"simulation_ligand": "glutathione",        "ligand_role": "product",   "basis": "GshB catalyses step 2 of GSH biosynthesis",                           "simulatable": True,  "note": ""},
    "gor":   {"simulation_ligand": "FAD",                "ligand_role": "cofactor",  "basis": "Glutathione reductase is a flavoprotein with FAD cofactor",            "simulatable": True,  "note": "GSH is also simulatable as substrate"},
    # NADH/NADPH metabolism
    "ndh":   {"simulation_ligand": "NADH",               "ligand_role": "substrate", "basis": "NADH dehydrogenase II oxidises NADH — class definition",              "simulatable": True,  "note": ""},
    "nuoA":  {"simulation_ligand": "NADH",               "ligand_role": "substrate", "basis": "Complex I catalytic core oxidises NADH",                              "simulatable": True,  "note": ""},
    "nuoB":  {"simulation_ligand": "NADH",               "ligand_role": "substrate", "basis": "Complex I NADH-input subunit",                                        "simulatable": True,  "note": ""},
    "nuoN":  {"simulation_ligand": "ubiquinone",         "ligand_role": "substrate", "basis": "Complex I proton-pumping module; ubiquinone is terminal electron acceptor","simulatable": True,"note": ""},
    "sdhA":  {"simulation_ligand": "FAD",                "ligand_role": "cofactor",  "basis": "SdhA is the FAD-containing catalytic subunit of complex II",           "simulatable": True,  "note": "Ubiquinone also present as electron acceptor"},
    "sdhB":  {"simulation_ligand": "ubiquinone",         "ligand_role": "substrate", "basis": "SdhB Fe-S cluster transfers electrons to ubiquinone",                 "simulatable": True,  "note": "NOT pyocyanin"},
    # Cytochromes and terminal oxidases
    "cioA":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "Cyanide-insensitive oxidase contains heme b cofactor",                "simulatable": True,  "note": ""},
    "cioB":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "CioB subunit contains heme cofactor",                                 "simulatable": True,  "note": ""},
    "cyoA":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "Cytochrome bo3 contains heme b and heme o",                           "simulatable": True,  "note": ""},
    "coxA":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "Cytochrome c oxidase contains heme a and Cu centres",                 "simulatable": True,  "note": ""},
    # Ferredoxins
    "fdxA":  {"simulation_ligand": "ferredoxin",         "ligand_role": "cofactor",  "basis": "Ferredoxin I contains [2Fe-2S] cluster as redox-active centre",       "simulatable": True,  "note": ""},
    "fdxB":  {"simulation_ligand": "ferredoxin",         "ligand_role": "cofactor",  "basis": "Ferredoxin II contains [2Fe-2S] cluster",                             "simulatable": True,  "note": ""},
    # Iron metabolism
    "fur":   {"simulation_ligand": "Fe2+",               "ligand_role": "effector",  "basis": "PMC5648857 — Fur uses Fe2+ as corepressor; Zn2+ is structural only",  "simulatable": True,  "note": "NOT heme; Fur is an iron homeostasis regulator, not a heme sensor"},
    "pvdA":  {"simulation_ligand": "FAD",                "ligand_role": "cofactor",  "basis": "PvdA is L-ornithine hydroxylase requiring FAD+NADPH+O2",              "simulatable": True,  "note": "FAD is the redox cofactor; 'flavin' was too vague in previous version"},
    "pvdD":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "NRPS megaenzyme; no single small-molecule ligand for ORCA simulation", "simulatable": False, "note": "Exclude from recommended pairs"},
    "pchA":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "Isochorismate synthase; no redox cofactor",                           "simulatable": False, "note": "Exclude from recommended pairs"},
    "bfrA":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "Bacterioferritin contains heme b as electron-transfer cofactor",      "simulatable": True,  "note": ""},
    "bfrB":  {"simulation_ligand": "heme",               "ligand_role": "cofactor",  "basis": "Bacterioferritin contains heme b as electron-transfer cofactor",      "simulatable": True,  "note": ""},
    # Phenazine biosynthesis — pyocyanin is their PRODUCT, not their ligand
    "phzA1": {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "PhzA1 is an isomerase in core biosynthesis; no redox cofactor",       "simulatable": False, "note": "Pyocyanin is downstream product, not a binding ligand"},
    "phzB1": {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "PhzB1 is an isomerase; no redox cofactor",                            "simulatable": False, "note": ""},
    "phzM":  {"simulation_ligand": "SAM",                "ligand_role": "substrate", "basis": "PhzM is SAM-dependent methyltransferase; SAM is methyl donor",        "simulatable": True,  "note": "Pyocyanin is the product; SAM is the substrate. NOT pyocyanin"},
    "phzS":  {"simulation_ligand": "FAD",                "ligand_role": "cofactor",  "basis": "PhzS is flavin-dependent monooxygenase; FAD+NADH are cofactors",      "simulatable": True,  "note": "Pyocyanin is the product; FAD is the cofactor. NOT pyocyanin"},
    "phzH":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "PhzH is glutamine amidotransferase; no redox cofactor",               "simulatable": False, "note": ""},
    # Efflux pumps — passive structural components excluded from simulation
    "mexA":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "MexA is membrane fusion protein; no ligand-binding pocket",           "simulatable": False, "note": ""},
    "mexB":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "MexB pump translocase; drug substrates not ORCA targets",             "simulatable": False, "note": ""},
    "oprM":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "OprM is passive outer membrane channel; no small-molecule binding",   "simulatable": False, "note": "Pyocyanin binds MexR/NalD upstream, NOT OprM itself"},
    "mexC":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "Membrane fusion protein",                                             "simulatable": False, "note": ""},
    "mexD":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "Pump translocase",                                                    "simulatable": False, "note": ""},
    "mexE":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "Membrane fusion protein",                                             "simulatable": False, "note": ""},
    "mexF":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "Pump translocase",                                                    "simulatable": False, "note": ""},
    "mexX":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "Pump translocase",                                                    "simulatable": False, "note": ""},
    "mexY":  {"simulation_ligand": "no_ligand",          "ligand_role": "no_ligand", "basis": "Pump translocase",                                                    "simulatable": False, "note": ""},
    # Quorum sensing
    "lasR":  {"simulation_ligand": "3OC12-HSL",          "ligand_role": "effector",  "basis": "PDB 2UV0, 3IX3 — LasR LBD with 3OC12-HSL co-crystal confirmed",      "simulatable": True,  "note": "Co-crystal in PDB; no docking needed"},
    "lasI":  {"simulation_ligand": "3OC12-HSL",          "ligand_role": "product",   "basis": "PubMed 15306017 — LasI synthesises 3OC12-HSL as product",             "simulatable": True,  "note": "Product is relevant for simulation context"},
    "rhlR":  {"simulation_ligand": "C4-HSL",             "ligand_role": "effector",  "basis": "PDB 3T5K, 7R3H, 8B4A — RhlR with C4-HSL co-crystal confirmed",       "simulatable": True,  "note": "Co-crystal in PDB; no docking needed"},
    "rhlI":  {"simulation_ligand": "C4-HSL",             "ligand_role": "product",   "basis": "Nature Comm s41467-023-43702-4 — RhlI synthesises C4-HSL",            "simulatable": True,  "note": ""},
    "mvfR":  {"simulation_ligand": "HHQ",                "ligand_role": "effector",  "basis": "Nature Comm s41467-023-43702-4 — PqsR binds HHQ and PQS, NOT pyocyanin","simulatable": True,"note": "HHQ = 2-heptyl-4-hydroxyquinoline; PQS = Pseudomonas quinolone signal"},
    "pqsA":  {"simulation_ligand": "CoA",                "ligand_role": "substrate", "basis": "PqsA is anthraniloyl-CoA ligase; CoA is the substrate",               "simulatable": True,  "note": ""},
}

# ---------------------------------------------------------------------------
# LIGAND_NOISE — globally prevalent ligands in Pa literature that appear in
# nearly every paper regardless of protein being studied.
# Root cause of pyocyanin contamination: it is a redox-active virulence factor
# and QS signal simultaneously, so it co-occurs with every gene in the dataset.
# Similarly, iron is discussed in nearly all Pa metabolism papers.
#
# These ligands are still recorded in ligands_noise_co_mentioned for
# transparency, but are separated from the signal co-mention counts and
# excluded from ligand diversity scoring in score_targets.py.
# ---------------------------------------------------------------------------
LIGAND_NOISE = {
    "pyocyanin",  # redox-active virulence factor; mentioned in >90% of Pa redox papers
    "phenazine",  # parent class of pyocyanin; same prevalence problem
    "iron",       # ubiquitous in Pa biology; not gene-specific
    "pyoverdine", # siderophore; mentioned broadly in AMR/virulence context
    "pyochelin",  # siderophore; same issue
}

# ---------------------------------------------------------------------------
# Ligands to detect in text (for co-mention counting only)
# ---------------------------------------------------------------------------
KNOWN_LIGANDS = [
    "NADH", "NADPH", "FAD", "FMN", "heme", "haem", "ubiquinone", "menaquinone",
    "pyocyanin", "pyoverdine", "pyochelin", "phenazine", "glutathione",
    "ferredoxin", "coenzyme A", "CoA", "ATP", "GTP", "iron", "Fe2+", "Fe3+",
    "zinc", "copper", "molybdenum", "flavin", "quinone", "SAM",
    "S-adenosylmethionine", "Mn2+", "HHQ", "PQS",
    "3OC12-HSL", "3-oxo-C12-HSL", "C4-HSL", "N-butanoyl-HSL",
    "H2O2", "hydrogen peroxide", "superoxide", "organic hydroperoxide",
]

# Synonym normalisation — applied after extraction so the JSON always uses
# the canonical name. This prevents duplicate entries in the Ligand Mentions
# sheet (e.g. "H2O2" and "hydrogen peroxide" appearing as separate rows).
LIGAND_SYNONYMS = {
    "hydrogen peroxide":   "H2O2",
    "S-adenosylmethionine":"SAM",
    "haem":                "heme",
    "N-butanoyl-HSL":      "C4-HSL",
    "3-oxo-C12-HSL":       "3OC12-HSL",
    "coenzyme A":          "CoA",
}

LOCUS_TAG_PATTERN = re.compile(r'\bPA\d{4}\b')


# ---------------------------------------------------------------------------
# Text extraction helpers
# ---------------------------------------------------------------------------
def get_text(paper: dict) -> str:
    return f"{paper.get('title', '')} {paper.get('abstract', '')}".lower()


def extract_genes_from_text(text: str) -> set:
    found = set()
    for gene in KNOWN_GENES:
        if re.search(rf'\b{re.escape(gene)}\b', text, re.IGNORECASE):
            found.add(gene)
    return found


def extract_ligands_from_text(text: str) -> set:
    found = set()
    for ligand in KNOWN_LIGANDS:
        if re.search(rf'\b{re.escape(ligand)}\b', text, re.IGNORECASE):
            # Normalise to canonical name before storing
            canonical = LIGAND_SYNONYMS.get(ligand.lower(), LIGAND_SYNONYMS.get(ligand, ligand))
            found.add(canonical)
    return found


def extract_novel_locus_tags(text: str) -> set:
    return set(LOCUS_TAG_PATTERN.findall(text.upper()))


# ---------------------------------------------------------------------------
# Paper loading / deduplication
# ---------------------------------------------------------------------------
def load_papers(*json_paths) -> list:
    papers = []
    for path in json_paths:
        if not path.exists():
            print(f"WARNING: {path} not found — skipping.")
            continue
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        papers.extend(data)
        print(f"  Loaded {len(data)} papers from {path.name}")
    return papers


def deduplicate_papers(papers: list) -> list:
    seen_pmids, seen_dois = set(), set()
    unique = []
    for p in papers:
        pmid = p.get("pmid", "").strip()
        doi  = p.get("doi",  "").strip().lower()
        if pmid and pmid in seen_pmids:
            continue
        if doi and doi in seen_dois:
            continue
        unique.append(p)
        if pmid: seen_pmids.add(pmid)
        if doi:  seen_dois.add(doi)
    return unique


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading raw paper files...")
    papers = load_papers(PUBMED_RAW, EPMC_RAW)
    papers = deduplicate_papers(papers)
    print(f"Total unique papers for extraction: {len(papers):,}\n")

    if len(papers) < 1000:
        print(f"WARNING: Only {len(papers)} papers. Run pubmed_fetch.py first.")

    gene_data   = defaultdict(lambda: {
        "count": 0, "papers": [],
        "ligands_co_mentioned": defaultdict(int),
        "ligands_noise_co_mentioned": defaultdict(int),
    })
    ligand_data = defaultdict(lambda: {"count": 0, "papers": []})
    novel_loci  = defaultdict(lambda: {"count": 0, "papers": []})

    for paper in papers:
        text      = get_text(paper)
        paper_ref = {"title": paper.get("title", "No title"), "pmid": paper.get("pmid", ""), "year": paper.get("year", "")}

        genes_found   = extract_genes_from_text(text)
        ligands_found = extract_ligands_from_text(text)
        loci_found    = extract_novel_locus_tags(text)

        for gene in genes_found:
            gd = gene_data[gene]
            gd["count"] += 1
            if len(gd["papers"]) < 5:
                gd["papers"].append(paper_ref)
            for ligand in ligands_found:
                if ligand in LIGAND_NOISE:
                    gd["ligands_noise_co_mentioned"][ligand] += 1
                else:
                    gd["ligands_co_mentioned"][ligand] += 1

        for ligand in ligands_found:
            ligand_data[ligand]["count"] += 1
            if len(ligand_data[ligand]["papers"]) < 3:
                ligand_data[ligand]["papers"].append(paper_ref)

        for locus in loci_found:
            novel_loci[locus]["count"] += 1
            if len(novel_loci[locus]["papers"]) < 3:
                novel_loci[locus]["papers"].append(paper_ref)

    # Build gene output — top_ligand from GENE_BIOCHEMISTRY ONLY
    gene_output    = {}
    missing_biochem = []

    for gene, data in sorted(gene_data.items(), key=lambda x: -x[1]["count"]):
        biochem = GENE_BIOCHEMISTRY.get(gene)
        if not biochem:
            missing_biochem.append(gene)
            biochem = {
                "simulation_ligand": "unknown",
                "ligand_role": "unknown",
                "basis": "MISSING from GENE_BIOCHEMISTRY — must add entry",
                "simulatable": False,
                "note": "",
            }
        gene_info = KNOWN_GENES.get(gene, {})
        gene_output[gene] = {
            "protein":      gene_info.get("protein",  ""),
            "function":     gene_info.get("function", ""),
            "category":     gene_info.get("category", ""),
            "paper_count":  data["count"],
            "sample_papers": data["papers"],
            # Authoritative biochemistry-derived fields
            "top_ligand":   biochem["simulation_ligand"],
            "ligand_role":  biochem["ligand_role"],
            "ligand_basis": biochem["basis"],
            "simulatable":  biochem["simulatable"],
            "ligand_note":  biochem.get("note", ""),
            # Co-mention data — informational only, never drives ligand assignment
            "ligands_co_mentioned":       dict(sorted(data["ligands_co_mentioned"].items(),       key=lambda x: -x[1])),
            "ligands_noise_co_mentioned": dict(sorted(data["ligands_noise_co_mentioned"].items(), key=lambda x: -x[1])),
        }

    if missing_biochem:
        print(f"\nACTION REQUIRED: {len(missing_biochem)} genes need GENE_BIOCHEMISTRY entries:")
        for g in missing_biochem:
            print(f"  - {g}")

    ligand_output = {
        k: {"count": v["count"], "is_noise": k in LIGAND_NOISE, "papers": v["papers"]}
        for k, v in sorted(ligand_data.items(), key=lambda x: -x[1]["count"])
    }
    novel_output = {
        k: v for k, v in sorted(novel_loci.items(), key=lambda x: -x[1]["count"])
        if k not in {g.upper() for g in KNOWN_GENES}
    }

    result = {
        "total_papers_processed": len(papers),
        "ligand_noise_filter": sorted(LIGAND_NOISE),
        "genes": gene_output,
        "ligands": ligand_output,
        "novel_locus_tags": novel_output,
    }

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print("\nExtraction complete.")
    print(f"  Genes detected:         {len(gene_output)}")
    print(f"  Ligands detected:       {len(ligand_output)}")
    print(f"  Novel locus tags found: {len(novel_output)}")
    simulatable_count = sum(1 for g in gene_output.values() if g["simulatable"])
    print(f"  Simulatable pairs:      {simulatable_count}")
    print(f"\nTop 10 genes by paper count:")
    for i, (gene, data) in enumerate(list(gene_output.items())[:10], 1):
        sim = "simulatable" if data["simulatable"] else "excluded"
        print(f"  {i:2}. {gene:<8} — {data['paper_count']:>4} papers — ligand: {data['top_ligand']:<22} ({sim})")
    print(f"\nSaved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
