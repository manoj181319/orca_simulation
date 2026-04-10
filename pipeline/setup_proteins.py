"""
setup_proteins.py  —  Phase 1, Pre-Step 0
Automatically creates the proteins/{protein_id}/ folder and config.json
for every simulatable protein found in the Phase 0 literature mining results.

All protein data (PDB IDs, ligand residue names, charge, multiplicity) is
embedded in the PROTEIN_REGISTRY below. Proteins marked with
review_multiplicity=True have metal centers or unusual spin states — these
MUST be confirmed with Shiva Sir before running simulations.

Usage:
    python setup_proteins.py                    # create all proteins
    python setup_proteins.py --dry-run          # preview without writing
    python setup_proteins.py --only katA_heme lasR_3OC12HSL
    python setup_proteins.py --force            # overwrite existing configs
"""

import argparse
import json
import sys
from pathlib import Path

PROTEINS_DIR = Path(__file__).parent.parent / "proteins"

# ─────────────────────────────────────────────────────────────────────────────
# PROTEIN REGISTRY
# Source: GENE_BIOCHEMISTRY table in extract_targets.py (Phase 0)
#         Only proteins with simulatable=True are included.
#         PDB IDs from GENE_REFERENCES in report_targets.py and score_targets.py
#
# pdb_id          : RCSB PDB accession. Empty string = needs manual lookup.
# ligand_residue  : HETATM residue name in the PDB file (3-letter code).
#                   Empty string = needs manual lookup in PDB viewer.
# pdb_notes       : How to find the correct PDB ID / residue name.
#
# Charge/multiplicity defaults:
#   Organic-only QM regions (heme without metal, FAD, NADPH, etc.):
#     neutral  → charge 0,  multiplicity 1  (closed-shell singlet)
#     radical  → charge +1, multiplicity 2  (doublet, one electron removed)
#     reduced  → charge -1, multiplicity 2  (doublet, one electron added)
#   Metal-containing systems (Fe, Mn, Cu centers):
#     marked review_multiplicity=True — confirm ALL values with Shiva Sir.
#     Placeholder multiplicity=1 is set but IS WRONG for high-spin metal centers.
#
# multiplicity_notes: Explains what the correct value should be and why.
# ─────────────────────────────────────────────────────────────────────────────

PROTEIN_REGISTRY = {

    # ── CATALASES (heme b cofactor) ──────────────────────────────────────────
    "katA_heme": {
        "gene": "katA",
        "protein_name": "Catalase KatA",
        "ligand_name": "heme",
        "pdb_id": "4E37",
        "ligand_residue": "HEM",
        "pdb_notes": "PDB 4E37 — Pa PAO1 KatA with heme b + NADPH confirmed.",
        "review_multiplicity": True,
        "multiplicity_notes": (
            "Catalase resting state has Fe(III)-heme. High-spin Fe(III) gives "
            "S=5/2 (mult=6) but the QM region spin depends on ligand field. "
            "Use the spin state appropriate for the oxidation state you want to model. "
            "Confirm with Shiva Sir before running."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "katB_heme": {
        "gene": "katB",
        "protein_name": "Catalase KatB",
        "ligand_name": "heme",
        "pdb_id": "4E37",
        "ligand_residue": "HEM",
        "pdb_notes": (
            "No Pa-specific PDB for KatB. Using KatA (4E37) as structural template "
            "per GENE_BIOCHEMISTRY. Confirm with Shiva Sir if a KatB-specific "
            "structure should be used instead."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": "Same as katA_heme — Fe(III) heme spin state must be confirmed.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "katE_heme": {
        "gene": "katE",
        "protein_name": "Catalase KatE",
        "ligand_name": "heme",
        "pdb_id": "",
        "ligand_residue": "HEM",
        "pdb_notes": (
            "No confirmed Pa-specific PDB for KatE. Search RCSB for "
            "'catalase Pseudomonas aeruginosa stationary phase'. "
            "E. coli HPII catalase (PDB 1IPH) may serve as template."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": "Same Fe(III) heme considerations as katA_heme.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── PEROXIREDOXIN SYSTEM ──────────────────────────────────────────────────
    "ahpC_OHP": {
        "gene": "ahpC",
        "protein_name": "Alkyl hydroperoxide reductase C",
        "ligand_name": "organic hydroperoxide",
        "pdb_id": "",
        "ligand_residue": "",
        "pdb_notes": (
            "AhpC substrate (organic hydroperoxide) is transient — no stable co-crystal. "
            "Search RCSB for Pa AhpC: try 'ahpC Pseudomonas aeruginosa'. "
            "Cumene hydroperoxide or tert-butyl hydroperoxide can be used as model substrate "
            "— confirm choice with Shiva Sir. Ligand must be manually added to PDB."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Standard organic ligand — default values should be correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "ahpF_FAD": {
        "gene": "ahpF",
        "protein_name": "Alkyl hydroperoxide reductase F",
        "ligand_name": "FAD",
        "pdb_id": "",
        "ligand_residue": "FAD",
        "pdb_notes": (
            "Search RCSB for Pa AhpF: try 'alkyl hydroperoxide reductase F "
            "Pseudomonas'. PDB 1KLU is AhpF from Salmonella — may serve as template."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "FAD is a standard organic cofactor. Default values should be correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── REDOX-SENSING REGULATORS ──────────────────────────────────────────────
    "oxyR_H2O2": {
        "gene": "oxyR",
        "protein_name": "OxyR",
        "ligand_name": "H2O2",
        "pdb_id": "",
        "ligand_residue": "",
        "pdb_notes": (
            "H2O2 reacts covalently with OxyR Cys199 — no stable co-crystal in PDB. "
            "Search for Pa OxyR apo structure. E. coli OxyR (PDB 1I6A) used as "
            "reference. H2O2 must be manually placed near Cys199 in the active site "
            "before running simulations."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "H2O2 is a standard small molecule. Default values should be correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "soxR_Fe2S2": {
        "gene": "soxR",
        "protein_name": "SoxR",
        "ligand_name": "[2Fe-2S]",
        "pdb_id": "",
        "ligand_residue": "SF4",
        "pdb_notes": (
            "No Pa-specific SoxR PDB — only E. coli SoxR structures exist (e.g. 2ZHH). "
            "These are NOT valid for Pa. Confirm with Shiva Sir whether to use E. coli "
            "structure as template or wait for Pa-specific data. "
            "[2Fe-2S] cluster in PDB is typically residue SF4 or FES."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "[2Fe-2S] cluster spin coupling is complex. Oxidized [2Fe-2S]2+ gives "
            "S=0 (antiferromagnetically coupled Fe3+-Fe3+, mult=1). "
            "Reduced [2Fe-2S]1+ gives S=1/2 (mult=2). "
            "Must confirm broken-symmetry DFT approach with Shiva Sir."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "mexR_H2O2": {
        "gene": "mexR",
        "protein_name": "MexR",
        "ligand_name": "H2O2",
        "pdb_id": "1LNW",
        "ligand_residue": "",
        "pdb_notes": (
            "PDB 1LNW is Pa MexR oxidised form (Cys30-Cys62 disulphide). "
            "No H2O2 co-crystal — place H2O2 near Cys30 manually. "
            "Also see PDB 4LFH for comparison."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "H2O2 effector — standard organic. Default values should be correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── SUPEROXIDE DISMUTASES (metal cofactor) ────────────────────────────────
    "sodA_Mn2": {
        "gene": "sodA",
        "protein_name": "Superoxide dismutase A (MnSOD)",
        "ligand_name": "Mn2+",
        "pdb_id": "",
        "ligand_residue": "MN",
        "pdb_notes": (
            "Search RCSB for Pa MnSOD (sodA). Try 'manganese superoxide dismutase "
            "Pseudomonas aeruginosa'. "
            "B. stearothermophilus MnSOD (1MSD) used as general reference."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "Mn(II) d5 high-spin: S=5/2, multiplicity=6. "
            "Mn(III) d4 high-spin: S=2, multiplicity=5. "
            "Set multiplicity based on the oxidation state you are modelling. "
            "Confirm with Shiva Sir."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 6},
            "radical": {"charge": 1,  "multiplicity": 5},
            "reduced": {"charge": -1, "multiplicity": 6},
        }
    },

    "sodB_Fe2": {
        "gene": "sodB",
        "protein_name": "Superoxide dismutase B (FeSOD)",
        "ligand_name": "Fe2+",
        "pdb_id": "",
        "ligand_residue": "FE2",
        "pdb_notes": (
            "Search RCSB for Pa FeSOD (sodB): try 'iron superoxide dismutase "
            "Pseudomonas aeruginosa'. PDB 1ISA is E. coli FeSOD — may serve as template."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "Fe(II) d6 high-spin: S=2, multiplicity=5. "
            "Fe(III) d5 high-spin: S=5/2, multiplicity=6. "
            "Confirm oxidation state and spin state with Shiva Sir."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 5},
            "radical": {"charge": 1,  "multiplicity": 6},
            "reduced": {"charge": -1, "multiplicity": 5},
        }
    },

    "sodM_Mn2": {
        "gene": "sodM",
        "protein_name": "Superoxide dismutase M (cambialistic)",
        "ligand_name": "Mn2+",
        "pdb_id": "",
        "ligand_residue": "MN",
        "pdb_notes": (
            "Cambialistic SOD — can use Mn2+ or Fe2+ depending on conditions. "
            "Search for Pa SodM. Confirm which metal form to simulate with Shiva Sir."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "Same considerations as sodA_Mn2 for Mn(II)/Mn(III). "
            "If simulating Fe form, same as sodB_Fe2. "
            "Confirm with Shiva Sir."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 6},
            "radical": {"charge": 1,  "multiplicity": 5},
            "reduced": {"charge": -1, "multiplicity": 6},
        }
    },

    # ── THIOREDOXIN SYSTEM ────────────────────────────────────────────────────
    "trxA_NADPH": {
        "gene": "trxA",
        "protein_name": "Thioredoxin 1",
        "ligand_name": "NADPH",
        "pdb_id": "",
        "ligand_residue": "NDP",
        "pdb_notes": (
            "TrxA does not directly bind NADPH — it is reduced by TrxB/NADPH. "
            "Consider simulating TrxB+NADPH instead (trxB_FAD). "
            "If simulating TrxA, search for Pa thioredoxin structure. "
            "NADPH in PDB is typically residue NDP or NAP."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "NADPH is a standard organic molecule. Default values should be correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "trxB_FAD": {
        "gene": "trxB",
        "protein_name": "Thioredoxin reductase",
        "ligand_name": "FAD",
        "pdb_id": "",
        "ligand_residue": "FAD",
        "pdb_notes": (
            "Search RCSB for Pa thioredoxin reductase. "
            "E. coli TrxB (PDB 1TDF) has FAD confirmed. "
            "Pa-specific structure: search 'thioredoxin reductase Pseudomonas aeruginosa'."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "FAD is a standard organic cofactor. Default values should be correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "trxC_NADPH": {
        "gene": "trxC",
        "protein_name": "Thioredoxin 2 (periplasmic)",
        "ligand_name": "NADPH",
        "pdb_id": "",
        "ligand_residue": "NDP",
        "pdb_notes": (
            "Periplasmic thioredoxin. Search for Pa TrxC structure. "
            "NADPH in PDB is typically residue NDP or NAP."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Standard organic molecule. Default values should be correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── GLUTAREDOXIN / GLUTATHIONE SYSTEM ────────────────────────────────────
    "grxA_GSH": {
        "gene": "grxA",
        "protein_name": "Glutaredoxin 1",
        "ligand_name": "glutathione",
        "pdb_id": "",
        "ligand_residue": "GSH",
        "pdb_notes": (
            "Search RCSB for Pa glutaredoxin 1. "
            "Glutathione in PDB is residue GSH (reduced) or GSS (oxidised)."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Glutathione is a standard organic tripeptide. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "grxB_GSH": {
        "gene": "grxB",
        "protein_name": "Glutaredoxin 2",
        "ligand_name": "glutathione",
        "pdb_id": "",
        "ligand_residue": "GSH",
        "pdb_notes": "Search RCSB for Pa glutaredoxin 2. Glutathione residue: GSH.",
        "review_multiplicity": False,
        "multiplicity_notes": "Standard organic molecule. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "gshA_GSH": {
        "gene": "gshA",
        "protein_name": "Glutamate-cysteine ligase",
        "ligand_name": "glutathione",
        "pdb_id": "",
        "ligand_residue": "GSH",
        "pdb_notes": (
            "GshA produces GSH but the product is the simulation target. "
            "Search for Pa GshA (glutamate-cysteine ligase). "
            "May need to manually place GSH in the active site."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Standard organic molecule. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "gshB_GSH": {
        "gene": "gshB",
        "protein_name": "Glutathione synthetase",
        "ligand_name": "glutathione",
        "pdb_id": "",
        "ligand_residue": "GSH",
        "pdb_notes": "Search for Pa GshB (glutathione synthetase). Place GSH manually.",
        "review_multiplicity": False,
        "multiplicity_notes": "Standard organic molecule. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "gor_FAD": {
        "gene": "gor",
        "protein_name": "Glutathione reductase",
        "ligand_name": "FAD",
        "pdb_id": "",
        "ligand_residue": "FAD",
        "pdb_notes": (
            "Search RCSB for Pa glutathione reductase (gor). "
            "Human GR (PDB 1GRA) has FAD+GSH confirmed — may serve as template. "
            "Also check for Pa-specific structures."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "FAD is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── NADH/NADPH METABOLISM ─────────────────────────────────────────────────
    "ndh_NADH": {
        "gene": "ndh",
        "protein_name": "NADH dehydrogenase II",
        "ligand_name": "NADH",
        "pdb_id": "",
        "ligand_residue": "NAI",
        "pdb_notes": (
            "Search RCSB for Pa NDH-II (type II NADH dehydrogenase). "
            "NADH in PDB is residue NAI (reduced) or NAD (oxidised). "
            "S. aureus NDH-II (PDB 4NWZ) has NADH — may serve as template."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "NADH is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "nuoA_NADH": {
        "gene": "nuoA",
        "protein_name": "NADH:ubiquinone oxidoreductase subunit A",
        "ligand_name": "NADH",
        "pdb_id": "",
        "ligand_residue": "NAI",
        "pdb_notes": (
            "Complex I subunit. Search RCSB for Pa Complex I or NuoA. "
            "T. thermophilus Complex I (PDB 4HEA) may serve as template."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "NADH is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "nuoB_NADH": {
        "gene": "nuoB",
        "protein_name": "NADH:ubiquinone oxidoreductase subunit B",
        "ligand_name": "NADH",
        "pdb_id": "",
        "ligand_residue": "NAI",
        "pdb_notes": "Complex I subunit NuoB. Same template notes as nuoA_NADH.",
        "review_multiplicity": False,
        "multiplicity_notes": "NADH is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "nuoN_UBQ": {
        "gene": "nuoN",
        "protein_name": "NADH:ubiquinone oxidoreductase subunit N",
        "ligand_name": "ubiquinone",
        "pdb_id": "",
        "ligand_residue": "UQ",
        "pdb_notes": (
            "Complex I proton-pumping module. Ubiquinone in PDB: try residue UQ, "
            "UBQ, or Q8 depending on chain length. Same template notes as nuoA_NADH."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Ubiquinone is a standard organic molecule. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "sdhA_FAD": {
        "gene": "sdhA",
        "protein_name": "Succinate dehydrogenase subunit A",
        "ligand_name": "FAD",
        "pdb_id": "",
        "ligand_residue": "FAD",
        "pdb_notes": (
            "Complex II FAD-binding subunit. Search for Pa SdhA. "
            "E. coli Complex II (PDB 1NEK) has FAD confirmed — may serve as template."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "FAD is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "sdhB_UBQ": {
        "gene": "sdhB",
        "protein_name": "Succinate dehydrogenase subunit B",
        "ligand_name": "ubiquinone",
        "pdb_id": "",
        "ligand_residue": "UQ",
        "pdb_notes": (
            "Complex II Fe-S subunit. Also contains [2Fe-2S], [3Fe-4S], [4Fe-4S] "
            "clusters. The simulation target is ubiquinone. "
            "Search for Pa SdhB. Same template as sdhA_FAD."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "SdhB contains multiple Fe-S clusters whose spin states affect the QM region. "
            "If Fe-S clusters are included in the QM region, multiplicity must be "
            "confirmed with Shiva Sir. If only ubiquinone is the focus, default values "
            "may be acceptable."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── CYTOCHROMES AND TERMINAL OXIDASES (heme cofactor) ────────────────────
    "cioA_heme": {
        "gene": "cioA",
        "protein_name": "Cyanide-insensitive oxidase subunit A",
        "ligand_name": "heme",
        "pdb_id": "",
        "ligand_residue": "HEM",
        "pdb_notes": (
            "Search for Pa CioA (cyanide-insensitive oxidase). "
            "Limited Pa-specific structures — check RCSB for homologs."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": "Heme-containing cytochrome — same Fe spin state considerations as katA_heme.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "cioB_heme": {
        "gene": "cioB",
        "protein_name": "Cyanide-insensitive oxidase subunit B",
        "ligand_name": "heme",
        "pdb_id": "",
        "ligand_residue": "HEM",
        "pdb_notes": "Same notes as cioA_heme.",
        "review_multiplicity": True,
        "multiplicity_notes": "Same Fe spin state considerations as katA_heme.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "cyoA_heme": {
        "gene": "cyoA",
        "protein_name": "Cytochrome bo3 oxidase subunit I",
        "ligand_name": "heme",
        "pdb_id": "",
        "ligand_residue": "HEM",
        "pdb_notes": (
            "Cytochrome bo3 contains heme b and heme o. "
            "Search for Pa CyoA or cytochrome bo3. "
            "E. coli CyoA (PDB 1FFT) may serve as template."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": "Heme-containing cytochrome — confirm Fe spin state with Shiva Sir.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "coxA_heme": {
        "gene": "coxA",
        "protein_name": "Cytochrome c oxidase subunit I",
        "ligand_name": "heme",
        "pdb_id": "",
        "ligand_residue": "HEM",
        "pdb_notes": (
            "Cytochrome c oxidase contains heme a and Cu centres. "
            "Search for Pa CoxA. Bovine CcO (PDB 1OCC) as general reference."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "Contains heme a/a3 AND Cu centers (CuA, CuB). "
            "The QM region spin state is highly complex. "
            "Confirm approach with Shiva Sir before running."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── FERREDOXINS ───────────────────────────────────────────────────────────
    "fdxA_Fdx": {
        "gene": "fdxA",
        "protein_name": "Ferredoxin I",
        "ligand_name": "ferredoxin ([2Fe-2S])",
        "pdb_id": "",
        "ligand_residue": "SF4",
        "pdb_notes": (
            "Ferredoxin contains [2Fe-2S] cluster. In PDB, [2Fe-2S] is residue SF4 or FES. "
            "Search for Pa FdxA (ferredoxin I). PDB 1FR5 (Azotobacter) as template reference."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "Same [2Fe-2S] spin coupling issues as soxR_Fe2S2. "
            "Oxidised [2Fe-2S]2+: S=0, mult=1. "
            "Reduced [2Fe-2S]1+: S=1/2, mult=2. "
            "Confirm broken-symmetry approach with Shiva Sir."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "fdxB_Fdx": {
        "gene": "fdxB",
        "protein_name": "Ferredoxin II",
        "ligand_name": "ferredoxin ([2Fe-2S])",
        "pdb_id": "",
        "ligand_residue": "SF4",
        "pdb_notes": "Same notes as fdxA_Fdx. Search for Pa FdxB (ferredoxin II).",
        "review_multiplicity": True,
        "multiplicity_notes": "Same [2Fe-2S] considerations as fdxA_Fdx.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── IRON METABOLISM ───────────────────────────────────────────────────────
    "fur_Fe2": {
        "gene": "fur",
        "protein_name": "Ferric uptake regulator Fur",
        "ligand_name": "Fe2+",
        "pdb_id": "",
        "ligand_residue": "FE2",
        "pdb_notes": (
            "Search for Pa Fur structure: try PDB 2IDA or 'Fur Pseudomonas aeruginosa'. "
            "Fe2+ is often removed during crystallisation — may need to add manually. "
            "Zn2+ structural site is also present but is NOT the simulation target. "
            "Fe2+ in PDB: residue FE2 or FE."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": (
            "Fe(II) in Fur is a regulatory corepressor. Fe(II) d6 high-spin: S=2, mult=5. "
            "Fe(II) low-spin: S=0, mult=1. Confirm spin state with Shiva Sir."
        ),
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 5},
            "radical": {"charge": 1,  "multiplicity": 6},
            "reduced": {"charge": -1, "multiplicity": 5},
        }
    },

    "pvdA_FAD": {
        "gene": "pvdA",
        "protein_name": "Pyoverdine biosynthesis protein A",
        "ligand_name": "FAD",
        "pdb_id": "",
        "ligand_residue": "FAD",
        "pdb_notes": (
            "PvdA is L-ornithine hydroxylase (FAD+NADPH+O2 dependent). "
            "Search for Pa PvdA: try PDB 3S0P or 'ornithine hydroxylase Pseudomonas'."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "FAD is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "bfrA_heme": {
        "gene": "bfrA",
        "protein_name": "Bacterioferritin A",
        "ligand_name": "heme",
        "pdb_id": "",
        "ligand_residue": "HEM",
        "pdb_notes": (
            "Search for Pa BfrA (bacterioferritin A). "
            "Note: bfrA in Pa encodes FtnA (bacterial ferritin), NOT a true bacterioferritin — "
            "it may not contain heme. Confirm the correct structure and heme status "
            "with Shiva Sir before running."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": "Heme-containing — confirm Fe spin state if heme is confirmed.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "bfrB_heme": {
        "gene": "bfrB",
        "protein_name": "Bacterioferritin B",
        "ligand_name": "heme",
        "pdb_id": "",
        "ligand_residue": "HEM",
        "pdb_notes": (
            "Search for Pa BfrB (bacterioferritin B). BfrB contains heme b confirmed. "
            "Try 'bacterioferritin Pseudomonas aeruginosa' in RCSB. PDB 3IS7 may match."
        ),
        "review_multiplicity": True,
        "multiplicity_notes": "Heme-containing — confirm Fe spin state with Shiva Sir.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── PHENAZINE BIOSYNTHESIS ─────────────────────────────────────────────────
    "phzM_SAM": {
        "gene": "phzM",
        "protein_name": "Phenazine methyltransferase",
        "ligand_name": "SAM",
        "pdb_id": "2RGF",
        "ligand_residue": "SAM",
        "pdb_notes": (
            "PDB 2RGF — Pa PhzM with SAM (S-adenosylmethionine) confirmed. "
            "SAM in PDB is residue SAM or SAH (hydrolysed product)."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "SAM is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "phzS_FAD": {
        "gene": "phzS",
        "protein_name": "Phenazine hydroxylase",
        "ligand_name": "FAD",
        "pdb_id": "2RGE",
        "ligand_residue": "FAD",
        "pdb_notes": (
            "PDB 2RGE — Pa PhzS crystal structure. FAD cofactor confirmed. "
            "Verify FAD residue name (FAD) in the PDB HETATM records."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "FAD is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    # ── QUORUM SENSING ────────────────────────────────────────────────────────
    "lasR_3OC12HSL": {
        "gene": "lasR",
        "protein_name": "LasR",
        "ligand_name": "3OC12-HSL",
        "pdb_id": "3IX3",
        "ligand_residue": "OHL",
        "pdb_notes": (
            "PDB 3IX3 — Pa LasR LBD with 3-oxo-C12-HSL co-crystal confirmed. "
            "Ligand residue is 'OHL' in the PDB HETATM records. "
            "Also see PDB 2UV0 for an alternative structure."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Organic autoinducer — standard molecule. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "lasI_3OC12HSL": {
        "gene": "lasI",
        "protein_name": "LasI",
        "ligand_name": "3OC12-HSL",
        "pdb_id": "",
        "ligand_residue": "OHL",
        "pdb_notes": (
            "LasI is the 3OC12-HSL synthase — 3OC12-HSL is its product. "
            "Search for Pa LasI structure: try 'LasI Pseudomonas aeruginosa'. "
            "If no Pa structure exists, may need to use homolog as template."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Organic autoinducer. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "rhlR_C4HSL": {
        "gene": "rhlR",
        "protein_name": "RhlR",
        "ligand_name": "C4-HSL",
        "pdb_id": "8B4A",
        "ligand_residue": "C4H",
        "pdb_notes": (
            "PDB 8B4A — Pa RhlR with C4-HSL co-crystal confirmed. "
            "Also see PDB 3T5K and 7R3H. "
            "Verify the exact ligand residue name (likely C4H or BHL) in the HETATM records."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Organic autoinducer. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "rhlI_C4HSL": {
        "gene": "rhlI",
        "protein_name": "RhlI",
        "ligand_name": "C4-HSL",
        "pdb_id": "",
        "ligand_residue": "C4H",
        "pdb_notes": (
            "RhlI is the C4-HSL synthase. Search for Pa RhlI structure. "
            "C4-HSL residue name is likely C4H — verify in PDB viewer."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "Organic autoinducer. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "mvfR_HHQ": {
        "gene": "mvfR",
        "protein_name": "MvfR (PqsR)",
        "ligand_name": "HHQ",
        "pdb_id": "",
        "ligand_residue": "HHQ",
        "pdb_notes": (
            "HHQ = 2-heptyl-4-hydroxyquinoline. PqsR/MvfR binds HHQ and PQS. "
            "Search RCSB for PqsR: try 'PqsR Pseudomonas aeruginosa HHQ'. "
            "HHQ residue name in PDB may vary — check HETATM records."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "HHQ is a standard organic molecule. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },

    "pqsA_CoA": {
        "gene": "pqsA",
        "protein_name": "PqsA",
        "ligand_name": "CoA",
        "pdb_id": "",
        "ligand_residue": "COA",
        "pdb_notes": (
            "PqsA is anthraniloyl-CoA ligase. CoA in PDB is residue COA. "
            "Search for Pa PqsA structure."
        ),
        "review_multiplicity": False,
        "multiplicity_notes": "CoA is a standard organic cofactor. Default values correct.",
        "states": {
            "neutral": {"charge": 0,  "multiplicity": 1},
            "radical": {"charge": 1,  "multiplicity": 2},
            "reduced": {"charge": -1, "multiplicity": 2},
        }
    },
}


def create_config(protein_id: str, entry: dict, force: bool, dry_run: bool) -> str:
    """Create config.json for a single protein. Returns status string."""
    folder = PROTEINS_DIR / protein_id
    config_path = folder / "config.json"

    if config_path.exists() and not force:
        return "exists"

    if dry_run:
        return "would_create"

    folder.mkdir(parents=True, exist_ok=True)

    config = {
        "protein_id":          protein_id,
        "gene":                entry["gene"],
        "protein_name":        entry["protein_name"],
        "ligand_name":         entry["ligand_name"],
        "pdb_id":              entry["pdb_id"],
        "ligand_residue":      entry["ligand_residue"],
        "pdb_notes":           entry["pdb_notes"],
        "review_multiplicity": entry["review_multiplicity"],
        "multiplicity_notes":  entry["multiplicity_notes"],
        "states":              entry["states"],
    }

    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    return "created"


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create protein folders and config.json files for all simulatable "
            "P. aeruginosa redox/AMR proteins."
        )
    )
    parser.add_argument(
        "--only", nargs="+",
        help="Create only specific protein IDs (e.g. --only katA_heme lasR_3OC12HSL)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Preview what would be created without writing any files"
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Overwrite existing config.json files"
    )
    args = parser.parse_args()

    target_ids = args.only if args.only else list(PROTEIN_REGISTRY.keys())

    # Validate --only args
    unknown = [pid for pid in target_ids if pid not in PROTEIN_REGISTRY]
    if unknown:
        print(f"ERROR: Unknown protein IDs: {unknown}")
        print(f"Valid IDs: {list(PROTEIN_REGISTRY.keys())}")
        sys.exit(1)

    if args.dry_run:
        print("\n[DRY RUN — no files will be written]\n")

    PROTEINS_DIR.mkdir(parents=True, exist_ok=True)

    created, existed, would_create = [], [], []
    needs_pdb    = []
    needs_review = []

    print(f"\nSetting up {len(target_ids)} protein(s)...\n")

    for protein_id in target_ids:
        entry  = PROTEIN_REGISTRY[protein_id]
        status = create_config(protein_id, entry, args.force, args.dry_run)

        if status == "created":
            created.append(protein_id)
            print(f"  [CREATED]  {protein_id:<22}  pdb={entry['pdb_id'] or '(not set)'}")
        elif status == "exists":
            existed.append(protein_id)
            print(f"  [EXISTS]   {protein_id:<22}  (use --force to overwrite)")
        elif status == "would_create":
            would_create.append(protein_id)
            print(f"  [WOULD CREATE]  {protein_id:<22}  pdb={entry['pdb_id'] or '(not set)'}")

        if not entry["pdb_id"]:
            needs_pdb.append(protein_id)
        if entry["review_multiplicity"]:
            needs_review.append(protein_id)

    print(f"\n{'='*60}")
    if args.dry_run:
        print(f"Dry run summary: {len(would_create)} would be created, {len(existed)} already exist.")
    else:
        print(f"Summary: {len(created)} created, {len(existed)} already existed.")

    if needs_pdb:
        print(f"\n{'─'*60}")
        print(f"ACTION REQUIRED — {len(needs_pdb)} proteins need a PDB ID set manually:")
        for pid in needs_pdb:
            print(f"  {pid}")
            print(f"    → {PROTEIN_REGISTRY[pid]['pdb_notes'][:80]}...")
        print(f"\n  Edit proteins/{{protein_id}}/config.json and fill in 'pdb_id'")
        print(f"  and 'ligand_residue'. Then run fetch_structures.py.")

    if needs_review:
        print(f"\n{'─'*60}")
        print(f"REVIEW REQUIRED — {len(needs_review)} proteins have metal centers:")
        print(f"  Charge/multiplicity defaults may be WRONG for these proteins.")
        print(f"  Confirm ALL values with Shiva Sir before running simulations:")
        for pid in needs_review:
            print(f"  {pid}")
        print(f"\n  See 'multiplicity_notes' in each config.json for details.")

    print(f"\n{'='*60}")
    print(f"Next step: run fetch_structures.py to download PDB files.")
    print(f"  For proteins with empty pdb_id: fill in config.json first.\n")


if __name__ == "__main__":
    main()
