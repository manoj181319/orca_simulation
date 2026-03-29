"""
score_targets.py  —  Phase 0, Step 3
Reads extracted targets and scores each gene/protein for simulation priority.

Scoring formula (weights defined as constants at top of file):
    score = (W_CITATION  * normalised_citation_frequency)
          + (W_LIGAND    * has_simulatable_ligand)
          + (W_PDB       * has_pdb_structure)
          + (W_CATEGORY  * category_weight)
          + (W_COMENTION * normalised_signal_ligand_diversity)

All terms are normalised to [0, 1]. Final score is always in [0, 1].

Ligand assignment:
    top_ligand comes directly from the 'top_ligand' field in extracted_targets.json,
    which is set by GENE_BIOCHEMISTRY in extract_targets.py (biochemistry-first).
    Co-mention statistics are used ONLY for the diversity bonus (W_COMENTION).
    Noise ligands (pyocyanin, iron, phenazine, etc.) are excluded from diversity.

Input:
    literature_mining/extracted/extracted_targets.json

Output:
    literature_mining/extracted/scored_targets.json

Usage:
    python score_targets.py
"""

import json
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_FILE = Path(__file__).parent / "extracted" / "extracted_targets.json"
OUTPUT_FILE = Path(__file__).parent / "extracted" / "scored_targets.json"

# ---------------------------------------------------------------------------
# Scoring weights — must sum to 1.0
# ---------------------------------------------------------------------------
W_CITATION   = 0.35   # How frequently the gene appears across papers
W_LIGAND     = 0.20   # Whether a known, simulatable ligand is co-mentioned
W_PDB        = 0.20   # Whether a crystal structure exists in PDB
W_CATEGORY   = 0.15   # Biological category relevance to redox AMR
W_COMENTION  = 0.10   # Diversity of ligands co-mentioned with this gene

assert abs(W_CITATION + W_LIGAND + W_PDB + W_CATEGORY + W_COMENTION - 1.0) < 1e-9, \
    "Weights must sum to 1.0"

# ---------------------------------------------------------------------------
# Category relevance weights [0.0 – 1.0]
# Based on how directly each category links redox biology to AMR
# ---------------------------------------------------------------------------
CATEGORY_WEIGHTS = {
    "oxidative_stress":   1.0,   # Direct ROS resistance = core AMR mechanism
    "electron_transport": 0.9,   # PMF drives efflux; respiratory adaptation = AMR
    "redox_enzyme":       0.85,  # NADPH/thiol balance — drug metabolism
    "efflux":             0.80,  # PMF-driven; most clinically relevant AMR class
    "iron_metabolism":    0.70,  # Fenton chemistry, siderophore competition
    "regulation":         0.65,  # Upstream controllers of redox/AMR genes
    "phenazine":          0.60,  # Electron shuttles — biofilm redox
    "quorum_sensing":     0.50,  # Indirect — coordinates redox virulence factors
}

# ---------------------------------------------------------------------------
# PDB availability — proteins with at least one crystal structure in PDB.
# NOTE: This only means the protein structure exists. It does NOT mean a
# docked complex with the simulation ligand is available.
# needs_docking is determined by HAS_COCRYSTAL_WITH_LIGAND below.
# ---------------------------------------------------------------------------
HAS_PDB_STRUCTURE = {
    # Confirmed P. aeruginosa PAO1 crystal structures in RCSB PDB (March 2026).
    # soxR excluded: only E. coli SoxR structures exist; Pa-specific not confirmed.
    "katA", "katB", "ahpC", "ahpF", "oxyR",
    "sodA", "sodB", "sodM",
    "trxA", "trxB", "grxA", "gor",
    "nuoA", "nuoB", "sdhA", "sdhB", "cioA", "fdxA",
    "fur", "pvdA", "pvdD",
    "mexA", "mexB", "oprM", "mexC", "mexD",
    "mexR",   # PDB 1LNW (oxidised), 4LFH — Pa MexR confirmed
    "lasR", "lasI", "rhlR", "rhlI", "mvfR",
    "phzM", "phzS",
}

# ---------------------------------------------------------------------------
# Co-crystal structures: (gene, ligand) pairs where a PDB entry contains the
# specific simulation ligand bound in the active site of the Pa protein.
# ONLY these pairs are marked "No docking needed" in the report.
# Conservative by design — when in doubt, always require docking.
# ---------------------------------------------------------------------------
HAS_COCRYSTAL_WITH_LIGAND = {
    ("katA",  "heme"),       # PDB 4E37 — Pa PAO1 KatA + HEM confirmed
    ("katA",  "NADPH"),      # PDB 4E37 — NADPH also present at peripheral site
    ("trxB",  "FAD"),        # TrxB flavoprotein structures with FAD confirmed
    ("gor",   "FAD"),        # Glutathione reductase flavoprotein with FAD
    ("gor",   "glutathione"),# Glutathione reductase with GSH substrate
    ("sdhA",  "FAD"),        # Complex II SdhA with FAD confirmed
    ("sdhA",  "ubiquinone"), # Complex II with ubiquinone
    ("sdhB",  "ubiquinone"), # Complex II Fe-S to ubiquinone
    ("lasR",  "3OC12-HSL"),  # PDB 2UV0, 3IX3 — LasR LBD with autoinducer
    ("rhlR",  "C4-HSL"),     # PDB 3T5K, 7R3H, 8B4A — RhlR with C4-HSL
}

# ---------------------------------------------------------------------------
# Simulatable ligands — small molecules suitable for ORCA DFT/TDDFT simulation.
# Used only to determine whether a gene gets the W_LIGAND score bonus.
# top_ligand itself always comes from GENE_BIOCHEMISTRY via the JSON.
# ---------------------------------------------------------------------------
SIMULATABLE_LIGANDS = {
    "NADH", "NADPH", "FAD", "FMN", "heme", "haem", "ubiquinone", "menaquinone",
    "glutathione", "ferredoxin", "CoA", "coenzyme A", "flavin", "quinone", "SAM",
    "3OC12-HSL", "3-oxo-C12-HSL", "C4-HSL", "N-butanoyl-HSL",
    "H2O2", "hydrogen peroxide", "organic hydroperoxide",
    "Fe2+", "Mn2+", "[2Fe-2S]", "HHQ", "PQS",
    # Note: pyocyanin, iron, phenazine intentionally excluded (LIGAND_NOISE)
}


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------
def normalise(values: list) -> list:
    if not values:
        return values
    lo, hi = min(values), max(values)
    if hi == lo:
        return [0.0] * len(values)
    return [(v - lo) / (hi - lo) for v in values]


def ligand_diversity_score(ligands_co_mentioned: dict) -> float:
    """Count distinct SIGNAL (non-noise) simulatable ligands co-mentioned (cap 5).
    Noise ligands are already separated in the JSON by extract_targets.py."""
    count = sum(1 for lig in ligands_co_mentioned if lig in SIMULATABLE_LIGANDS)
    return min(count, 5) / 5.0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if not INPUT_FILE.exists():
        print(f"ERROR: {INPUT_FILE} not found. Run extract_targets.py first.")
        return

    with open(INPUT_FILE, encoding="utf-8") as f:
        data = json.load(f)

    genes = data.get("genes", {})

    if not genes:
        print("No genes found in extracted_targets.json. Exiting.")
        return

    print(f"Scoring {len(genes)} genes...\n")

    # --- Raw values for normalisation ---
    citation_counts    = [d["paper_count"] for d in genes.values()]
    # Use signal co-mentions only (noise ligands already separated in JSON)
    ligand_diversity_raw = [
        ligand_diversity_score(d.get("ligands_co_mentioned", {}))
        for d in genes.values()
    ]

    # Normalise citation counts
    norm_citations = normalise(citation_counts)

    scored = []
    for i, (gene, data) in enumerate(genes.items()):
        category = data.get("category", "")
        # top_ligand and simulatable come from GENE_BIOCHEMISTRY via the JSON.
        # Co-mention stats are used only for the diversity bonus (W_COMENTION).
        best_ligand  = data.get("top_ligand", "")
        is_simulatable = data.get("simulatable", False)
        # Signal co-mentions only (noise ligands already separated by extract_targets)
        signal_ligands = data.get("ligands_co_mentioned", {})

        citation_score  = norm_citations[i]
        ligand_score    = 1.0 if is_simulatable else 0.0
        pdb_score       = 1.0 if gene in HAS_PDB_STRUCTURE else 0.0
        category_score  = CATEGORY_WEIGHTS.get(category, 0.3)
        comention_score = ligand_diversity_raw[i]

        final_score = (
            W_CITATION  * citation_score +
            W_LIGAND    * ligand_score +
            W_PDB       * pdb_score +
            W_CATEGORY  * category_score +
            W_COMENTION * comention_score
        )

        scored.append({
            "gene":         gene,
            "protein":      data.get("protein", ""),
            "function":     data.get("function", ""),
            "category":     category,
            "paper_count":  data["paper_count"],
            "has_pdb":      gene in HAS_PDB_STRUCTURE,
            "simulatable":  is_simulatable,
            "top_ligand":   best_ligand,
            "ligand_role":  data.get("ligand_role", ""),
            "ligand_basis": data.get("ligand_basis", ""),
            "ligand_note":  data.get("ligand_note", ""),
            "ligands_co_mentioned":       signal_ligands,
            "ligands_noise_co_mentioned": data.get("ligands_noise_co_mentioned", {}),
            "sample_papers": data.get("sample_papers", []),
            "score_breakdown": {
                "citation":  round(W_CITATION  * citation_score,  4),
                "ligand":    round(W_LIGAND    * ligand_score,     4),
                "pdb":       round(W_PDB       * pdb_score,        4),
                "category":  round(W_CATEGORY  * category_score,   4),
                "comention": round(W_COMENTION * comention_score,   4),
            },
            "final_score": round(final_score, 4),
        })

    # Sort by final score descending
    scored.sort(key=lambda x: -x["final_score"])
    for rank, item in enumerate(scored, 1):
        item["rank"] = rank

    # Recommended pairs: only genes marked simulatable=True in GENE_BIOCHEMISTRY.
    # needs_docking = True unless a confirmed co-crystal exists for (gene, ligand).
    recommended_pairs = [
        {
            "rank":             item["rank"],
            "gene":             item["gene"],
            "protein":          item["protein"],
            "suggested_ligand": item["top_ligand"],
            "ligand_role":      item["ligand_role"],
            "ligand_basis":     item["ligand_basis"],
            "has_pdb":          item["has_pdb"],
            "score":            item["final_score"],
            "needs_docking":    (item["gene"], item["top_ligand"]) not in HAS_COCRYSTAL_WITH_LIGAND,
        }
        for item in scored
        if item["simulatable"]
    ][:20]

    output = {
        "total_genes_scored": len(scored),
        "scoring_weights": {
            "citation":  W_CITATION,
            "ligand":    W_LIGAND,
            "pdb":       W_PDB,
            "category":  W_CATEGORY,
            "comention": W_COMENTION,
        },
        "ranked_genes": scored,
        "recommended_pairs": recommended_pairs,
    }

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f"Scoring complete. Top 15 candidates:\n")
    print(f"  {'Rank':<5} {'Gene':<8} {'Score':<7} {'PDB':<5} {'Sim':<5} {'Category':<20} {'Ligand':<22} {'Protein'}")
    print(f"  {'-'*100}")
    for item in scored[:15]:
        pdb = "Yes" if item["has_pdb"] else "No"
        sim = "Yes" if item["simulatable"] else "No"
        print(
            f"  {item['rank']:<5} {item['gene']:<8} {item['final_score']:<7.4f} "
            f"{pdb:<5} {sim:<5} {item['category']:<20} {item['top_ligand']:<22} {item['protein']}"
        )

    print(f"\nSaved to: {OUTPUT_FILE}")


if __name__ == '__main__':
    main()
