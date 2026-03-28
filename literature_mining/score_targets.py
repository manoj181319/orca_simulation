"""
score_targets.py  —  Phase 0, Step 3
Reads extracted targets and scores each gene/protein for simulation priority.

Scoring formula (weights defined as constants at top of file):
    score = (W_CITATION  * normalised_citation_frequency)
          + (W_LIGAND    * has_known_ligand)
          + (W_PDB       * has_pdb_structure)
          + (W_CATEGORY  * category_weight)
          + (W_COMENTION * normalised_ligand_comention_count)

All terms are normalised to [0, 1] before weighting so the final score
is always in [0, 1].

Category weights reflect biological relevance to redox-driven AMR:
    oxidative_stress > electron_transport > redox_enzyme > efflux >
    iron_metabolism > regulation > phenazine > quorum_sensing

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
# PDB availability — manually curated for known P. aeruginosa redox proteins
# This is a static list; update when new structures are deposited.
# Source: PDB search performed March 2026
# ---------------------------------------------------------------------------
HAS_PDB_STRUCTURE = {
    "katA", "katB", "ahpC", "ahpF", "oxyR", "soxR",
    "sodA", "sodB", "sodM",
    "trxA", "trxB",
    "grxA",
    "gor",
    "nuoA", "nuoB",
    "sdhA", "sdhB",
    "cioA",
    "fdxA",
    "fur",
    "pvdA", "pvdD",
    "mexA", "mexB", "oprM",
    "mexC", "mexD",
    "lasR", "lasI",
    "rhlR", "rhlI",
    "mvfR",
    "phzM", "phzS",
}

# ---------------------------------------------------------------------------
# Ligands that are simulatable with ORCA (confirmed by Shiva Sir or standard)
# A gene co-mentioning one of these gets the ligand bonus
# ---------------------------------------------------------------------------
SIMULATABLE_LIGANDS = {
    "NADH", "NADPH", "FAD", "FMN", "heme", "haem", "ubiquinone",
    "pyocyanin", "phenazine", "glutathione", "ferredoxin",
    "coenzyme A", "CoA", "flavin", "quinone",
}


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------
def normalise(values: list) -> list:
    """Min-max normalise a list of floats to [0, 1]."""
    if not values:
        return values
    min_v = min(values)
    max_v = max(values)
    if max_v == min_v:
        return [0.0] * len(values)
    return [(v - min_v) / (max_v - min_v) for v in values]


def has_simulatable_ligand(ligands_co_mentioned: dict) -> bool:
    return any(lig in SIMULATABLE_LIGANDS for lig in ligands_co_mentioned)


def top_ligand(ligands_co_mentioned: dict) -> str:
    if not ligands_co_mentioned:
        return ""
    return max(ligands_co_mentioned, key=ligands_co_mentioned.get)


def ligand_diversity_score(ligands_co_mentioned: dict) -> float:
    """Number of distinct simulatable ligands co-mentioned (capped at 5)."""
    count = sum(
        1 for lig in ligands_co_mentioned if lig in SIMULATABLE_LIGANDS
    )
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

    # --- Raw values ---
    citation_counts = [d["paper_count"] for d in genes.values()]
    ligand_diversity_raw = [
        ligand_diversity_score(d.get("ligands_co_mentioned", {}))
        for d in genes.values()
    ]

    # Normalise citation counts
    norm_citations = normalise(citation_counts)

    scored = []
    for i, (gene, data) in enumerate(genes.items()):
        category = data.get("category", "")
        ligands_co = data.get("ligands_co_mentioned", {})

        citation_score   = norm_citations[i]
        ligand_score     = 1.0 if has_simulatable_ligand(ligands_co) else 0.0
        pdb_score        = 1.0 if gene in HAS_PDB_STRUCTURE else 0.0
        category_score   = CATEGORY_WEIGHTS.get(category, 0.3)
        comention_score  = ligand_diversity_raw[i]

        final_score = (
            W_CITATION  * citation_score +
            W_LIGAND    * ligand_score +
            W_PDB       * pdb_score +
            W_CATEGORY  * category_score +
            W_COMENTION * comention_score
        )

        scored.append({
            "gene": gene,
            "protein": data.get("protein", ""),
            "function": data.get("function", ""),
            "category": category,
            "paper_count": data["paper_count"],
            "has_pdb": gene in HAS_PDB_STRUCTURE,
            "has_simulatable_ligand": has_simulatable_ligand(ligands_co),
            "top_ligand": top_ligand(ligands_co),
            "ligands_co_mentioned": ligands_co,
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

    # Recommended protein-ligand pairs:
    # Top candidates that have both a known PDB structure and a simulatable ligand
    recommended_pairs = [
        {
            "rank": item["rank"],
            "gene": item["gene"],
            "protein": item["protein"],
            "suggested_ligand": item["top_ligand"],
            "has_pdb": item["has_pdb"],
            "score": item["final_score"],
            "needs_docking": not item["has_pdb"],
        }
        for item in scored
        if item["has_simulatable_ligand"]
    ][:20]   # Top 20 recommended pairs for the report

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
    print(f"  {'Rank':<5} {'Gene':<8} {'Score':<7} {'PDB':<5} {'Ligand':<6} {'Category':<20} {'Protein'}")
    print(f"  {'-'*90}")
    for item in scored[:15]:
        pdb = "Yes" if item["has_pdb"] else "No"
        lig = "Yes" if item["has_simulatable_ligand"] else "No"
        print(
            f"  {item['rank']:<5} {item['gene']:<8} {item['final_score']:<7.4f} "
            f"{pdb:<5} {lig:<6} {item['category']:<20} {item['protein']}"
        )

    print(f"\nSaved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
