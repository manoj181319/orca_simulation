"""
extract_targets.py  —  Phase 0, Step 2
Reads raw JSON files from pubmed_fetch.py and extracts gene names, protein names,
and ligand mentions relevant to P. aeruginosa redox biology and AMR.

Strategy:
- Uses a curated known-gene list as the primary backbone (high precision)
- Supplements with regex patterns for novel/unlisted protein names (recall)
- Counts paper frequency per target and stores sample paper titles

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
INPUT_DIR = Path(__file__).parent / "raw_papers"
OUTPUT_DIR = Path(__file__).parent / "extracted"

PUBMED_RAW = INPUT_DIR / "pubmed_raw.json"
EPMC_RAW = INPUT_DIR / "epmc_raw.json"
OUTPUT_FILE = OUTPUT_DIR / "extracted_targets.json"

# ---------------------------------------------------------------------------
# Curated known-gene list for P. aeruginosa redox biology and AMR
#
# Format: { "gene_name": { "protein_name": str, "function": str, "category": str } }
#
# Categories:
#   oxidative_stress  — direct ROS scavenging / sensing
#   electron_transport — respiratory chain components
#   redox_enzyme      — NADH/NADPH-linked or thiol-redox enzymes
#   efflux            — PMF-driven drug efflux
#   iron_metabolism   — iron/metal redox homeostasis
#   regulation        — redox-sensing transcription factors
#   phenazine         — pyocyanin and phenazine biosynthesis/redox
#   quorum_sensing    — QS linked to redox state
# ---------------------------------------------------------------------------
KNOWN_GENES = {
    # Catalases and peroxidases
    "katA":  {"protein": "Catalase KatA",                   "function": "H2O2 decomposition (primary catalase)",          "category": "oxidative_stress"},
    "katB":  {"protein": "Catalase KatB",                   "function": "H2O2 decomposition (inducible catalase)",        "category": "oxidative_stress"},
    "katE":  {"protein": "Catalase KatE",                   "function": "H2O2 decomposition (stationary phase)",          "category": "oxidative_stress"},
    "ahpC":  {"protein": "Alkyl hydroperoxide reductase C",  "function": "Organic peroxide reduction",                    "category": "oxidative_stress"},
    "ahpF":  {"protein": "Alkyl hydroperoxide reductase F",  "function": "AhpC reductase subunit",                        "category": "oxidative_stress"},
    "oxyR":  {"protein": "OxyR",                             "function": "H2O2-sensing transcription factor",             "category": "regulation"},
    "soxR":  {"protein": "SoxR",                             "function": "Superoxide-sensing transcription factor",       "category": "regulation"},
    "soxS":  {"protein": "SoxS",                             "function": "SoxR-regulated AMR gene activator",             "category": "regulation"},

    # Superoxide dismutases
    "sodA":  {"protein": "Superoxide dismutase A (MnSOD)",   "function": "O2•- dismutation",                              "category": "oxidative_stress"},
    "sodB":  {"protein": "Superoxide dismutase B (FeSOD)",   "function": "O2•- dismutation",                              "category": "oxidative_stress"},
    "sodM":  {"protein": "Superoxide dismutase M (Mn/FeSOD)","function": "O2•- dismutation (manganese/iron)",             "category": "oxidative_stress"},

    # Thioredoxin system
    "trxA":  {"protein": "Thioredoxin 1",                    "function": "Thiol-disulfide oxidoreductase",                "category": "redox_enzyme"},
    "trxB":  {"protein": "Thioredoxin reductase",            "function": "NADPH-dependent Trx regeneration",              "category": "redox_enzyme"},
    "trxC":  {"protein": "Thioredoxin 2",                    "function": "Thiol-disulfide oxidoreductase (periplasmic)",  "category": "redox_enzyme"},

    # Glutaredoxin / glutathione system
    "grxA":  {"protein": "Glutaredoxin 1",                   "function": "Glutathione-dependent oxidoreductase",          "category": "redox_enzyme"},
    "grxB":  {"protein": "Glutaredoxin 2",                   "function": "Glutathione-dependent oxidoreductase",          "category": "redox_enzyme"},
    "gshA":  {"protein": "Glutamate-cysteine ligase",         "function": "Glutathione biosynthesis step 1",              "category": "redox_enzyme"},
    "gshB":  {"protein": "Glutathione synthetase",            "function": "Glutathione biosynthesis step 2",              "category": "redox_enzyme"},
    "gor":   {"protein": "Glutathione reductase",             "function": "NADPH-dependent GSH regeneration",             "category": "redox_enzyme"},

    # NADH/NADPH metabolism
    "ndh":   {"protein": "NADH dehydrogenase II",             "function": "Non-proton-pumping NADH oxidase",              "category": "electron_transport"},
    "nuoA":  {"protein": "NADH:ubiquinone oxidoreductase (complex I) subunit A", "function": "Proton-pumping NADH oxidase", "category": "electron_transport"},
    "nuoB":  {"protein": "NADH:ubiquinone oxidoreductase subunit B", "function": "Complex I — electron input",            "category": "electron_transport"},
    "nuoN":  {"protein": "NADH:ubiquinone oxidoreductase subunit N", "function": "Complex I — proton translocation",      "category": "electron_transport"},
    "sdhA":  {"protein": "Succinate dehydrogenase subunit A", "function": "Complex II — succinate oxidation",             "category": "electron_transport"},
    "sdhB":  {"protein": "Succinate dehydrogenase subunit B", "function": "Complex II — Fe-S cluster",                   "category": "electron_transport"},

    # Cytochromes and terminal oxidases
    "cioA":  {"protein": "Cyanide-insensitive oxidase subunit A", "function": "Terminal oxidase (antibiotic tolerance)",  "category": "electron_transport"},
    "cioB":  {"protein": "Cyanide-insensitive oxidase subunit B", "function": "Terminal oxidase",                         "category": "electron_transport"},
    "cyoA":  {"protein": "Cytochrome bo3 oxidase subunit I",  "function": "Terminal oxidase (aerobic)",                  "category": "electron_transport"},
    "coxA":  {"protein": "Cytochrome c oxidase subunit I",    "function": "Terminal oxidase (microaerobic)",              "category": "electron_transport"},

    # Ferredoxins
    "fdxA":  {"protein": "Ferredoxin I",                      "function": "Electron carrier",                             "category": "redox_enzyme"},
    "fdxB":  {"protein": "Ferredoxin II",                     "function": "Electron carrier",                             "category": "redox_enzyme"},

    # Iron metabolism
    "fur":   {"protein": "Ferric uptake regulator Fur",       "function": "Iron-sensing transcription factor",            "category": "iron_metabolism"},
    "pvdA":  {"protein": "Pyoverdine biosynthesis protein A", "function": "Siderophore biosynthesis",                     "category": "iron_metabolism"},
    "pvdD":  {"protein": "Pyoverdine synthetase D",           "function": "Siderophore biosynthesis",                     "category": "iron_metabolism"},
    "pchA":  {"protein": "Salicylate biosynthesis protein",   "function": "Pyochelin siderophore biosynthesis",           "category": "iron_metabolism"},
    "bfrA":  {"protein": "Bacterioferritin A",                "function": "Iron storage (redox-linked)",                  "category": "iron_metabolism"},
    "bfrB":  {"protein": "Bacterioferritin B",                "function": "Iron storage (redox-linked)",                  "category": "iron_metabolism"},

    # Phenazines
    "phzA1": {"protein": "Phenazine biosynthesis protein A1", "function": "Phenazine biosynthesis (electron shuttle)",    "category": "phenazine"},
    "phzB1": {"protein": "Phenazine biosynthesis protein B1", "function": "Phenazine biosynthesis",                       "category": "phenazine"},
    "phzM":  {"protein": "Phenazine methyltransferase",       "function": "Pyocyanin biosynthesis",                       "category": "phenazine"},
    "phzS":  {"protein": "Phenazine hydroxylase",             "function": "Pyocyanin biosynthesis",                       "category": "phenazine"},
    "phzH":  {"protein": "Phenazine amide synthetase",        "function": "Phenazine-1-carboxamide biosynthesis",         "category": "phenazine"},

    # Efflux pumps (PMF-driven = redox-linked)
    "mexA":  {"protein": "MexA (MexAB-OprM efflux)",          "function": "Drug efflux membrane fusion protein",         "category": "efflux"},
    "mexB":  {"protein": "MexB (MexAB-OprM efflux)",          "function": "RND pump — multidrug efflux",                 "category": "efflux"},
    "oprM":  {"protein": "OprM (MexAB-OprM efflux)",          "function": "Outer membrane channel",                      "category": "efflux"},
    "mexC":  {"protein": "MexC (MexCD-OprJ efflux)",          "function": "RND pump membrane fusion",                    "category": "efflux"},
    "mexD":  {"protein": "MexD (MexCD-OprJ efflux)",          "function": "RND pump — fluoroquinolone efflux",           "category": "efflux"},
    "mexE":  {"protein": "MexE (MexEF-OprN efflux)",          "function": "RND pump membrane fusion",                    "category": "efflux"},
    "mexF":  {"protein": "MexF (MexEF-OprN efflux)",          "function": "RND pump — redox-regulated efflux",           "category": "efflux"},
    "mexX":  {"protein": "MexX (MexXY-OprM efflux)",          "function": "RND pump — aminoglycoside efflux",            "category": "efflux"},
    "mexY":  {"protein": "MexY (MexXY-OprM efflux)",          "function": "RND pump",                                    "category": "efflux"},

    # Quorum sensing (redox-linked)
    "lasR":  {"protein": "LasR",                              "function": "3OC12-HSL-sensing QS regulator",              "category": "quorum_sensing"},
    "lasI":  {"protein": "LasI",                              "function": "3OC12-HSL synthase",                          "category": "quorum_sensing"},
    "rhlR":  {"protein": "RhlR",                              "function": "C4-HSL-sensing QS regulator",                 "category": "quorum_sensing"},
    "rhlI":  {"protein": "RhlI",                              "function": "C4-HSL synthase",                             "category": "quorum_sensing"},
    "mvfR":  {"protein": "MvfR (MvfR/PqsR)",                 "function": "PQS QS regulator — redox-linked virulence",   "category": "quorum_sensing"},
    "pqsA":  {"protein": "PqsA",                              "function": "PQS biosynthesis — CoA ligase",               "category": "quorum_sensing"},
}

# ---------------------------------------------------------------------------
# Ligands to detect in text
# ---------------------------------------------------------------------------
KNOWN_LIGANDS = [
    "NADH", "NADPH", "FAD", "FMN", "heme", "haem", "ubiquinone", "menaquinone",
    "pyocyanin", "pyoverdine", "pyochelin", "phenazine", "glutathione", "thioredoxin",
    "ferredoxin", "coenzyme A", "CoA", "ATP", "GTP", "iron", "Fe2+", "Fe3+",
    "zinc", "copper", "molybdenum", "flavin", "quinone",
]

# Regex: protein names that look like PA gene identifiers or conventional naming
# Matches: PA0001-style locus tags, or gene names 3-6 chars starting with lowercase
LOCUS_TAG_PATTERN = re.compile(r'\bPA\d{4}\b')
GENE_NAME_PATTERN = re.compile(
    r'\b(?:PA|pa)[a-zA-Z0-9]{2,5}\b|'      # PA-prefixed identifiers
    r'\b[a-z]{3}[A-Z0-9]{0,3}\b'           # conventional gene names like trxA, nuoB
)


# ---------------------------------------------------------------------------
# Text extraction helpers
# ---------------------------------------------------------------------------
def get_text(paper: dict) -> str:
    """Combine title and abstract into a single searchable string."""
    return f"{paper.get('title', '')} {paper.get('abstract', '')}".lower()


def extract_genes_from_text(text: str) -> set:
    """Find known gene names in text (case-insensitive)."""
    found = set()
    for gene in KNOWN_GENES:
        # Match whole word, case-insensitive
        if re.search(rf'\b{re.escape(gene)}\b', text, re.IGNORECASE):
            found.add(gene)
    return found


def extract_ligands_from_text(text: str) -> set:
    """Find known ligand names in text (case-insensitive)."""
    found = set()
    for ligand in KNOWN_LIGANDS:
        if re.search(rf'\b{re.escape(ligand)}\b', text, re.IGNORECASE):
            found.add(ligand)
    return found


def extract_novel_locus_tags(text: str) -> set:
    """Extract PA#### locus tags not in the known-gene list."""
    found = set(LOCUS_TAG_PATTERN.findall(text.upper()))
    return found


# ---------------------------------------------------------------------------
# Main extraction logic
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
    """Deduplicate combined list by PMID then DOI."""
    seen_pmids = set()
    seen_dois = set()
    unique = []
    for p in papers:
        pmid = p.get("pmid", "").strip()
        doi = p.get("doi", "").strip().lower()
        if pmid and pmid in seen_pmids:
            continue
        if doi and doi in seen_dois:
            continue
        unique.append(p)
        if pmid:
            seen_pmids.add(pmid)
        if doi:
            seen_dois.add(doi)
    return unique


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading raw paper files...")
    papers = load_papers(PUBMED_RAW, EPMC_RAW)
    papers = deduplicate_papers(papers)
    print(f"Total unique papers for extraction: {len(papers):,}\n")

    if len(papers) < 1000:
        print(
            f"WARNING: Only {len(papers)} papers loaded. "
            "Run pubmed_fetch.py first and ensure it completed successfully."
        )

    # Accumulators
    # gene_data[gene] = { "count": int, "papers": [ {title, pmid, year} ], "ligands_co_mentioned": Counter }
    gene_data = defaultdict(lambda: {
        "count": 0,
        "papers": [],
        "ligands_co_mentioned": defaultdict(int),
        "protein": "",
        "function": "",
        "category": "",
    })

    ligand_data = defaultdict(lambda: {"count": 0, "papers": []})

    novel_loci = defaultdict(lambda: {"count": 0, "papers": []})

    for paper in papers:
        text = get_text(paper)
        title = paper.get("title", "No title")
        pmid = paper.get("pmid", "")
        year = paper.get("year", "")
        paper_ref = {"title": title, "pmid": pmid, "year": year}

        genes_found = extract_genes_from_text(text)
        ligands_found = extract_ligands_from_text(text)
        loci_found = extract_novel_locus_tags(text)

        for gene in genes_found:
            gd = gene_data[gene]
            gd["count"] += 1
            gd["protein"] = KNOWN_GENES[gene]["protein"]
            gd["function"] = KNOWN_GENES[gene]["function"]
            gd["category"] = KNOWN_GENES[gene]["category"]
            if len(gd["papers"]) < 5:   # store up to 5 sample titles
                gd["papers"].append(paper_ref)
            for ligand in ligands_found:
                gd["ligands_co_mentioned"][ligand] += 1

        for ligand in ligands_found:
            ligand_data[ligand]["count"] += 1
            if len(ligand_data[ligand]["papers"]) < 3:
                ligand_data[ligand]["papers"].append(paper_ref)

        for locus in loci_found:
            novel_loci[locus]["count"] += 1
            if len(novel_loci[locus]["papers"]) < 3:
                novel_loci[locus]["papers"].append(paper_ref)

    # Convert defaultdicts to plain dicts for JSON serialisation
    gene_output = {}
    for gene, data in sorted(gene_data.items(), key=lambda x: -x[1]["count"]):
        gene_output[gene] = {
            "protein": data["protein"],
            "function": data["function"],
            "category": data["category"],
            "paper_count": data["count"],
            "sample_papers": data["papers"],
            "ligands_co_mentioned": dict(
                sorted(data["ligands_co_mentioned"].items(), key=lambda x: -x[1])
            ),
        }

    ligand_output = {
        k: v for k, v in sorted(ligand_data.items(), key=lambda x: -x[1]["count"])
    }

    novel_output = {
        k: v for k, v in sorted(novel_loci.items(), key=lambda x: -x[1]["count"])
        if k not in {g.upper() for g in KNOWN_GENES}
    }

    result = {
        "total_papers_processed": len(papers),
        "genes": gene_output,
        "ligands": ligand_output,
        "novel_locus_tags": novel_output,
    }

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)

    print(f"Extraction complete.")
    print(f"  Known genes detected:   {len(gene_output)}")
    print(f"  Ligands detected:       {len(ligand_output)}")
    print(f"  Novel locus tags found: {len(novel_output)}")
    print(f"\nTop 10 genes by paper count:")
    for i, (gene, data) in enumerate(list(gene_output.items())[:10], 1):
        print(f"  {i:2}. {gene:<8} — {data['paper_count']:>4} papers — {data['protein']}")
    print(f"\nSaved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
