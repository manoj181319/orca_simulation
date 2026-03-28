"""
report_targets.py  —  Phase 0, Step 4
Generates a formatted Excel report for Shiva Sir with three sheets:
    1. Redox Protein Targets  — ranked table with gene, protein, paper count,
                                 category, score, known ligands, UniProt link
    2. Ligand Mentions        — frequency table of all ligands found
    3. Recommended Pairs      — top protein-ligand pairs with 'Needs Docking?' column

Input:
    literature_mining/extracted/scored_targets.json

Output:
    literature_mining/redox_amr_targets_Pa.xlsx

Usage:
    python report_targets.py
"""

import json
from pathlib import Path
import openpyxl
from openpyxl.styles import (
    Font, PatternFill, Alignment, Border, Side
)
from openpyxl.utils import get_column_letter

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
INPUT_FILE = Path(__file__).parent / "extracted" / "scored_targets.json"
OUTPUT_FILE = Path(__file__).parent / "redox_amr_targets_Pa.xlsx"

# ---------------------------------------------------------------------------
# UniProt base URL for P. aeruginosa PAO1 gene searches
# ---------------------------------------------------------------------------
UNIPROT_BASE = "https://www.uniprot.org/uniprotkb?query={gene}+pseudomonas+aeruginosa&reviewed=true"

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
HEADER_BG      = "1F4E79"   # Dark blue
HEADER_FG      = "FFFFFF"   # White
SUBHEADER_BG   = "2E75B6"   # Medium blue
ALT_ROW_BG     = "DEEAF1"   # Light blue
HIGHLIGHT_BG   = "E2EFDA"   # Light green (high score rows)
WARN_BG        = "FFF2CC"   # Light yellow (needs docking)
BORDER_COLOR   = "BFBFBF"

# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------
def header_font(size=11):
    return Font(name="Calibri", bold=True, color="FFFFFF", size=size)

def body_font(size=10):
    return Font(name="Calibri", size=size)

def header_fill(color=HEADER_BG):
    return PatternFill("solid", fgColor=color)

def row_fill(color):
    return PatternFill("solid", fgColor=color)

def center():
    return Alignment(horizontal="center", vertical="center", wrap_text=True)

def left():
    return Alignment(horizontal="left", vertical="center", wrap_text=True)

def thin_border():
    s = Side(style="thin", color=BORDER_COLOR)
    return Border(left=s, right=s, top=s, bottom=s)

def set_col_widths(ws, widths: dict):
    for col_letter, width in widths.items():
        ws.column_dimensions[col_letter].width = width

def apply_header_row(ws, row_num, values, bg=HEADER_BG):
    for col, val in enumerate(values, 1):
        cell = ws.cell(row=row_num, column=col, value=val)
        cell.font = header_font()
        cell.fill = header_fill(bg)
        cell.alignment = center()
        cell.border = thin_border()

def apply_body_row(ws, row_num, values, bg=None):
    for col, val in enumerate(values, 1):
        cell = ws.cell(row=row_num, column=col, value=val)
        cell.font = body_font()
        if bg:
            cell.fill = row_fill(bg)
        cell.alignment = left()
        cell.border = thin_border()

def add_hyperlink(ws, row, col, url, display_text):
    cell = ws.cell(row=row, column=col)
    cell.value = display_text
    cell.hyperlink = url
    cell.font = Font(name="Calibri", size=10, color="0563C1", underline="single")
    cell.alignment = left()
    cell.border = thin_border()

def freeze_and_filter(ws, freeze_cell="A2"):
    ws.freeze_panes = freeze_cell
    ws.auto_filter.ref = ws.dimensions


# ---------------------------------------------------------------------------
# Sheet 1: Redox Protein Targets
# ---------------------------------------------------------------------------
def build_sheet1(wb, ranked_genes):
    ws = wb.active
    ws.title = "Redox Protein Targets"
    ws.sheet_view.showGridLines = False
    ws.row_dimensions[1].height = 30

    headers = [
        "Rank", "Gene", "Protein Name", "Function", "Category",
        "Papers", "Score", "Has PDB", "Top Ligand", "UniProt Search"
    ]
    apply_header_row(ws, 1, headers)

    for i, gene_data in enumerate(ranked_genes, 2):
        rank = gene_data["rank"]
        score = gene_data["final_score"]
        bg = HIGHLIGHT_BG if score >= 0.70 else (ALT_ROW_BG if i % 2 == 0 else None)

        row_values = [
            rank,
            gene_data["gene"],
            gene_data["protein"],
            gene_data["function"],
            gene_data["category"].replace("_", " ").title(),
            gene_data["paper_count"],
            round(score, 4),
            "Yes" if gene_data["has_pdb"] else "No",
            gene_data.get("top_ligand", ""),
            "",   # hyperlink added separately
        ]
        apply_body_row(ws, i, row_values, bg=bg)

        # Hyperlink in column J
        uniprot_url = UNIPROT_BASE.format(gene=gene_data["gene"])
        add_hyperlink(ws, i, 10, uniprot_url, f"Search {gene_data['gene']}")

        ws.row_dimensions[i].height = 22

    set_col_widths(ws, {
        "A": 6, "B": 9, "C": 32, "D": 45, "E": 18,
        "F": 8, "G": 8, "H": 9, "I": 16, "J": 18
    })
    freeze_and_filter(ws)


# ---------------------------------------------------------------------------
# Sheet 2: Ligand Mentions
# ---------------------------------------------------------------------------
def build_sheet2(wb, scored_data):
    ws = wb.create_sheet("Ligand Mentions")
    ws.sheet_view.showGridLines = False
    ws.row_dimensions[1].height = 28

    # Aggregate ligand counts from all genes
    ligand_totals = {}
    for gene_data in scored_data.get("ranked_genes", []):
        for ligand, count in gene_data.get("ligands_co_mentioned", {}).items():
            ligand_totals[ligand] = ligand_totals.get(ligand, 0) + count

    sorted_ligands = sorted(ligand_totals.items(), key=lambda x: -x[1])

    headers = ["Rank", "Ligand / Cofactor", "Co-mention Count", "Simulatable?"]
    apply_header_row(ws, 1, headers, bg=SUBHEADER_BG)

    SIMULATABLE = {
        "NADH", "NADPH", "FAD", "FMN", "heme", "haem", "ubiquinone",
        "pyocyanin", "phenazine", "glutathione", "ferredoxin",
        "coenzyme A", "CoA", "flavin", "quinone",
    }

    for i, (ligand, count) in enumerate(sorted_ligands, 2):
        sim = "Yes" if ligand in SIMULATABLE else "No"
        bg = HIGHLIGHT_BG if sim == "Yes" else (ALT_ROW_BG if i % 2 == 0 else None)
        apply_body_row(ws, i, [i - 1, ligand, count, sim], bg=bg)
        ws.row_dimensions[i].height = 20

    set_col_widths(ws, {"A": 6, "B": 24, "C": 18, "D": 14})
    freeze_and_filter(ws)


# ---------------------------------------------------------------------------
# Sheet 3: Recommended Pairs
# ---------------------------------------------------------------------------
def build_sheet3(wb, recommended_pairs):
    ws = wb.create_sheet("Recommended Pairs")
    ws.sheet_view.showGridLines = False
    ws.row_dimensions[1].height = 28

    headers = [
        "Priority", "Gene", "Protein Name", "Suggested Ligand",
        "Score", "Has PDB Structure", "Needs Docking?", "Notes"
    ]
    apply_header_row(ws, 1, headers, bg=SUBHEADER_BG)

    for i, pair in enumerate(recommended_pairs, 2):
        needs_docking = "YES — request from Shiva Sir" if pair["needs_docking"] else "No — pull from PDB"
        bg = WARN_BG if pair["needs_docking"] else (ALT_ROW_BG if i % 2 == 0 else None)

        row_values = [
            pair["rank"],
            pair["gene"],
            pair["protein"],
            pair["suggested_ligand"],
            round(pair["score"], 4),
            "Yes" if pair["has_pdb"] else "No",
            needs_docking,
            "",   # Notes column — left blank for coordinator to fill
        ]
        apply_body_row(ws, i, row_values, bg=bg)
        ws.row_dimensions[i].height = 22

    set_col_widths(ws, {
        "A": 8, "B": 9, "C": 32, "D": 18,
        "E": 8, "F": 18, "G": 28, "H": 30
    })
    freeze_and_filter(ws)


# ---------------------------------------------------------------------------
# Title / metadata sheet (injected as first sheet)
# ---------------------------------------------------------------------------
def build_cover(wb):
    ws = wb.create_sheet("About", 0)
    ws.sheet_view.showGridLines = False

    meta = [
        ("Project",     "UV-Induced Excitation in Bacterial Proteins"),
        ("Organism",    "Pseudomonas aeruginosa"),
        ("Focus",       "Redox-linked antimicrobial resistance"),
        ("Phase",       "Phase 0 — Literature Mining"),
        ("Coordinator", "Shiva Sir"),
        ("Generated by","orca_simulation pipeline (report_targets.py)"),
        ("Description", (
            "This report identifies protein and gene candidates from "
            "P. aeruginosa involved in redox biology and antimicrobial "
            "resistance. Candidates are ranked by citation frequency, "
            "PDB structure availability, and co-mention with simulatable ligands."
        )),
    ]

    ws.column_dimensions["A"].width = 16
    ws.column_dimensions["B"].width = 70

    for row, (key, val) in enumerate(meta, 2):
        cell_key = ws.cell(row=row, column=1, value=key)
        cell_key.font = Font(name="Calibri", bold=True, size=11)
        cell_key.alignment = left()

        cell_val = ws.cell(row=row, column=2, value=val)
        cell_val.font = Font(name="Calibri", size=11)
        cell_val.alignment = left()
        ws.row_dimensions[row].height = 22

    # Title
    ws.row_dimensions[1].height = 36
    title_cell = ws.cell(row=1, column=1, value="Redox AMR Targets — Pseudomonas aeruginosa")
    title_cell.font = Font(name="Calibri", bold=True, size=14, color=HEADER_BG)
    title_cell.alignment = left()
    ws.merge_cells("A1:B1")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if not INPUT_FILE.exists():
        print(f"ERROR: {INPUT_FILE} not found. Run score_targets.py first.")
        return

    with open(INPUT_FILE, encoding="utf-8") as f:
        scored_data = json.load(f)

    ranked_genes = scored_data.get("ranked_genes", [])
    recommended_pairs = scored_data.get("recommended_pairs", [])

    if not ranked_genes:
        print("No ranked genes found in scored_targets.json. Exiting.")
        return

    print(f"Building Excel report for {len(ranked_genes)} genes "
          f"and {len(recommended_pairs)} recommended pairs...")

    wb = openpyxl.Workbook()

    # Sheet order: About, Redox Protein Targets, Ligand Mentions, Recommended Pairs
    build_cover(wb)           # adds "About" at index 0
    build_sheet1(wb, ranked_genes)        # active sheet (index 1 after cover)
    build_sheet2(wb, scored_data)
    build_sheet3(wb, recommended_pairs)

    # Remove default empty sheet if it exists
    if "Sheet" in wb.sheetnames:
        del wb["Sheet"]

    # Set "Redox Protein Targets" as the active sheet on open
    wb.active = wb["Redox Protein Targets"]

    wb.save(OUTPUT_FILE)

    print(f"\nReport saved to: {OUTPUT_FILE}")
    print(f"Sheets: {wb.sheetnames}")
    print(f"\nNext step: send {OUTPUT_FILE.name} to Shiva Sir for target confirmation.")


if __name__ == "__main__":
    main()
