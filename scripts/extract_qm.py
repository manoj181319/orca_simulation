"""
extract_qm.py
=============
Script 1 of 4 — QM Region Extraction

Pipeline position:
    00_raw_protein/docked_complex.pdb
            ↓
    [this script]
            ↓
    01_qm_region/qm_region.xyz

Method: Distance cutoff around ligand
    - Locates ligand atoms by residue name in PDB
    - Selects all atoms within CUTOFF angstroms of any ligand atom
    - Writes selected atoms to XYZ format

Usage:
    Run from orca_simulation/ root directory:
    python scripts/extract_qm.py
"""

import os
from Bio import PDB
from Bio.PDB import NeighborSearch

# =============================================================
# CONFIGURATION
# Edit these values before running
# =============================================================

LIGAND_NAME = "LIG"       # 3-letter residue name of your ligand as it
                           # appears in the PDB file (e.g. LIG, HEM, FAD)
                           # Open docked_complex.pdb and check HETATM lines

CUTOFF = 5.0               # Distance cutoff in angstroms
                           # All atoms within this distance from the ligand
                           # will be included in the QM region
                           # Typical range: 5.0 to 8.0 — confirm with guide

# =============================================================
# PATHS
# Relative to orca_simulation/ root — do not edit
# =============================================================

INPUT_PDB  = os.path.join("00_raw_protein", "docked_complex.pdb")
OUTPUT_DIR = os.path.join("01_qm_region")
OUTPUT_XYZ = os.path.join(OUTPUT_DIR, "qm_region.xyz")

# =============================================================
# HELPER FUNCTIONS
# =============================================================

def get_element_symbol(atom):
    """
    Returns clean element symbol from Biopython atom object.
    XYZ format requires element symbol (C, N, O, H etc.)
    not ATOM names (CA, CB, OG1 etc.)
    """
    element = atom.element.strip()
    if element:
        return element.capitalize()
    # Fallback if element field is empty in PDB
    name = atom.get_name().strip()
    for char in name:
        if char.isalpha():
            return char.capitalize()
    return "X"


def load_pdb(input_pdb):
    """
    Load PDB file using Biopython parser.
    Returns Biopython structure object.
    """
    print(f"\n[Step 1/4] Loading PDB file")
    print(f"           Path : {input_pdb}")

    if not os.path.exists(input_pdb):
        raise FileNotFoundError(
            f"\nERROR: PDB file not found at '{input_pdb}'\n"
            f"Make sure you are running this script from the orca_simulation/ directory.\n"
            f"Command: python scripts/extract_qm.py"
        )

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("docked_complex", input_pdb)

    total_atoms = len(list(structure.get_atoms()))
    print(f"           Total atoms in structure : {total_atoms}")

    return structure


def get_ligand_atoms(structure, ligand_name):
    """
    Find and return all atoms belonging to the specified ligand.
    Ligand is identified by its 3-letter residue name in the PDB.
    """
    print(f"\n[Step 2/4] Locating ligand")
    print(f"           Ligand name : {ligand_name}")

    ligand_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                if resname == ligand_name:
                    atoms_in_residue = list(residue.get_atoms())
                    ligand_atoms.extend(atoms_in_residue)
                    print(f"           Found in chain {chain.id}, "
                          f"residue {residue.get_id()[1]} "
                          f"— {len(atoms_in_residue)} atoms")

    if not ligand_atoms:
        raise ValueError(
            f"\nERROR: Ligand '{ligand_name}' not found in PDB file.\n"
            f"Steps to fix:\n"
            f"  1. Open 00_raw_protein/docked_complex.pdb\n"
            f"  2. Search for HETATM lines\n"
            f"  3. Column 4 contains the residue name — use that as LIGAND_NAME\n"
        )

    print(f"           Total ligand atoms : {len(ligand_atoms)}")
    return ligand_atoms


def select_qm_region(structure, ligand_atoms, cutoff):
    """
    Select all atoms within cutoff distance from any ligand atom.
    Uses Biopython NeighborSearch for efficient distance calculation.
    Ligand atoms are always included regardless of distance.
    """
    print(f"\n[Step 3/4] Selecting QM region")
    print(f"           Cutoff distance : {cutoff} Å")

    all_atoms = list(structure.get_atoms())
    neighbor_search = NeighborSearch(all_atoms)

    # Use a set to avoid duplicate atoms
    qm_atom_set = set()

    # Always include ligand atoms
    for atom in ligand_atoms:
        qm_atom_set.add(atom)

    # Find all atoms within cutoff of each ligand atom
    for lig_atom in ligand_atoms:
        nearby_atoms = neighbor_search.search(
            center = lig_atom.get_vector(),
            radius = cutoff,
            level  = "A"       # A = atom level search
        )
        for atom in nearby_atoms:
            qm_atom_set.add(atom)

    qm_atoms = list(qm_atom_set)
    print(f"           Ligand atoms            : {len(ligand_atoms)}")
    print(f"           Surrounding atoms found : {len(qm_atoms) - len(ligand_atoms)}")
    print(f"           Total QM region atoms   : {len(qm_atoms)}")

    return qm_atoms


def write_xyz(qm_atoms, output_xyz, input_pdb, cutoff):
    """
    Write QM region atoms to XYZ format file.

    XYZ format (required by ORCA):
        Line 1  : number of atoms
        Line 2  : comment line
        Line 3+ : element  x  y  z
    """
    print(f"\n[Step 4/4] Writing XYZ file")
    print(f"           Output : {output_xyz}")

    os.makedirs(os.path.dirname(output_xyz), exist_ok=True)

    n_atoms = len(qm_atoms)
    comment = (
        f"QM region | source: {input_pdb} | "
        f"ligand: {LIGAND_NAME} | cutoff: {cutoff}A | "
        f"atoms: {n_atoms}"
    )

    with open(output_xyz, "w") as f:
        f.write(f"{n_atoms}\n")
        f.write(f"{comment}\n")
        for atom in qm_atoms:
            element = get_element_symbol(atom)
            x, y, z = atom.get_vector().get_array()
            f.write(f"{element:<4s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")

    print(f"           {n_atoms} atoms written successfully")


# =============================================================
# MAIN
# =============================================================

def main():
    print("=" * 55)
    print("  extract_qm.py — QM Region Extraction")
    print("  Script 1 of 4")
    print("=" * 55)
    print(f"  Input  : {INPUT_PDB}")
    print(f"  Output : {OUTPUT_XYZ}")
    print(f"  Ligand : {LIGAND_NAME}")
    print(f"  Cutoff : {CUTOFF} Å")

    structure    = load_pdb(INPUT_PDB)
    ligand_atoms = get_ligand_atoms(structure, LIGAND_NAME)
    qm_atoms     = select_qm_region(structure, ligand_atoms, CUTOFF)
    write_xyz(qm_atoms, OUTPUT_XYZ, INPUT_PDB, CUTOFF)

    print("\n" + "=" * 55)
    print("  Extraction complete.")
    print(f"  QM region saved to : {OUTPUT_XYZ}")
    print("  Next step          : run scripts/generate_states.py")
    print("=" * 55)


if __name__ == "__main__":
    main()
