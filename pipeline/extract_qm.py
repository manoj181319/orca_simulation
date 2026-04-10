"""
extract_qm.py  —  Phase 1, Step 1
Extracts the QM region from a downloaded PDB crystal structure.

The QM region is defined as:
    - All atoms of the ligand (identified by ligand_residue in config.json)
    - All protein residues that have at least one atom within --cutoff angstroms
      of any ligand atom (Biopython NeighborSearch)

Writes the QM region as an .xyz file for use by generate_states.py.

Usage:
    python extract_qm.py --protein katA_heme
    python extract_qm.py --protein katA_heme --cutoff 6.0
"""

import argparse
import json
import sys
from pathlib import Path

try:
    from Bio.PDB import PDBParser, NeighborSearch, Selection
except ImportError:
    print("ERROR: Biopython not installed. Run: pip install biopython")
    sys.exit(1)

PROTEINS_DIR = Path(__file__).parent.parent / "proteins"
QM_DIR       = Path(__file__).parent.parent / "qm_regions"

ELEMENT_FALLBACK = {
    "CA": "C", "CB": "C", "CG": "C", "CD": "C", "CE": "C", "CZ": "C",
    "NA": "N", "NB": "N", "ND": "N", "NE": "N", "NH": "N", "NZ": "N",
    "OA": "O", "OB": "O", "OD": "O", "OE": "O", "OG": "O", "OH": "O",
    "SA": "S", "SB": "S", "SD": "S", "SG": "S",
    "FE": "Fe", "MN": "Mn", "ZN": "Zn", "CU": "Cu", "MG": "Mg",
}


def load_config(protein_id: str) -> dict:
    config_path = PROTEINS_DIR / protein_id / "config.json"
    if not config_path.exists():
        print(f"ERROR: config.json not found at {config_path}")
        sys.exit(1)
    with open(config_path) as f:
        config = json.load(f)
    for field in ("pdb_id", "ligand_residue"):
        if field not in config:
            print(f"ERROR: '{field}' missing in {config_path}")
            sys.exit(1)
    return config


def get_element(atom) -> str:
    """Get element symbol, falling back to first char of atom name."""
    el = atom.element
    if el and el.strip() and el.strip() != "X":
        return el.strip().capitalize()
    name = atom.get_name().strip()
    return ELEMENT_FALLBACK.get(name[:2].upper(),
           ELEMENT_FALLBACK.get(name[:1].upper(), name[0].capitalize()))


def extract_qm_region(pdb_path: Path, ligand_resname: str, cutoff: float):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]

    # Collect all ligand atoms (HETATM with matching residue name)
    ligand_atoms = []
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip() == ligand_resname.strip().upper():
                ligand_atoms.extend(residue.get_atoms())

    if not ligand_atoms:
        print(f"ERROR: Ligand residue '{ligand_resname}' not found in {pdb_path}")
        print("       Check the HETATM records in the PDB file for the correct residue name.")
        sys.exit(1)

    print(f"  Ligand '{ligand_resname}': {len(ligand_atoms)} atoms found")

    # Neighbor search across all protein atoms
    all_atoms = list(model.get_atoms())
    ns = NeighborSearch(all_atoms)

    # Find residues within cutoff of any ligand atom
    nearby_residues = set()
    for lat in ligand_atoms:
        hits = ns.search(lat.get_vector().get_array(), cutoff, level="R")
        for res in hits:
            nearby_residues.add(res)

    print(f"  Protein residues within {cutoff} Å: {len(nearby_residues)}")

    # Collect all atoms: ligand + nearby residues (excluding water HOH/WAT)
    qm_atoms = []
    for atom in ligand_atoms:
        coord = atom.get_vector().get_array()
        el = get_element(atom)
        qm_atoms.append((el, coord[0], coord[1], coord[2]))

    for residue in nearby_residues:
        resname = residue.get_resname().strip()
        if resname in ("HOH", "WAT", "H2O"):
            continue
        if resname == ligand_resname.strip().upper():
            continue   # already added above
        for atom in residue.get_atoms():
            coord = atom.get_vector().get_array()
            el = get_element(atom)
            qm_atoms.append((el, coord[0], coord[1], coord[2]))

    return qm_atoms


def write_xyz(atoms: list, out_path: Path, comment: str = ""):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for el, x, y, z in atoms:
            f.write(f"{el:<4}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract QM region (ligand + neighboring residues) from PDB structure."
    )
    parser.add_argument("--protein", required=True,
                        help="Protein-ligand pair ID (e.g. katA_heme)")
    parser.add_argument("--cutoff", type=float, default=5.0,
                        help="Distance cutoff in angstroms (default: 5.0)")
    args = parser.parse_args()

    config = load_config(args.protein)
    pdb_path = PROTEINS_DIR / args.protein / "docked.pdb"
    out_path = QM_DIR / f"{args.protein}.xyz"

    if not pdb_path.exists():
        print(f"ERROR: {pdb_path} not found. Run fetch_structures.py first.")
        sys.exit(1)

    print(f"\nExtracting QM region for [{args.protein}]")
    print(f"  PDB file:        {pdb_path}")
    print(f"  Ligand residue:  {config['ligand_residue']}")
    print(f"  Cutoff:          {args.cutoff} Å")

    atoms = extract_qm_region(pdb_path, config["ligand_residue"], args.cutoff)

    comment = (f"QM region: {args.protein} | ligand={config['ligand_residue']} "
               f"| cutoff={args.cutoff}A | {len(atoms)} atoms")
    write_xyz(atoms, out_path, comment)

    print(f"  Total QM atoms:  {len(atoms)}")
    print(f"  Written to:      {out_path}\n")


if __name__ == "__main__":
    main()
