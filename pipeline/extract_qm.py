"""
extract_qm.py  —  Phase 1, Step 1
Extracts the QM region from a PDB crystal structure.

The QM region is defined as:
    - All atoms of the ligand (HETATM residue matching ligand_residue in config.json)
      from the specified chain only
    - All protein residues in that same chain within --cutoff angstroms of the ligand
    - Water molecules excluded

IMPORTANT — always specify --chain for multimers.
Many PDB structures are tetramers, dimers etc. Without --chain the script captures
the ligand from every subunit, giving a QM region far too large for DFT.
Default chain is A.

Recommended settings for desktop DFT (i5, 16 GB):
    --chain A --cutoff 3.5    →  typically 80–150 atoms  (recommended)
    --chain A --cutoff 5.0    →  typically 200–350 atoms (slow, may fail)

Usage:
    python extract_qm.py --protein katA_heme --chain A --cutoff 3.5
    python extract_qm.py --protein katA_heme --list-chains
"""

import argparse
import json
import sys
from pathlib import Path

try:
    from Bio.PDB import PDBParser, NeighborSearch
except ImportError:
    print("ERROR: Biopython not installed. Run: pip install biopython")
    sys.exit(1)

PROTEINS_DIR = Path(__file__).parent.parent / "proteins"
QM_DIR       = Path(__file__).parent.parent / "qm_regions"

WARN_ATOMS   = 200
ABORT_ATOMS  = 500

SKIP_RESIDUES = {"HOH", "WAT", "H2O", "DOD"}

ELEMENT_FALLBACK = {
    "CA": "C",  "CB": "C",  "CG": "C",  "CD": "C",  "CE": "C",  "CZ": "C",
    "NA": "N",  "NB": "N",  "ND": "N",  "NE": "N",  "NH": "N",  "NZ": "N",
    "OA": "O",  "OB": "O",  "OD": "O",  "OE": "O",  "OG": "O",  "OH": "O",
    "SA": "S",  "SB": "S",  "SD": "S",  "SG": "S",
    "FE": "Fe", "MN": "Mn", "ZN": "Zn", "CU": "Cu", "MG": "Mg",
}


def load_config(protein_id):
    config_path = PROTEINS_DIR / protein_id / "config.json"
    if not config_path.exists():
        print(f"ERROR: config.json not found at {config_path}")
        print("       Run setup_proteins.py first.")
        sys.exit(1)
    with open(config_path) as f:
        config = json.load(f)
    if not config.get("ligand_residue", "").strip():
        print(f"ERROR: 'ligand_residue' is empty in config.json for [{protein_id}].")
        print("       Open docked.pdb, find HETATM records, note 3-letter code,")
        print("       update ligand_residue in config.json.")
        sys.exit(1)
    return config


def get_element(atom):
    el = atom.element
    if el and el.strip() and el.strip() not in ("X", ""):
        return el.strip().capitalize()
    name = atom.get_name().strip()
    return ELEMENT_FALLBACK.get(
        name[:2].upper(),
        ELEMENT_FALLBACK.get(name[:1].upper(), name[0].capitalize())
    )


def extract_qm_region(pdb_path, ligand_resname, chain_id, cutoff):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]

    chain_ids = [c.id for c in model]
    if chain_id not in chain_ids:
        print(f"ERROR: Chain '{chain_id}' not found in {pdb_path.name}")
        print(f"       Available chains: {chain_ids}")
        print(f"       Re-run with --list-chains to see what's available.")
        sys.exit(1)

    chain = model[chain_id]
    ligand_resname_upper = ligand_resname.strip().upper()

    # Collect ligand atoms from this chain only
    ligand_atoms = []
    for residue in chain:
        if residue.get_resname().strip() == ligand_resname_upper:
            ligand_atoms.extend(residue.get_atoms())

    if not ligand_atoms:
        # Report what HETATM residues ARE in this chain to help debug
        hetatm_found = sorted({
            r.get_resname().strip()
            for r in chain
            if r.id[0] != " " and r.get_resname().strip() not in SKIP_RESIDUES
        })
        print(f"ERROR: Ligand '{ligand_resname}' not found in chain {chain_id}")
        print(f"       HETATM residues in chain {chain_id}: {hetatm_found}")
        print(f"       Available chains: {chain_ids}")
        print(f"       Run --list-chains to see all chains and their ligands.")
        sys.exit(1)

    print(f"  Ligand '{ligand_resname}' in chain {chain_id}: {len(ligand_atoms)} atoms")

    # Neighbor search restricted to this chain only
    chain_atoms = list(chain.get_atoms())
    ns = NeighborSearch(chain_atoms)

    nearby_residues = set()
    for lat in ligand_atoms:
        hits = ns.search(lat.get_vector().get_array(), cutoff, level="R")
        for res in hits:
            nearby_residues.add(res)

    print(f"  Residues within {cutoff} A (chain {chain_id}): {len(nearby_residues)}")

    # Build QM atom list
    qm_atoms = []
    for atom in ligand_atoms:
        coord = atom.get_vector().get_array()
        qm_atoms.append((get_element(atom), coord[0], coord[1], coord[2]))

    for residue in sorted(nearby_residues, key=lambda r: r.id[1]):
        resname = residue.get_resname().strip()
        if resname in SKIP_RESIDUES or resname == ligand_resname_upper:
            continue
        for atom in residue.get_atoms():
            coord = atom.get_vector().get_array()
            qm_atoms.append((get_element(atom), coord[0], coord[1], coord[2]))

    return qm_atoms


def write_xyz(atoms, out_path, comment=""):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for el, x, y, z in atoms:
            f.write(f"{el:<4}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract QM region (ligand + neighboring residues, one chain) from PDB."
    )
    parser.add_argument("--protein", required=True,
                        help="Protein-ligand pair ID (e.g. katA_heme)")
    parser.add_argument("--chain", default="A",
                        help="PDB chain ID to extract from (default: A)")
    parser.add_argument("--cutoff", type=float, default=3.5,
                        help="Distance cutoff in angstroms (default: 3.5)")
    parser.add_argument("--list-chains", action="store_true",
                        help="List all chains and their HETATM residues then exit")
    args = parser.parse_args()

    config   = load_config(args.protein)
    pdb_path = PROTEINS_DIR / args.protein / "docked.pdb"

    if not pdb_path.exists():
        print(f"ERROR: {pdb_path} not found. Run fetch_structures.py first.")
        sys.exit(1)

    if args.list_chains:
        pdb_parser = PDBParser(QUIET=True)
        model = pdb_parser.get_structure("p", str(pdb_path))[0]
        print(f"\nChains in {pdb_path.name}:")
        for chain in model:
            hetatm = sorted({
                r.get_resname().strip()
                for r in chain
                if r.id[0] != " " and r.get_resname().strip() not in SKIP_RESIDUES
            })
            n_res = sum(1 for r in chain if r.id[0] == " ")
            print(f"  Chain {chain.id}  —  {n_res} residues  —  HETATM: {hetatm if hetatm else '(none)'}")
        print()
        return

    print(f"\nExtracting QM region for [{args.protein}]")
    print(f"  PDB:      {pdb_path.name}")
    print(f"  Chain:    {args.chain}")
    print(f"  Ligand:   {config['ligand_residue']}")
    print(f"  Cutoff:   {args.cutoff} A")

    atoms = extract_qm_region(
        pdb_path, config["ligand_residue"], args.chain, args.cutoff
    )

    n = len(atoms)
    print(f"\n  Total QM atoms: {n}")

    if n > ABORT_ATOMS:
        print(f"\n  ERROR: {n} atoms is too large for desktop DFT.")
        print(f"  ORCA will likely crash or take weeks on an i5 with 16 GB RAM.")
        print(f"  Fix: reduce --cutoff to 3.0 or 3.5 and re-run.")
        print(f"  The .xyz file has NOT been written.")
        sys.exit(1)
    elif n > WARN_ATOMS:
        print(f"  WARNING: {n} atoms is on the large side. Each ORCA job may take")
        print(f"  several hours. Consider --cutoff 3.0 for faster runs.")

    out_path = QM_DIR / f"{args.protein}.xyz"
    comment  = (f"QM region: {args.protein} | ligand={config['ligand_residue']} "
                f"| chain={args.chain} | cutoff={args.cutoff}A | {n} atoms")
    write_xyz(atoms, out_path, comment)
    print(f"  Written to: {out_path}\n")


if __name__ == "__main__":
    main()
