"""
fetch_structures.py  —  Phase 1, Pre-Step
Downloads crystal structures from RCSB PDB for each protein-ligand pair.
Reads pdb_id from proteins/{protein_id}/config.json and saves the structure
to proteins/{protein_id}/docked.pdb (the path expected by all downstream scripts).

Usage:
    python fetch_structures.py                        # fetch all proteins in proteins/
    python fetch_structures.py --protein katA_heme    # fetch one specific protein
    python fetch_structures.py --force                # re-download even if file exists
"""

import argparse
import json
import sys
import time
import requests
from pathlib import Path

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PROTEINS_DIR = Path(__file__).parent.parent / "proteins"


def load_config(protein_id: str) -> dict:
    config_path = PROTEINS_DIR / protein_id / "config.json"
    if not config_path.exists():
        print(f"  ERROR: config.json not found at {config_path}")
        sys.exit(1)
    with open(config_path) as f:
        config = json.load(f)
    if "pdb_id" not in config:
        print(f"  ERROR: 'pdb_id' field missing in {config_path}")
        sys.exit(1)
    return config


def fetch_pdb(pdb_id: str, out_path: Path) -> bool:
    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    print(f"  Downloading {pdb_id.upper()} from RCSB...")
    try:
        resp = requests.get(url, timeout=60)
        if resp.status_code == 404:
            print(f"  ERROR: PDB ID {pdb_id.upper()} not found on RCSB.")
            return False
        resp.raise_for_status()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(resp.text, encoding="utf-8")
        size_kb = out_path.stat().st_size / 1024
        print(f"  Saved to {out_path}  ({size_kb:.1f} KB)")
        return True
    except Exception as e:
        print(f"  ERROR downloading {pdb_id}: {e}")
        return False


def get_all_protein_ids() -> list:
    if not PROTEINS_DIR.exists():
        print(f"ERROR: proteins/ directory not found at {PROTEINS_DIR}")
        sys.exit(1)
    ids = [
        d.name for d in PROTEINS_DIR.iterdir()
        if d.is_dir() and (d / "config.json").exists()
    ]
    if not ids:
        print("ERROR: No protein folders with config.json found in proteins/")
        sys.exit(1)
    return sorted(ids)


def main():
    parser = argparse.ArgumentParser(
        description="Download PDB crystal structures from RCSB for simulation."
    )
    parser.add_argument(
        "--protein",
        help="Protein-ligand pair ID (e.g. katA_heme). If omitted, fetches all."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if docked.pdb already exists."
    )
    args = parser.parse_args()

    protein_ids = [args.protein] if args.protein else get_all_protein_ids()

    print(f"\nFetching structures for {len(protein_ids)} protein(s)...\n")

    success, skipped, failed = [], [], []

    for protein_id in protein_ids:
        print(f"[{protein_id}]")
        config = load_config(protein_id)
        pdb_id = config["pdb_id"]
        out_path = PROTEINS_DIR / protein_id / "docked.pdb"

        if out_path.exists() and not args.force:
            print(f"  SKIP: docked.pdb already exists. Use --force to re-download.")
            skipped.append(protein_id)
            continue

        ok = fetch_pdb(pdb_id, out_path)
        if ok:
            success.append(protein_id)
        else:
            failed.append(protein_id)
        time.sleep(0.5)

    print(f"\n{'='*50}")
    print(f"Done. Success: {len(success)}  Skipped: {len(skipped)}  Failed: {len(failed)}")
    if failed:
        print(f"Failed proteins: {failed}")
        print("Check PDB IDs in their config.json files.")
    print(f"{'='*50}\n")


if __name__ == "__main__":
    main()
