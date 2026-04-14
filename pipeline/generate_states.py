"""
generate_states.py  —  Phase 1, Step 2
Generates ORCA geometry optimisation input files for each molecular state
(neutral, radical, reduced) and runs ORCA on each.

Reads charge and multiplicity from proteins/{protein_id}/config.json.
Input geometry from qm_regions/{protein_id}.xyz.

Output per state:
    results/{protein_id}/{state}/geometry_opt/{protein_id}_{state}_opt.inp
    results/{protein_id}/{state}/geometry_opt/{protein_id}_{state}_opt.xyz  (after ORCA run)

Usage:
    python generate_states.py --protein katA_heme
    python generate_states.py --protein katA_heme --states neutral radical
    python generate_states.py --protein katA_heme --no-run   # generate inputs without running ORCA
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

PROTEINS_DIR = Path(__file__).parent.parent / "proteins"
QM_DIR       = Path(__file__).parent.parent / "qm_regions"
RESULTS_DIR  = Path(__file__).parent.parent / "results"

ALL_STATES = ["neutral", "radical", "reduced"]

ORCA_TEMPLATE = """\
# ORCA geometry optimisation — {protein_id} | state: {state}
! B3LYP def2-TZVP Opt

%pal
  nprocs 4
end

%maxcore 2048

* xyzfile {charge} {mult} {xyz_path}
"""


def load_config(protein_id: str) -> dict:
    config_path = PROTEINS_DIR / protein_id / "config.json"
    if not config_path.exists():
        print(f"ERROR: config.json not found at {config_path}")
        sys.exit(1)
    with open(config_path) as f:
        config = json.load(f)
    return config


def get_charge_mult(config: dict, state: str):
    states = config.get("states", {})
    if state not in states:
        print(f"ERROR: State '{state}' not found in config.json.")
        print(f"       Available states: {list(states.keys())}")
        sys.exit(1)
    entry = states[state]
    charge = entry.get("charge")
    mult   = entry.get("multiplicity")
    if charge is None or mult is None:
        print(f"ERROR: 'charge' or 'multiplicity' missing for state '{state}' in config.json")
        sys.exit(1)
    return int(charge), int(mult)


def get_orca_path() -> str:
    """Return full absolute path to ORCA executable.
    ORCA requires full pathname for parallel (nprocs > 1) runs."""
    import shutil
    path = shutil.which("orca")
    if not path:
        print("ERROR: 'orca' executable not found in PATH.")
        print("       Add your ORCA installation directory to PATH and retry.")
        sys.exit(1)
    return path


def run_orca(inp_path: Path) -> bool:
    out_path = inp_path.with_suffix(".out")
    orca_exe = get_orca_path()
    print(f"    Running ORCA on {inp_path.name} ...")
    print(f"    ORCA path: {orca_exe}")
    try:
        result = subprocess.run(
            [orca_exe, str(inp_path.resolve())],
            stdout=open(out_path, "w"),
            stderr=subprocess.STDOUT,
            cwd=inp_path.parent
        )
        if result.returncode != 0:
            print(f"    WARNING: ORCA exited with code {result.returncode}")
            return False
        print(f"    Completed. Output: {out_path}")
        return True
    except Exception as e:
        print(f"    ERROR running ORCA: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Generate ORCA geometry optimisation inputs for each molecular state."
    )
    parser.add_argument("--protein", required=True,
                        help="Protein-ligand pair ID (e.g. katA_heme)")
    parser.add_argument("--states", nargs="+", default=ALL_STATES,
                        choices=ALL_STATES,
                        help="States to process (default: all three)")
    parser.add_argument("--no-run", action="store_true",
                        help="Generate input files without running ORCA")
    args = parser.parse_args()

    config   = load_config(args.protein)
    xyz_path = QM_DIR / f"{args.protein}.xyz"

    if not xyz_path.exists():
        print(f"ERROR: {xyz_path} not found. Run extract_qm.py first.")
        sys.exit(1)

    print(f"\nGenerating geometry optimisation inputs for [{args.protein}]")
    print(f"  States: {args.states}")
    print(f"  Input geometry: {xyz_path}\n")

    results = {}

    for state in args.states:
        charge, mult = get_charge_mult(config, state)
        out_dir = RESULTS_DIR / args.protein / state / "geometry_opt"
        out_dir.mkdir(parents=True, exist_ok=True)

        inp_name = f"{args.protein}_{state}_opt.inp"
        inp_path = out_dir / inp_name

        # Write ORCA input — use absolute path to xyz so ORCA finds it
        inp_content = ORCA_TEMPLATE.format(
            protein_id=args.protein,
            state=state,
            charge=charge,
            mult=mult,
            xyz_path=str(xyz_path.resolve())
        )
        inp_path.write_text(inp_content)
        print(f"  [{state}] charge={charge} mult={mult}")
        print(f"    Input:  {inp_path}")

        if args.no_run:
            print(f"    Skipping ORCA run (--no-run flag set)")
            results[state] = "skipped"
            continue

        ok = run_orca(inp_path)
        results[state] = "ok" if ok else "failed"

        # ORCA writes the final optimised geometry to {input_stem}.xyz in the
        # same folder as the input. Our input is {protein_id}_{state}_opt.inp
        # so ORCA produces {protein_id}_{state}_opt.xyz directly — no rename needed.
        opt_xyz = out_dir / f"{args.protein}_{state}_opt.xyz"
        if ok:
            if opt_xyz.exists():
                print(f"    Optimised geometry: {opt_xyz}")
            else:
                # Scan for any .xyz that is not a trajectory file
                xyz_candidates = [
                    f for f in out_dir.glob("*.xyz")
                    if not f.name.endswith("_trj.xyz")
                ]
                if xyz_candidates:
                    xyz_candidates[0].rename(opt_xyz)
                    print(f"    Optimised geometry (renamed): {opt_xyz}")
                else:
                    print(f"    WARNING: No optimised .xyz found in {out_dir}")
                    print(f"    ORCA may have failed silently. Check the .out file.")
                    results[state] = "failed"

    print(f"\n{'='*50}")
    print("Geometry optimisation summary:")
    for state, status in results.items():
        print(f"  {state:<10} {status}")
    print(f"{'='*50}\n")

    if not args.no_run and any(v == "failed" for v in results.values()):
        print("Some states failed. Check .out files in results/{protein}/*/geometry_opt/")
        sys.exit(1)


if __name__ == "__main__":
    main()
