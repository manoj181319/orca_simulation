"""
generate_inputs.py  —  Phase 1, Step 3
Generates all ORCA input files for the full simulation matrix:
    3 states × 4 methods × 5 epsilon values = 60 inputs per protein

Each input file is placed in its own dedicated folder:
    results/{protein_id}/{state}/{method}/eps_{value}/
        {protein_id}_{state}_{method}_eps{value}.inp

Reads:
    - Optimised geometry from results/{protein_id}/{state}/geometry_opt/{protein_id}_{state}_opt.xyz
    - Charge and multiplicity from proteins/{protein_id}/config.json

Usage:
    python generate_inputs.py --protein katA_heme
    python generate_inputs.py --protein katA_heme --states neutral
    python generate_inputs.py --protein katA_heme --methods DFT TDDFT
    python generate_inputs.py --protein katA_heme --eps 10 20
"""

import argparse
import json
import sys
from pathlib import Path

PROTEINS_DIR = Path(__file__).parent.parent / "proteins"
RESULTS_DIR  = Path(__file__).parent.parent / "results"

ALL_STATES  = ["neutral", "radical", "reduced"]
ALL_METHODS = ["DFT", "CDFT", "SDFT", "TDDFT"]
ALL_EPSILON = [10, 20, 30, 40, 80]

# ─── ORCA input templates per method ────────────────────────────────────────

DFT_TEMPLATE = """\
# ORCA DFT single point — {protein_id} | {state} | eps={eps}
! B3LYP def2-TZVP

%pal
  nprocs 4
end

%maxcore 2048

%cpcm
  epsilon {eps}
end

* xyzfile {charge} {mult} {xyz_path}
"""

# NOTE ON CDFT: Before running CDFT calculations, you MUST edit each generated
# .inp file to specify correct atom indices for the constrained region.
# The placeholder "atoms_alpha 1 end" constrains only atom 1 — this is wrong
# for your system. Confirm the CDFT atom selection with Shiva Sir first.
CDFT_TEMPLATE = """# ORCA Constrained DFT single point — {protein_id} | {state} | eps={eps}
# WARNING: Edit %cdft atom indices before running. See NOTE ON CDFT above.
! B3LYP def2-TZVP

%pal
  nprocs 4
end

%maxcore 2048

%cpcm
  epsilon {eps}
end

%cdft
  charge_cons 1
  atoms_alpha 1 end
  atoms_beta  1 end
end

* xyzfile {charge} {mult} {xyz_path}
"""

SDFT_TEMPLATE = """\
# ORCA Spin-DFT single point — {protein_id} | {state} | eps={eps}
! UKS B3LYP def2-TZVP

%pal
  nprocs 4
end

%maxcore 2048

%cpcm
  epsilon {eps}
end

* xyzfile {charge} {mult} {xyz_path}
"""

TDDFT_TEMPLATE = """\
# ORCA TDDFT single point — {protein_id} | {state} | eps={eps}
! B3LYP def2-TZVP

%pal
  nprocs 4
end

%maxcore 2048

%cpcm
  epsilon {eps}
end

%tddft
  nroots 10
  maxdim  5
end

* xyzfile {charge} {mult} {xyz_path}
"""

METHOD_TEMPLATES = {
    "DFT":   DFT_TEMPLATE,
    "CDFT":  CDFT_TEMPLATE,
    "SDFT":  SDFT_TEMPLATE,
    "TDDFT": TDDFT_TEMPLATE,
}


def load_config(protein_id: str) -> dict:
    config_path = PROTEINS_DIR / protein_id / "config.json"
    if not config_path.exists():
        print(f"ERROR: config.json not found at {config_path}")
        sys.exit(1)
    with open(config_path) as f:
        return json.load(f)


def get_charge_mult(config: dict, state: str):
    states = config.get("states", {})
    if state not in states:
        print(f"ERROR: State '{state}' not found in config.json")
        sys.exit(1)
    entry = states[state]
    charge = entry.get("charge")
    mult   = entry.get("multiplicity")
    if charge is None or mult is None:
        print(f"ERROR: 'charge' or 'multiplicity' missing for state '{state}' in config.json")
        sys.exit(1)
    return int(charge), int(mult)


def validate_epsilon(eps_list: list) -> list:
    invalid = [e for e in eps_list if e not in ALL_EPSILON]
    if invalid:
        print(f"ERROR: Invalid epsilon value(s): {invalid}")
        print(f"       Valid values are: {ALL_EPSILON}")
        sys.exit(1)
    return eps_list


def main():
    parser = argparse.ArgumentParser(
        description="Generate all ORCA input files for the full simulation matrix."
    )
    parser.add_argument("--protein", required=True,
                        help="Protein-ligand pair ID (e.g. katA_heme)")
    parser.add_argument("--states", nargs="+", default=ALL_STATES,
                        choices=ALL_STATES,
                        help=f"States to generate (default: all — {ALL_STATES})")
    parser.add_argument("--methods", nargs="+", default=ALL_METHODS,
                        choices=ALL_METHODS,
                        help=f"Methods to generate (default: all — {ALL_METHODS})")
    parser.add_argument("--eps", nargs="+", type=int, default=ALL_EPSILON,
                        help=f"Epsilon values (default: all — {ALL_EPSILON})")
    args = parser.parse_args()

    eps_list = validate_epsilon(args.eps)
    config   = load_config(args.protein)

    total = len(args.states) * len(args.methods) * len(eps_list)
    print(f"\nGenerating inputs for [{args.protein}]")
    print(f"  States:   {args.states}")
    print(f"  Methods:  {args.methods}")
    print(f"  Epsilon:  {eps_list}")
    print(f"  Total inputs to create: {total}\n")

    created = 0
    skipped = 0

    for state in args.states:
        charge, mult = get_charge_mult(config, state)

        # Optimised geometry path
        xyz_path = (RESULTS_DIR / args.protein / state /
                    "geometry_opt" / f"{args.protein}_{state}_opt.xyz")

        if not xyz_path.exists():
            print(f"  WARNING: Optimised geometry not found for state '{state}':")
            print(f"           {xyz_path}")
            print(f"           Run generate_states.py first. Skipping {state}.\n")
            skipped += len(args.methods) * len(eps_list)
            continue

        for method in args.methods:
            template = METHOD_TEMPLATES[method]

            for eps in eps_list:
                out_dir = (RESULTS_DIR / args.protein / state /
                           method / f"eps_{eps}")
                out_dir.mkdir(parents=True, exist_ok=True)

                inp_name = f"{args.protein}_{state}_{method}_eps{eps}.inp"
                inp_path = out_dir / inp_name

                content = template.format(
                    protein_id=args.protein,
                    state=state,
                    eps=eps,
                    charge=charge,
                    mult=mult,
                    xyz_path=str(xyz_path.resolve())
                )
                inp_path.write_text(content)
                created += 1

        print(f"  [{state}] charge={charge} mult={mult} — inputs created")

    print(f"\n{'='*50}")
    print(f"Input generation complete.")
    print(f"  Created: {created}")
    print(f"  Skipped: {skipped} (missing geometry_opt .xyz)")
    print(f"\nFolder structure:")
    print(f"  results/{args.protein}/{{state}}/{{method}}/eps_{{value}}/")
    print(f"{'='*50}\n")

    if skipped > 0:
        print("Run generate_states.py for missing states before running simulations.\n")


if __name__ == "__main__":
    main()
