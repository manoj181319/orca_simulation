"""
generate_states.py

Script 2 — State Generation
-----------------------------
Reads:   01_qm_region/qm_region.xyz
Creates: ORCA geometry optimization input files for 3 states
Runs:    ORCA optimization for each state
Saves:   neutral_opt.xyz, radical_opt.xyz, reduced_opt.xyz
         into 02_state_generation/<state>/

States:
    Neutral  → charge  0, multiplicity 1
    Radical  → charge +1, multiplicity 2
    Reduced  → charge -1, multiplicity 1
"""

import os
import subprocess
import shutil

# ─────────────────────────────────────────────
# CONFIGURATION — edit these before running
# ─────────────────────────────────────────────

QM_XYZ      = "01_qm_region/qm_region.xyz"       # input from Script 1
STATE_DIR   = "02_state_generation"               # output base folder
ORCA_CMD    = "orca"                              # ORCA is in system PATH

FUNCTIONAL  = "B3LYP"                            # DFT functional
BASIS_SET   = "def2-TZVP"                        # basis set
NPROCS      = 4                                  # number of CPU cores
MAXCORE     = 2048                               # memory per core in MB

# State definitions — adjust charge/multiplicity if guide specifies otherwise
STATES = [
    {"name": "neutral", "charge":  0, "multiplicity": 1},
    {"name": "radical", "charge": +1, "multiplicity": 2},
    {"name": "reduced", "charge": -1, "multiplicity": 1},
]

# ─────────────────────────────────────────────
# STEP 1 — Read QM region XYZ file
# ─────────────────────────────────────────────

def read_xyz(xyz_file):
    """
    Read XYZ file and return:
        - n_atoms: number of atoms
        - comment: comment line
        - coords: list of (element, x, y, z) tuples
    """
    print(f"[1] Reading QM region: {xyz_file}")
    if not os.path.exists(xyz_file):
        raise FileNotFoundError(
            f"QM region file not found: {xyz_file}\n"
            f"Run extract_qm.py first."
        )

    with open(xyz_file, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    coords  = []

    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        if len(parts) == 4:
            element = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            coords.append((element, x, y, z))

    print(f"     Atoms read: {n_atoms}")
    print(f"     Comment   : {comment}")
    return n_atoms, comment, coords

# ─────────────────────────────────────────────
# STEP 2 — Generate ORCA input file
# ─────────────────────────────────────────────

def generate_orca_input(state, coords, input_path):
    """
    Write ORCA geometry optimization input file for a given state.
    
    Input file structure:
        - keywords line (functional, basis, OPT, CPCM)
        - parallel cores block
        - memory block
        - coordinate block (charge, multiplicity, xyz)
    """
    name         = state["name"]
    charge       = state["charge"]
    multiplicity = state["multiplicity"]

    print(f"     Generating input: {input_path}")

    with open(input_path, "w") as f:
        # Keywords line
        f.write(f"! {FUNCTIONAL} {BASIS_SET} OPT CPCM\n")
        f.write(f"\n")

        # Parallel block
        f.write(f"%pal\n")
        f.write(f"  nprocs {NPROCS}\n")
        f.write(f"end\n")
        f.write(f"\n")

        # Memory block
        f.write(f"%maxcore {MAXCORE}\n")
        f.write(f"\n")

        # Coordinate block
        f.write(f"* xyz {charge} {multiplicity}\n")
        for (element, x, y, z) in coords:
            f.write(f"  {element:<4s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
        f.write(f"*\n")

    print(f"     Input written : {input_path}")

# ─────────────────────────────────────────────
# STEP 3 — Run ORCA optimization
# ─────────────────────────────────────────────

def run_orca(input_path, output_dir, state_name):
    """
    Run ORCA on the given input file.
    Working directory is set to output_dir so all ORCA files
    land in the correct output folder.
    """
    print(f"     Running ORCA for: {state_name}")
    print(f"     Output folder  : {output_dir}")

    input_abs  = os.path.abspath(input_path)
    output_abs = os.path.abspath(output_dir)

    # ORCA must be run with full path to input file
    # stdout is redirected to .out file inside output folder
    out_file = os.path.join(output_abs, f"{state_name}.out")

    with open(out_file, "w") as f_out:
        result = subprocess.run(
            [ORCA_CMD, input_abs],
            stdout=f_out,
            stderr=subprocess.STDOUT,
            cwd=output_abs   # ORCA writes all files here
        )

    if result.returncode != 0:
        print(f"     ⚠ WARNING: ORCA returned non-zero exit code for {state_name}")
        print(f"       Check output: {out_file}")
    else:
        print(f"     ✓ ORCA completed: {state_name}")

    return out_file

# ─────────────────────────────────────────────
# STEP 4 — Extract optimized XYZ from ORCA output
# ─────────────────────────────────────────────

def extract_optimized_xyz(output_dir, state_name, save_path):
    """
    ORCA writes the optimized geometry to a file named:
        <inputname>.xyz  inside the output directory

    This function locates that file and copies it to the
    state generation folder as <state>_opt.xyz
    """
    print(f"     Extracting optimized geometry for: {state_name}")

    # ORCA names the xyz file after the input file basename
    expected_xyz = os.path.join(output_dir, f"{state_name}.xyz")

    if not os.path.exists(expected_xyz):
        print(f"     ⚠ WARNING: Optimized XYZ not found: {expected_xyz}")
        print(f"       ORCA may have failed. Check the .out file.")
        return False

    shutil.copy(expected_xyz, save_path)
    print(f"     ✓ Optimized geometry saved: {save_path}")
    return True

# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    print("=" * 55)
    print("  State Generation — generate_states.py")
    print("=" * 55)

    # Step 1: Read QM region
    n_atoms, comment, coords = read_xyz(QM_XYZ)

    # Process each state
    for state in STATES:
        name = state["name"]

        print(f"\n{'─' * 55}")
        print(f"  Processing state: {name.upper()}")
        print(f"  Charge: {state['charge']}  |  Multiplicity: {state['multiplicity']}")
        print(f"{'─' * 55}")

        # Define paths for this state
        input_dir  = os.path.join(STATE_DIR, name, "input")
        output_dir = os.path.join(STATE_DIR, name, "output")
        input_path = os.path.join(input_dir, f"{name}.inp")
        opt_xyz    = os.path.join(STATE_DIR, name, f"{name}_opt.xyz")

        # Ensure folders exist
        os.makedirs(input_dir,  exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        # Step 2: Generate ORCA input file
        generate_orca_input(state, coords, input_path)

        # Step 3: Run ORCA
        run_orca(input_path, output_dir, name)

        # Step 4: Extract optimized geometry
        extract_optimized_xyz(output_dir, name, opt_xyz)

    # Summary
    print(f"\n{'=' * 55}")
    print("  State generation complete.")
    print(f"  Optimized geometries saved in: {STATE_DIR}/")
    print(f"{'=' * 55}")

    # Check which states completed successfully
    print("\n  Summary:")
    for state in STATES:
        name     = state["name"]
        opt_path = os.path.join(STATE_DIR, name, f"{name}_opt.xyz")
        status   = "✓ OK" if os.path.exists(opt_path) else "✗ MISSING"
        print(f"    {name:<10s} → {status}")

if __name__ == "__main__":
    main()
