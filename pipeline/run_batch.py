"""
run_batch.py  —  Phase 1, Step 4
Runs ORCA on all .inp files matching the specified filters.
Team members each run different slices of the simulation matrix.

Usage examples:
    # Run all DFT neutral simulations for katA_heme
    python run_batch.py --protein katA_heme --state neutral --method DFT

    # Run only epsilon 10 and 20 for radical TDDFT
    python run_batch.py --protein katA_heme --state radical --method TDDFT --eps-start 10 --eps-end 20

    # Run all methods for one state (useful if one team member owns a full state)
    python run_batch.py --protein katA_heme --state reduced

    # Dry run — list jobs without running
    python run_batch.py --protein katA_heme --state neutral --method DFT --dry-run
"""

import argparse
import subprocess
import sys
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "results"

ALL_STATES  = ["neutral", "radical", "reduced"]
ALL_METHODS = ["DFT", "CDFT", "SDFT", "TDDFT"]
ALL_EPSILON = [10, 20, 30, 40, 80]


def get_eps_subset(eps_start: int, eps_end: int) -> list:
    subset = [e for e in ALL_EPSILON if eps_start <= e <= eps_end]
    if not subset:
        print(f"ERROR: No valid epsilon values between {eps_start} and {eps_end}.")
        print(f"       Valid values are: {ALL_EPSILON}")
        sys.exit(1)
    return subset


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


def run_orca(inp_path: Path) -> str:
    """Run ORCA on a single input file. Returns 'ok' or 'failed'."""
    out_path = inp_path.with_suffix(".out")
    orca_exe = get_orca_path()
    try:
        result = subprocess.run(
            [orca_exe, str(inp_path.resolve())],
            stdout=open(out_path, "w"),
            stderr=subprocess.STDOUT,
            cwd=inp_path.parent
        )
        return "ok" if result.returncode == 0 else "failed"
    except Exception as e:
        print(f"\nERROR running ORCA on {inp_path.name}: {e}")
        return "failed"


def collect_jobs(protein: str, states: list, methods: list, eps_list: list) -> list:
    """Collect all .inp paths matching the filter criteria."""
    jobs = []
    for state in states:
        for method in methods:
            for eps in eps_list:
                inp_dir = RESULTS_DIR / protein / state / method / f"eps_{eps}"
                inp_name = f"{protein}_{state}_{method}_eps{eps}.inp"
                inp_path = inp_dir / inp_name
                jobs.append((state, method, eps, inp_path))
    return jobs


def main():
    parser = argparse.ArgumentParser(
        description="Run ORCA batch simulations for a protein-ligand pair.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Team assignment example (5 members, katA_heme):
  Member 1: --state neutral --method DFT
  Member 2: --state neutral --method TDDFT
  Member 3: --state radical --method DFT --method CDFT
  Member 4: --state reduced
  Member 5: --state neutral --method CDFT --method SDFT
        """
    )
    parser.add_argument("--protein", required=True,
                        help="Protein-ligand pair ID (e.g. katA_heme)")
    parser.add_argument("--state", nargs="+", default=ALL_STATES,
                        choices=ALL_STATES,
                        help="Molecular state(s) to run (default: all)")
    parser.add_argument("--method", nargs="+", default=ALL_METHODS,
                        choices=ALL_METHODS,
                        help="Method(s) to run (default: all)")
    parser.add_argument("--eps-start", type=int, default=10,
                        help="Start of epsilon range (default: 10)")
    parser.add_argument("--eps-end", type=int, default=80,
                        help="End of epsilon range inclusive (default: 80)")
    parser.add_argument("--dry-run", action="store_true",
                        help="List jobs without running ORCA")
    args = parser.parse_args()

    eps_list = get_eps_subset(args.eps_start, args.eps_end)
    jobs     = collect_jobs(args.protein, args.state, args.method, eps_list)

    # Check all input files exist before starting
    missing = [(s, m, e, p) for s, m, e, p in jobs if not p.exists()]
    if missing:
        print(f"\nERROR: {len(missing)} input file(s) not found.")
        print("       Run generate_inputs.py first.\n")
        for s, m, e, p in missing[:10]:
            print(f"  Missing: {p}")
        if len(missing) > 10:
            print(f"  ... and {len(missing)-10} more.")
        sys.exit(1)

    print(f"\nORCA Batch Runner — [{args.protein}]")
    print(f"  States:   {args.state}")
    print(f"  Methods:  {args.method}")
    print(f"  Epsilon:  {eps_list}")
    print(f"  Jobs:     {len(jobs)}")
    if args.dry_run:
        print(f"  Mode:     DRY RUN (no ORCA execution)\n")
    else:
        print()

    ok_list     = []
    failed_list = []

    for i, (state, method, eps, inp_path) in enumerate(jobs, 1):
        label = f"[{i}/{len(jobs)}] {state}/{method}/eps_{eps}"

        if args.dry_run:
            print(f"  {label}  →  {inp_path}")
            continue

        print(f"  {label} ...", end="", flush=True)
        status = run_orca(inp_path)

        if status == "ok":
            print(" done")
            ok_list.append(inp_path)
        else:
            print(" FAILED")
            failed_list.append(inp_path)

    if args.dry_run:
        print(f"\nDry run complete. {len(jobs)} jobs listed.")
        return

    print(f"\n{'='*55}")
    print(f"Batch complete.")
    print(f"  Successful: {len(ok_list)}/{len(jobs)}")
    print(f"  Failed:     {len(failed_list)}/{len(jobs)}")
    if failed_list:
        print("\nFailed jobs:")
        for p in failed_list:
            print(f"  {p}")
        print("\nCheck the corresponding .out files for ORCA error messages.")
    print(f"{'='*55}\n")

    if failed_list:
        sys.exit(1)


if __name__ == "__main__":
    main()
