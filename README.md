# orca_simulation

A two-phase computational pipeline for studying UV-induced electronic excitation in redox-active proteins of *Pseudomonas aeruginosa* PAO1 — a clinically critical multidrug-resistant pathogen.

**Phase 0** mines 6,000+ research papers to identify and rank protein-ligand simulation targets.  
**Phase 1** extracts quantum mechanical cluster models from PDB crystal structures and runs DFT/TDDFT simulations using ORCA.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Phase 0 — Literature Mining](#phase-0--literature-mining)
- [Phase 1 — Simulation Pipeline](#phase-1--simulation-pipeline)
- [Configuration](#configuration)
- [Simulation Parameters](#simulation-parameters)
- [Troubleshooting](#troubleshooting)

---

## Project Overview

*Pseudomonas aeruginosa* resists antibiotics through redox-linked mechanisms — ROS detoxification enzymes, efflux pumps driven by proton motive force, and quorum sensing-regulated virulence. This project uses Density Functional Theory (DFT) to compute electronic excitation properties of 43 redox protein active sites across:

- **3 molecular states** — neutral, radical, reduced
- **4 methods** — DFT, Constrained DFT (CDFT), Spin-DFT (SDFT), Time-Dependent DFT (TDDFT)
- **5 dielectric constants** — ε = 10, 20, 30, 40, 80 (CPCM implicit solvent)

Total: **2,580 ORCA single-point calculations** across 43 protein-ligand pairs.

All calculations use **B3LYP/def2-TZVP** as implemented in ORCA 6.x.

---

## Repository Structure

```
orca_simulation/
│
├── literature_mining/              # Phase 0 — literature mining
│   ├── pubmed_fetch.py             # Step 1: fetch papers from PubMed + EuropePMC
│   ├── extract_targets.py          # Step 2: extract genes, proteins, ligands
│   ├── score_targets.py            # Step 3: rank proteins by simulation priority
│   ├── report_targets.py           # Step 4: generate Excel report
│   ├── raw_papers/                 # gitignored — raw JSON from APIs
│   ├── extracted/                  # gitignored — intermediate JSON files
│   └── redox_amr_targets_Pa.xlsx   # output report (committed)
│
├── pipeline/                       # Phase 1 — simulation pipeline
│   ├── setup_proteins.py           # Pre-step: create protein folders + config.json
│   ├── fetch_structures.py         # Step 1: download PDB structures from RCSB
│   ├── extract_qm.py               # Step 2: extract QM cluster from PDB
│   ├── generate_states.py          # Step 3: generate geometry optimisation inputs
│   ├── generate_inputs.py          # Step 4: generate all 60 ORCA inputs per protein
│   └── run_batch.py                # Step 5: run ORCA simulations
│
├── proteins/                       # one folder per protein-ligand pair
│   └── {protein_id}/
│       ├── config.json             # PDB ID, ligand residue, charge/multiplicity
│       └── docked.pdb              # crystal structure (downloaded by fetch_structures.py)
│
├── qm_regions/                     # gitignored — extracted QM cluster .xyz files
├── results/                        # ORCA inputs committed; outputs gitignored
│   └── {protein_id}/
│       └── {state}/
│           ├── geometry_opt/       # geometry optimisation files
│           └── {method}/
│               └── eps_{value}/    # single-point input and output files
│
├── .gitignore
└── README.md
```

---

## Dependencies

### Python packages

```bash
pip install requests biopython openpyxl
```

| Package | Version | Used by |
|---|---|---|
| requests | any | pubmed_fetch.py, fetch_structures.py |
| biopython | ≥1.79 | extract_qm.py (PDBParser, NeighborSearch) |
| openpyxl | ≥3.0 | report_targets.py |

### ORCA

ORCA 5.x or 6.x must be installed and accessible from your terminal:

```bash
orca --version
```

If this fails, add your ORCA installation directory to PATH:

```bash
# Linux / macOS
export PATH="/path/to/orca:$PATH"
# Add to ~/.bashrc or ~/.zshrc to make permanent

# Windows (PowerShell)
$env:PATH += ";C:\path\to\orca"
```

> **Important:** When `nprocs > 1`, ORCA requires a full absolute path to its executable. The pipeline handles this automatically using `shutil.which("orca")` — no manual configuration needed as long as ORCA is in your PATH.

---

## Phase 0 — Literature Mining

Run these four scripts in order, once, from the `literature_mining/` directory.

### Step 1 — Fetch papers

```bash
cd literature_mining
python pubmed_fetch.py
```

Queries PubMed and EuropePMC with 42 search strings covering all aspects of *P. aeruginosa* redox biology and AMR. Deduplicates by PMID then DOI.

**Output:** `raw_papers/pubmed_raw.json`, `raw_papers/epmc_raw.json`  
**Expected result:** ~6,000 unique papers  
**Runtime:** 15–30 minutes (rate-limited to 3 req/s for NCBI)

Optional arguments:
```bash
python pubmed_fetch.py --max-per-query 200   # increase papers per query (default: 150)
```

---

### Step 2 — Extract targets

```bash
python extract_targets.py
```

Identifies gene mentions, assigns simulation ligands from the `GENE_BIOCHEMISTRY` table (biochemistry-derived, not statistics-derived), and separates noise ligands (pyocyanin, iron, phenazine, pyoverdine, pyochelin) from signal co-mentions.

**Output:** `extracted/extracted_targets.json`

---

### Step 3 — Score targets

```bash
python score_targets.py
```

Ranks all detected genes using a weighted formula:

| Component | Weight | Description |
|---|---|---|
| Citation frequency | 35% | Normalised paper count |
| Simulatable ligand | 20% | Has a DFT-suitable cofactor |
| PDB availability | 20% | Crystal structure exists |
| Category relevance | 15% | How directly linked to AMR |
| Ligand diversity | 10% | Signal co-mention breadth |

**Output:** `extracted/scored_targets.json`

---

### Step 4 — Generate report

```bash
python report_targets.py
```

Produces a formatted 4-sheet Excel workbook with ranked proteins, ligand mentions, and recommended simulation pairs with clickable PubMed/PDB hyperlinks.

**Output:** `redox_amr_targets_Pa.xlsx`

---

## Phase 1 — Simulation Pipeline

All scripts are run from the **repository root**, not from inside `pipeline/`.

### Pre-step — Create protein folders

```bash
python pipeline/setup_proteins.py
```

Auto-creates `proteins/{protein_id}/config.json` for all 43 simulatable proteins. The `PROTEIN_REGISTRY` inside this script contains PDB IDs, ligand residue names, and charge/multiplicity defaults for each protein.

```bash
# Preview without writing
python pipeline/setup_proteins.py --dry-run

# Create specific proteins only
python pipeline/setup_proteins.py --only katA_heme lasR_3OC12HSL

# Overwrite existing configs
python pipeline/setup_proteins.py --force
```

After running, the script prints two lists:
- Proteins with **empty `pdb_id`** — fill these in `config.json` before Step 1
- Proteins with **`review_multiplicity: true`** — metal-containing systems where default charge/multiplicity must be verified against experimental data before running

---

### Step 1 — Download structures

```bash
python pipeline/fetch_structures.py --protein katA_heme
```

Downloads the PDB crystal structure from RCSB and saves it to `proteins/{protein_id}/docked.pdb`. Requires `pdb_id` to be set in `config.json`.

```bash
# Download all proteins (only those with pdb_id set)
python pipeline/fetch_structures.py

# Re-download even if file exists
python pipeline/fetch_structures.py --protein katA_heme --force

# Show download status of all proteins without downloading
python pipeline/fetch_structures.py --report
```

---

### Step 2 — Extract QM region

```bash
python pipeline/extract_qm.py --protein katA_heme --chain A --cutoff 3.5
```

Extracts the ligand and all protein residues within `--cutoff` angstroms from a single chain, writing an XYZ coordinate file for the QM cluster.

```bash
# List all chains and HETATM residues in the PDB (run this first)
python pipeline/extract_qm.py --protein katA_heme --list-chains

# Extract with default settings (chain A, 3.5 Å)
python pipeline/extract_qm.py --protein katA_heme
```

**Output:** `qm_regions/{protein_id}.xyz`

> **Always use `--list-chains` first** to confirm your PDB has the expected chain and ligand residue name. Many PDB structures are multimers (dimers, tetramers) — without `--chain`, all subunits are extracted, producing a QM region too large for DFT.

Recommended cutoff values:

| Cutoff | Typical atom count | Suitable for |
|---|---|---|
| 3.5 Å | 100–180 atoms | Desktop workstation (recommended) |
| 4.0 Å | 150–250 atoms | Workstation with ≥32 GB RAM |
| 5.0 Å | 250–400+ atoms | HPC cluster only |

The script aborts if the QM region exceeds 500 atoms and warns above 200 atoms.

---

### Step 3 — Geometry optimisation (optional)

```bash
# Generate input files and review them first
python pipeline/generate_states.py --protein katA_heme --no-run

# Then run geometry optimisation
python pipeline/generate_states.py --protein katA_heme
```

> **Note on geometry optimisation:** For large QM clusters (>100 atoms) with B3LYP/def2-TZVP on desktop hardware, geometry optimisation can take days to weeks. It is standard practice in cluster-model DFT to use crystal structure coordinates directly as input for single-point calculations. To skip geometry optimisation, copy the QM region xyz file to the geometry_opt folder manually:

```bash
# Linux / macOS
mkdir -p results/katA_heme/neutral/geometry_opt
mkdir -p results/katA_heme/radical/geometry_opt
mkdir -p results/katA_heme/reduced/geometry_opt
cp qm_regions/katA_heme.xyz results/katA_heme/neutral/geometry_opt/katA_heme_neutral_opt.xyz
cp qm_regions/katA_heme.xyz results/katA_heme/radical/geometry_opt/katA_heme_radical_opt.xyz
cp qm_regions/katA_heme.xyz results/katA_heme/reduced/geometry_opt/katA_heme_reduced_opt.xyz

# Windows (PowerShell)
New-Item -ItemType Directory -Force results\katA_heme\neutral\geometry_opt
New-Item -ItemType Directory -Force results\katA_heme\radical\geometry_opt
New-Item -ItemType Directory -Force results\katA_heme\reduced\geometry_opt
Copy-Item qm_regions\katA_heme.xyz results\katA_heme\neutral\geometry_opt\katA_heme_neutral_opt.xyz
Copy-Item qm_regions\katA_heme.xyz results\katA_heme\radical\geometry_opt\katA_heme_radical_opt.xyz
Copy-Item qm_regions\katA_heme.xyz results\katA_heme\reduced\geometry_opt\katA_heme_reduced_opt.xyz
```

Optional arguments:
```bash
# Run only specific states
python pipeline/generate_states.py --protein katA_heme --states neutral radical
```

---

### Step 4 — Generate simulation inputs

```bash
python pipeline/generate_inputs.py --protein katA_heme
```

Creates all 60 ORCA input files for every combination of state × method × epsilon. Requires the geometry_opt xyz files from Step 3.

```bash
# Generate only DFT inputs first
python pipeline/generate_inputs.py --protein katA_heme --methods DFT

# Generate only neutral state
python pipeline/generate_inputs.py --protein katA_heme --states neutral

# Generate specific epsilon values only
python pipeline/generate_inputs.py --protein katA_heme --eps 10 20
```

Valid epsilon values: `10 20 30 40 80` — any other value is rejected.

**Output structure:**
```
results/katA_heme/
  neutral/
    DFT/eps_10/katA_heme_neutral_DFT_eps10.inp
    DFT/eps_20/katA_heme_neutral_DFT_eps20.inp
    ...
    TDDFT/eps_80/katA_heme_neutral_TDDFT_eps80.inp
  radical/  ...
  reduced/  ...
```

> **CDFT note:** The CDFT input files contain a placeholder atom index (`atoms_alpha 1 end`) that constrains only atom 1. Before running any CDFT calculation, open the generated `.inp` file and replace `1` with the actual atom index range of your ligand.

---

### Step 5 — Run simulations

```bash
# Always do a dry run first to verify job list
python pipeline/run_batch.py --protein katA_heme --method DFT --dry-run

# Run all DFT jobs
python pipeline/run_batch.py --protein katA_heme --method DFT

# Run everything for one protein
python pipeline/run_batch.py --protein katA_heme

# Filter by state and epsilon range
python pipeline/run_batch.py --protein katA_heme --state neutral --method TDDFT --eps-start 10 --eps-end 40
```

ORCA output files (`.out`, `.gbw`) are written to the same folder as the `.inp` file.

**Verifying results:**
```bash
# Check if a job completed successfully
grep "ORCA TERMINATED NORMALLY" results/katA_heme/neutral/DFT/eps_10/katA_heme_neutral_DFT_eps10.out
grep "FINAL SINGLE POINT ENERGY"  results/katA_heme/neutral/DFT/eps_10/katA_heme_neutral_DFT_eps10.out
```

---

## Configuration

### config.json format

Each protein has a `proteins/{protein_id}/config.json`:

```json
{
  "protein_id": "katA_heme",
  "gene": "katA",
  "protein_name": "Catalase KatA",
  "ligand_name": "heme",
  "pdb_id": "4E37",
  "ligand_residue": "HEM",
  "pdb_notes": "PDB 4E37 — Pa PAO1 KatA with heme b + NADPH confirmed.",
  "review_multiplicity": true,
  "multiplicity_notes": "Fe(III) high-spin d5, S=5/2 — verify spin state.",
  "states": {
    "neutral": { "charge": 1,  "multiplicity": 6 },
    "radical": { "charge": 2,  "multiplicity": 5 },
    "reduced": { "charge": 0,  "multiplicity": 5 }
  }
}
```

### Charge and multiplicity

ORCA enforces a strict parity rule: **even number of electrons requires odd multiplicity** and vice versa. Always verify before running:

1. Count total electrons: sum of atomic numbers across all atoms in the QM cluster minus the charge
2. If electrons are even → multiplicity must be odd (1, 3, 5, 7...)
3. If electrons are odd → multiplicity must be even (2, 4, 6...)

For metal-containing systems (proteins flagged `review_multiplicity: true`), the spin state depends on the metal oxidation state and ligand field. Refer to the `multiplicity_notes` field in each `config.json` for guidance.

### nprocs and maxcore

Edit the template sections inside `pipeline/generate_states.py` and `pipeline/generate_inputs.py`:

```python
%pal
  nprocs 4        # set to number of physical CPU cores
end

%maxcore 2048    # RAM per process in MB
```

Safe formula for `maxcore`:
```
maxcore = (Total RAM in MB × 0.75) / nprocs
```

| RAM | nprocs | maxcore |
|---|---|---|
| 16 GB | 4 | 3072 MB |
| 32 GB | 4 | 6144 MB |
| 64 GB | 8 | 6144 MB |

Do not set `nprocs` higher than the number of **physical** cores — hyperthreading does not benefit ORCA.

---

## Simulation Parameters

| Parameter | Value |
|---|---|
| DFT functional | B3LYP |
| Basis set | def2-TZVP |
| Implicit solvent | CPCM |
| Epsilon values | 10, 20, 30, 40, 80 |
| Molecular states | neutral, radical, reduced |
| Methods | DFT, CDFT, SDFT, TDDFT |
| TDDFT excited states | 10 roots |
| Geometry | Crystal structure coordinates |
| QM region cutoff | 3.5 Å from ligand (recommended) |

---

## Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| `orca: command not found` | ORCA not in PATH | Add ORCA directory to PATH |
| `ORCA has to be called with full pathname` | PATH resolution fails in MPI | Ensure `shutil.which("orca")` returns a full path — check PATH |
| `multiplicity (N) is even and number of electrons (M) is even → impossible` | Parity violation | Recount electrons and adjust charge or multiplicity to satisfy parity rule |
| `config.json not found` | setup_proteins.py not run | Run `python pipeline/setup_proteins.py` |
| `Ligand residue 'HEM' not found in chain A` | Wrong ligand_residue in config.json | Run `--list-chains` to see actual HETATM residue names |
| QM region has 1000+ atoms | Multimer PDB, all chains extracted | Always specify `--chain A` with `extract_qm.py` |
| Geometry optimisation running for hours | Normal for large clusters — use single-point instead | Copy xyz from qm_regions/ to geometry_opt/ directly (see Step 3) |
| `pdb_id is empty` | PDB not set in config.json | Search rcsb.org, fill `pdb_id` and `ligand_residue` in config.json |

---

## .gitignore

The following are excluded from version control:

```
# ORCA output files (large)
*.out
*.gbw
*.densities
*.ges
*.prop
*.engrad
*_trj.xyz
*.hess
*.cis

# Literature mining raw data
literature_mining/raw_papers/
literature_mining/extracted/
```

`.inp` files, `.xyz` geometry files, `config.json`, and `docked.pdb` are all committed.
