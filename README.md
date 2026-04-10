# ORCA Simulation Pipeline

**Project:** UV-Induced Excitation in Bacterial Proteins  
**Target Organism:** *Pseudomonas aeruginosa*  

## Overview
This repository contains a comprehensive computational pipeline for studying redox biology and antimicrobial resistance in *Pseudomonas aeruginosa*. The pipeline is divided into two main phases:
1. **Literature Mining (Phase 0):** Automated extraction and scoring of protein-ligand targets from scientific literature (PubMed and EuropePMC).
2. **Simulation Pipeline (Phase 1):** Automated generation, execution, and management of quantum mechanical (QM) simulations using ORCA. 

The simulation phase evaluates proteins across multiple molecular states (neutral, radical, reduced) and methodologies (DFT, CDFT, SDFT, TDDFT) under various dielectric constants (epsilon values).

---

## Directory Structure

```text
orca_simulation/
  ├── literature_mining/       # Phase 0: Target discovery and scoring
  ├── pipeline/                # Phase 1: Automation scripts for ORCA setup and execution
  ├── proteins/                # Target-specific configuration and structural data (.pdb)
  ├── qm_regions/              # Extracted QM regions (.xyz) [Git ignored]
  └── results/                 # ORCA output files organized by state/method/eps [Git ignored]
```

---

## Prerequisites and Dependencies

The pipeline relies on a mix of standard Python libraries and a few external packages.

**External Packages:**
- `requests` (for fetching literature and PDB structures)
- `openpyxl` (for generating Excel reports)
- `biopython` (for extracting QM regions from PDB files)

Install them via pip:
```bash
pip install requests openpyxl biopython
```

**Standard Library Modules Used:**
- `json`, `re`, `pathlib`, `argparse`, `subprocess`, `collections`

**Computational Software:**
- **ORCA** (Quantum Chemistry software) must be installed and accessible in your system path.

---

## Pipeline Execution Order

### Phase 0: Literature Mining
Run these scripts to discover and rank potential protein targets.

1. Fetch literature data:
   ```bash
   python literature_mining/pubmed_fetch.py
   ```
2. Extract targets:
   ```bash
   python literature_mining/extract_targets.py
   ```
3. Score the extracted targets:
   ```bash
   python literature_mining/score_targets.py
   ```
4. Generate the final report:
   ```bash
   python literature_mining/report_targets.py
   ```

### Phase 1: Simulation (Per Protein-Ligand Pair)
For a given protein (e.g., `katA_heme`), execute the following steps in order:

1. **Setup & Fetch:** Ensure the configuration exists and fetch the structure.
   ```bash
   python pipeline/setup_proteins.py
   python pipeline/fetch_structures.py --protein <protein_id>
   ```
   *(Alternatively, place a custom docked `.pdb` file at `proteins/<protein_id>/docked.pdb`)*

2. **Extract QM Region:**
   ```bash
   python pipeline/extract_qm.py --protein <protein_id>
   ```

3. **Generate Optimized Geometries (States):**
   ```bash
   python pipeline/generate_states.py --protein <protein_id>
   ```

4. **Generate ORCA Inputs:**
   ```bash
   python pipeline/generate_inputs.py --protein <protein_id>
   ```

5. **Run Batch Simulations:**
   Execute the simulations for a specific state and method.
   ```bash
   python pipeline/run_batch.py --protein <protein_id> --state neutral --method DFT
   ```

---

## Global Simulation Parameters

The pipeline enforces consistent parameters across all generated ORCA inputs:

- **Hardware:** 4 processors (`nprocs`), 2048 MB memory per core (`maxcore`)
- **Method/Basis:** B3LYP functional with def2-TZVP basis set
- **Simulated States:** `neutral`, `radical`, `reduced`
- **Simulation Methods:** `DFT`, `CDFT`, `SDFT`, `TDDFT`
- **Dielectric Constants (Epsilon):** `10`, `20`, `30`, `40`, `80`

Each specific protein requires a `config.json` inside its respective `proteins/<protein_id>/` folder to define the required **charge** and **multiplicity** for each state.