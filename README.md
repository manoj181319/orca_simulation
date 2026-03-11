# ⚗️ ORCA QM Simulation Pipeline

> Automation pipeline for QM region extraction, input file generation, and batch execution of ORCA quantum chemistry simulations for UV-induced excitation studies in bacterial proteins.

---

## 📌 Project Overview

This project automates the complete workflow for running quantum chemistry simulations using ORCA. The study focuses on **UV-induced excitation in bacterial proteins**, analyzing how electron behavior changes as molecules transition from ground state to excited state.

The pipeline handles:
- Extracting the QM region from a docked protein complex
- Generating optimized geometries for three molecular states
- Creating 450+ ORCA input files automatically
- Running batch simulations across multiple methods and solvent environments

---

## 🧪 Simulation Parameters

| Parameter | Value |
|---|---|
| Molecular States | Neutral, Radical, Reduced |
| Methods (Phase 1) | DFT, CDFT, SDFT |
| Methods (Phase 2) | TDDFT *(remote machine, later)* |
| Epsilon Range | 1 → 50 |
| Phase 1 Simulations | **450** (3 states × 3 methods × 50 ε) |
| Functional | B3LYP |
| Basis Set | def2-TZVP |
| Solvent Model | CPCM |
| ORCA nprocs | 4 |
| ORCA maxcore | 2048 MB |

---

## 🗂️ Folder Structure

```
orca_simulation/
│
├── 00_raw_protein/              # Input protein & docked complex files
│   └── docked_complex.pdb
│
├── 01_qm_region/                # Extracted QM region coordinates
│   └── qm_region.xyz
│
├── 02_state_generation/         # Geometry optimization for each state
│   ├── neutral/
│   │   ├── input/               # ORCA optimization input
│   │   ├── output/              # ORCA optimization output
│   │   └── neutral_opt.xyz      # Final optimized geometry
│   ├── radical/
│   │   ├── input/
│   │   ├── output/
│   │   └── radical_opt.xyz
│   └── reduced/
│       ├── input/
│       ├── output/
│       └── reduced_opt.xyz
│
├── 03_simulations/              # Production simulation runs
│   ├── neutral/
│   │   ├── DFT/
│   │   │   ├── input/           # 50 input files (eps_01.inp ... eps_50.inp)
│   │   │   └── output/          # eps_01/ ... eps_50/ (auto-created at runtime)
│   │   ├── TDDFT/
│   │   ├── CDFT/
│   │   └── SDFT/
│   ├── radical/                 # Same structure as neutral
│   └── reduced/                 # Same structure as neutral
│
├── scripts/                     # All automation scripts
│   ├── extract_qm.py
│   ├── generate_states.py
│   ├── generate_inputs.py
│   └── run_batch.py
│
└── logs/                        # Runtime logs
```

---

## 🔄 Pipeline Flow

```
Docked Complex (PDB)
        │
        ▼
┌─────────────────────┐
│   extract_qm.py     │  Extract QM region using distance cutoff (Biopython)
└─────────────────────┘
        │
        ▼
  qm_region.xyz
        │
        ▼
┌─────────────────────┐
│ generate_states.py  │  Run ORCA geometry optimization for 3 states
└─────────────────────┘
        │
        ▼
  neutral_opt.xyz
  radical_opt.xyz
  reduced_opt.xyz
        │
        ▼
┌─────────────────────┐
│ generate_inputs.py  │  Generate 450 ORCA input files
└─────────────────────┘  (3 states × 3 methods × 50 epsilon)
        │
        ▼
  03_simulations/**/input/
        │
        ▼
┌─────────────────────┐
│   run_batch.py      │  Execute ORCA for all input files
└─────────────────────┘
        │
        ▼
  03_simulations/**/output/
```

---

## 📜 Script Details

### `extract_qm.py`
Reads the docked complex PDB file and extracts the QM region using a **distance cutoff method** around the ligand.

- Library: `Biopython` (`Bio.PDB`, `NeighborSearch`)
- Default cutoff: `5.0 Å` (configurable)
- Output: `01_qm_region/qm_region.xyz`

```bash
python scripts/extract_qm.py
```

---

### `generate_states.py`
Creates three optimized geometries by running ORCA geometry optimization with different charge and multiplicity settings.

| State | Charge | Multiplicity |
|---|---|---|
| Neutral | 0 | 1 |
| Radical | +1 | 2 |
| Reduced | -1 | 1 |

```bash
python scripts/generate_states.py
```

---

### `generate_inputs.py`
Generates all ORCA input files by looping over states, methods, and epsilon values. All 450 files are saved into the correct folders automatically.

```bash
python scripts/generate_inputs.py
```

---

### `run_batch.py`
Runs ORCA on all generated input files. Supports **CLI arguments** for parallel execution across team members.

```bash
# Run everything
python scripts/run_batch.py

# Run specific state
python scripts/run_batch.py --state neutral

# Run specific state and method
python scripts/run_batch.py --state radical --method DFT

# Run specific epsilon range
python scripts/run_batch.py --state reduced --method CDFT --eps-start 1 --eps-end 25
```

---

## ⚡ ORCA Input File Template

All input files follow this structure — only method keywords, charge, multiplicity, and epsilon change:

```
! B3LYP def2-TZVP OPT CPCM

%pal
  nprocs 4
end

%maxcore 2048

%cpcm
  epsilon 20
end

* xyz 0 1
  C    x.xxxxxx   y.yyyyyy   z.zzzzzz
  N    x.xxxxxx   y.yyyyyy   z.zzzzzz
  ...
*
```

| Method | Key Difference |
|---|---|
| DFT | Standard keywords only |
| CDFT | Additional `%cdft` constraints block |
| SDFT | `UKS` keyword + multiplicity change |
| TDDFT | `%tddft` block with `nroots` *(Phase 2)* |

---

## 👥 Parallel Execution Plan

Simulations are distributed by molecular state so team members can run independently:

| Person | State | Simulations |
|---|---|---|
| Person 1 | Neutral | 150 (3 methods × 50 ε) |
| Person 2 | Radical | 150 (3 methods × 50 ε) |
| Person 3 | Reduced | 150 (3 methods × 50 ε) |
| Person 4 & 5 | TDDFT (Phase 2) | 150 on remote machine |

---

## 🖥️ System Requirements

| Component | Spec |
|---|---|
| CPU | Intel i5-10300H (4 physical cores) |
| RAM | 16 GB |
| GPU | Not used (ORCA is CPU-only) |
| OS | Windows / Linux |
| Python | 3.8+ |
| ORCA | Must be installed and in system PATH |
| Biopython | `pip install biopython` |

---

## 📊 Current Status

- [x] Folder structure created
- [x] Repository initialized
- [x] `extract_qm.py` — complete
- [x] `generate_states.py` — complete
- [ ] `generate_inputs.py` — in progress
- [ ] `run_batch.py` — pending
- [ ] Protein structure files from coordinator — pending
- [ ] CDFT fragment definitions — pending confirmation

---

## ⚠️ Important Notes

- **Never commit output files** — `.gbw`, `.out`, `.tmp` files are large and gitignored
- **Always verify ORCA input files** before batch runs — wrong functional = entire month wasted
- **TDDFT is excluded from Phase 1** — will run on remote machine later
- **Do not use GPU settings** — ORCA does not support GPU acceleration

---

*This pipeline is part of an ongoing research project on UV-induced excitation in bacterial proteins.*
