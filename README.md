# ⚗️ ORCA QM Simulation Pipeline

> Automation pipeline for QM region extraction, input file generation, and batch execution of ORCA quantum chemistry simulations.

---

## 📌 Overview

This pipeline automates the complete workflow for running quantum chemistry simulations using ORCA. It handles extraction of the QM region from a docked protein complex, geometry optimization for multiple molecular states, automatic generation of hundreds of input files, and batch execution across multiple computational methods and solvent environments.

---

## 🔄 Pipeline Flow

```
Docked Complex (PDB)
        │
        ▼
┌─────────────────────┐
│   extract_qm.py     │  Extract QM region using distance cutoff
└─────────────────────┘
        │
        ▼
  qm_region.xyz
        │
        ▼
┌─────────────────────┐
│ generate_states.py  │  Run ORCA geometry optimization for each state
└─────────────────────┘
        │
        ▼
  neutral_opt.xyz
  radical_opt.xyz
  reduced_opt.xyz
        │
        ▼
┌─────────────────────┐
│ generate_inputs.py  │  Generate ORCA input files
└─────────────────────┘  (states × methods × epsilon values)
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

## 🗂️ Folder Structure

```
orca_simulation/
│
├── 00_raw_protein/              # Input protein & docked complex files
│   ├── protein.pdb
│   └── docked_complex.pdb
│
├── 01_qm_region/                # Extracted QM region coordinates
│   └── qm_region.xyz
│
├── 02_state_generation/         # Geometry optimization for each state
│   ├── neutral/
│   │   ├── input/
│   │   ├── output/
│   │   └── neutral_opt.xyz
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
│   │   │   ├── input/           # eps_01.inp ... eps_50.inp
│   │   │   └── output/          # eps_01/ ... eps_50/ (auto-created)
│   │   ├── TDDFT/
│   │   ├── CDFT/
│   │   └── SDFT/
│   ├── radical/
│   └── reduced/
│
├── scripts/
│   ├── extract_qm.py
│   ├── generate_states.py
│   ├── generate_inputs.py
│   └── run_batch.py
│
└── logs/
```

---

## 📜 Scripts

### `extract_qm.py`
Reads the docked complex PDB file and extracts the QM region using a distance cutoff around the ligand.

- Library: `Biopython`
- Default cutoff: `5.0 Å` (configurable)
- Output: `01_qm_region/qm_region.xyz`

```bash
python scripts/extract_qm.py
```

Configure before running:
```python
LIGAND_NAME = "LIG"   # 3-letter ligand residue name from PDB
CUTOFF = 5.0          # distance cutoff in angstroms
```

---

### `generate_states.py`
Runs ORCA geometry optimization for three molecular states using the extracted QM region.

| State | Charge | Multiplicity |
|---|---|---|
| Neutral | 0 | 1 |
| Radical | +1 | 2 |
| Reduced | -1 | 1 |

```bash
python scripts/generate_states.py
```

Configure before running:
```python
FUNCTIONAL = "B3LYP"
BASIS_SET  = "def2-TZVP"
NPROCS     = 4
MAXCORE    = 2048
```

---

### `generate_inputs.py`
Generates all ORCA input files by iterating over states, methods, and epsilon values. Files are saved automatically into the correct folders.

```bash
python scripts/generate_inputs.py
```

---

### `run_batch.py`
Runs ORCA on all generated input files. Supports CLI arguments for selective execution.

```bash
# Run all simulations
python scripts/run_batch.py

# Run specific state
python scripts/run_batch.py --state neutral

# Run specific state and method
python scripts/run_batch.py --state radical --method DFT

# Run specific epsilon range
python scripts/run_batch.py --state reduced --method CDFT --eps-start 1 --eps-end 25
```

---

## ⚡ ORCA Input File Structure

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
| TDDFT | `%tddft` block with `nroots` |

---

## 🖥️ Requirements

- Python 3.8+
- Biopython — `pip install biopython`
- ORCA installed and available in system PATH

---

## ⚠️ Notes

- Output files (`.gbw`, `.out`, `.tmp`) are gitignored — they are large and should stay local
- Always verify input files before running batch simulations
- ORCA does not use GPU — CPU cores and RAM are the relevant specs
- Output subfolders (`eps_01/`, `eps_02/`, ...) are created automatically at runtime by `run_batch.py`

*This pipeline is part of an ongoing research project on UV-induced excitation in bacterial proteins.*
