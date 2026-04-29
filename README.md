# OpenMM System Preparation Pipeline

A lightweight Python script to prepare solvated molecular systems for simulation using **OpenMM**, with optional ligand parametrization via **Espaloma, SMIRNOFF, or GAFF**.

This tool converts a receptor (PDB) and optional ligand (SDF) into fully parameterized systems ready for **OpenMM**, **AMBER**, and **GROMACS** workflows.

---

## Features

-  Protein preparation using `pdbfixer`
  - Missing atoms & hydrogens added
  - Non-standard residues replaced
-  Optional ligand support (SDF input)
-  Multiple ligand force fields:
  - Espaloma (default)
  - SMIRNOFF (OpenFF)
  - GAFF
-  Automatic solvation and ion placement
-  Automatic platform selection (CUDA/OpenCL/CPU)
-  Multi-format export:
  - OpenMM (`system.xml`, checkpoint)
  - AMBER (`.prmtop`, `.rst7`)
  - GROMACS (`.gro`, `.top`)
  - PDB (solvated system)

---

## Installation

Recommended: use a micromamba/conda environment.

```bash
micromamba create -n openmm_env python=3.10
micromamba activate openmm_env

micromamba install -c conda-forge     openmm     openmmtools     pdbfixer     parmed     rdkit     openff-toolkit     openmmforcefields
```

For Espaloma:

```bash
pip install espaloma
```

---

## Usage

### Protein only (apo system)

```bash
python openMM_prepare.py -r receptor.pdb
```

### Protein + ligand

```bash
python openMM_prepare.py -r receptor.pdb -l ligand.sdf
```

### Optional logging level

```bash
python openMM_prepare.py -r receptor.pdb -l ligand.sdf -L debug
```

---

##  Output

All outputs are written to a directory named after the ligand (or receptor if no ligand is provided):

```
<system_name>/
│
├── <system_name>_solvated.pdb
├── <system_name>_solvated.prmtop
├── <system_name>_solvated.rst7
├── <system_name>_solvated.gro
├── <system_name>_solvated.top
├── system.xml
├── <system_name>_solvated.chk
└── <receptor>_fixed.pdb
```

Logs are written to:

```
logs/<system_name>.log
```

---

##  Default Simulation Settings

| Parameter          | Value              |
|------------------|-------------------|
| Force field       | AMBER ff14SB      |
| Water model       | TIP3P-FB          |
| Box shape         | Cube              |
| Padding           | 2 nm              |
| Ionic strength    | 0.15 M            |
| Cutoff            | 1.0 nm            |
| Constraints       | HBonds            |
| Hydrogen mass     | 4 amu (HMR)       |
| Timestep          | 4 fs              |

---

## Ligand Force Fields

Controlled internally via:

```python
ligand_ff = 'espaloma'
```

Options:

- `espaloma`
- `SMIRNOFF`
- `GAFF`

---

##  Notes

- Ligand input must be **SDF format** with valid 3D coordinates.
- Undefined stereochemistry is allowed by default.
- If no ligand is provided, the script runs a **protein-only (apo) setup**.
- The force field includes:

```
amber/tip3p_HFE_multivalent.xml
```

This improves ion behavior but may influence systems with unusual ligand charge distributions.

---

## Internals

Pipeline overview:

1. Protein preparation (`pdbfixer`)
2. Ligand loading (RDKit → OpenFF)
3. Ligand parametrization (OpenMMForceFields)
4. System assembly (Modeller)
5. Solvation (TIP3P-FB + ions)
6. System creation (OpenMM)
7. Export (OpenMM / AMBER / GROMACS)

---

## 🔧 Troubleshooting

### Ligand not loading

- Ensure SDF is valid and contains coordinates

### CUDA not detected

- Check:
  ```bash
  nvidia-smi
  ```

### Espaloma issues

```bash
pip install espaloma
```

---

## ‍ Author

Joe Loeffler, Monica Fernandez-Quintero 

---

## License

MIT License

