# OpenMM System Preparation Pipeline

A lightweight Python program to prepare solvated molecular systems for simulation using **OpenMM**, with optional ligand parametrization via **Espaloma, SMIRNOFF, or GAFF**.

This tool converts a receptor (PDB) and optional ligand (SDF) into fully parameterized systems ready for **OpenMM**, **AMBER**, and **GROMACS** workflows.  

It can be run from the command-line, given that it was properly installed, e.g., via pip. Besides that, we support usage in python scripts and notebooks (this feature is under development).

---

## Features

-  Protein preparation using `pdbfixer`
  - Missing atoms & hydrogens added
  - Non-standard residues replaced
-  Optional ligand support (SDF input)
-  Multiple ligand force fields:
   - Espaloma
   - SMIRNOFF (OpenFF)
   - GAFF (default)
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
Either you already have an environment with all the openmm tools, or you create a new environment just for this tool

### Step 1 -- install dependencies
```bash
micromamba create -n simprepper python=3.10
micromamba activate simprepper

# mandatory
micromamba install -c conda-forge openmm   openmmtools   pdbfixer   openff-toolkit   openmmforcefields   parmed   rdkit     
# optional
micromamba install -c conda-forge espaloma
```
or alternatively:

```bash
pip install espaloma
```

### Step 2 -- install `simprepper`
After installing the dependencies, actually install `simprepper`:
```bash
pip install .
# or for developers:
pip install -e .
```
Afterwards, to check your installation, try:
```bash
simprepper --help
```

---

## Usage

Checkout the subdirectory [examples](examples/), for comprehensive examples

### Protein only (apo system)

```bash
simprepper -r receptor.pdb
```

### Protein + ligand

```bash
simprepper -r receptor.pdb -l ligand.sdf
```

### Optional logging level

```bash
simprepper -r receptor.pdb -l ligand.sdf -L debug
```

### Specifying simulation settings

Simulation settings can be provided through a configuration file:

```bash
simprepper -r receptor.pdb -i config.ini
```

When a config file is supplied, the simulation parameters are read from this file rather than from the command line.

To generate an example configuration file containing the default simulation settings, run:

```bash
simprepper-example-config -o default_config.ini
```

This will create a template that can be modified and reused for future simulations.


The configuration file is organized into sections:
```bash
[system]
membrane_protein = False

[simulation]
nb_cutoff = 1.0  # nanometer
hydrogen_mass = 4.0  # g/mol
timestep = 0.004  # picosecond
temperature = 300.0  # kelvin
ionic_strength = 0.15  # molar
ph = 7.4

[simulation box]
box_shape = cube
padding = 3.0  # nanometer

[forcefields]
ligand_ff = GAFF
protein_ff = amber14/protein.ff14SB.xml
water_ff = amber14/tip3pfb.xml
ion_ff = amber/tip3p_HFE_multivalent.xml
lipid_ff = None
```

Numerical values are stored in a machine-readable format. Unit annotations are included as comments for readability and are ignored during parsing.


### Membrane protein systems

**Important**: The input PDB file for a membrane protein system must be pre-oriented with respect to the membrane, for example using the [PPM server](https://opm.phar.umich.edu/ppm_server) or by obtaining the structure from the [OPM database](https://opm.phar.umich.edu/). In addition, any `DUM` residues (created by PPM) must be removed from the structure before running `simprepper`.

An example of a membrane protein system configuration can be generated using a `--membrane` flag:
```bash
simprepper-example-config -o default_config.ini --membrane
```

The `membrane_protein` flag should be set to `True` in the config file. The setup for membrane systems includes an additional section called `[membrane]` and does not contain the `[simulation box]` section since OpenMM's `addMembrane()` function does not take box parameters.

```bash
[system]
membrane_protein = True

[membrane]
lipid_type = POPC
membrane_center_z = 0.0
minimum_padding = 1.0  # nanometer

[simulation]
nb_cutoff = 1.0  # nanometer
hydrogen_mass = 4.0  # g/mol
timestep = 0.004  # picosecond
temperature = 300.0  # kelvin
ionic_strength = 0.15  # molar
ph = 7.4

[forcefields]
ligand_ff = GAFF
protein_ff = amber14/protein.ff14SB.xml
water_ff = amber14/tip3pfb.xml
ion_ff = amber/tip3p_HFE_multivalent.xml
lipid_ff = amber14/lipid17.xml
```


## Testing

For developers, and to check, if the installation worked out, check the subdirectory [examples](examples/), 
which currently contains three different use-cases.


---

##  Output

All outputs are written to a directory named after the ligand (or receptor if no ligand is provided):

```bash
<system_name>/
│
├── <receptor>_fixed.pdb
├── <system_name>_solvated.pdb
├── gmx (optional)
│   ├── <system_name>_solvated.gro
│   └── <system_name>_solvated.top
├── amber (optional)
│   ├── <system_name>_solvated.rst7
│   └── <system_name>_solvated.prmtop
└── openmm (optional)
    ├── <system_name>_solvated.chk
    └── <system_name>_system.xml
```

Logs are written to:

```
simprepper_logs/<system_name>.log
```

---

##  Default Simulation Settings

| Parameter             | Value                 |
|-----------------------|-----------------------|
| Protein Force field   | AMBER ff14SB          |
| Water model           | TIP3P-FB              |
| Ion Parameters        | HFE multivalent       |
| Lipid Force Field     | None                  |
| Box shape             | Cube                  |
| Padding               | 3 nm                  |
| Ionic strength        | 0.15 M                |
| pH                    | 7.4                   |
| Cutoff                | 1.0 nm                |
| Constraints           | HBonds                |
| Hydrogen mass         | 4 amu (HMR)           |
| Timestep              | 4 fs (0.004 ps)       |
| Temperature           | 300 K                 |

---

## Force Fields

###  Available force fields

The list of available force fields can be printed by running:
```bash
simprepper-forcefields
```

### Ligand Force Fields

Controlled internally via:

```python
--ligand_ff = 'espaloma'
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
5. Solvation (TIP3P-FB + ions or lipids + TIP3P-FB + ions)
6. System creation (OpenMM)
7. Export (OpenMM / AMBER / GROMACS)

---

## 🔧 Troubleshooting

### Ligand not loading

- Ensure SDF is valid and contains coordinates as well as hydrogens!

### CUDA not detected

- Check:
  ```bash
  nvidia-smi -L
  ```
This should list the NVIDIA GPUs in your machine.

### Espaloma issues

Ensure that you have installed espaloma, if you want to use it as a ligand force field.
```bash
pip install espaloma
```

---

## ‍ Author

Joe Loeffler, Monica Fernandez-Quintero, Patrick K. Quoika, Julia Bandera

---

## License

MIT License

