#!/usr/bin/env python
# general imports
import argparse
import os
import sys
import pdbfixer
import parmed
import logging
import numpy as np

# OpenMM imports
from openmm import *
from openmm.app import *
from openmm.unit import *

from openmmtools.utils import get_fastest_platform

# OpenFF-toolkit imports
from openff.toolkit import Molecule
from openff.toolkit import Topology as offTopology
from openff.units.openmm import to_openmm as offquantity_to_openmm

# RDKit imports
from rdkit.Chem import SDMolSupplier

# cython-optimized pairwise distance function
# profiled to run in ~50% of the time of pdist
# much less memory hungry. No storage of all distances.
try:
    # add to PYTHONPATH current workdir and script dir
    sys.path.insert(0, os.curdir)
    sys.path.insert(0, os.path.dirname(__file__))

    from _pwdistance import pw_dist
except ImportError:
    logging.warning("Using numpy/scipy (slower) pwdist routine for simulation setup.")
    from scipy.spatial.distance import pdist

    def pw_dist(xyz_array):
        return np.amax(pdist(xyz_array, "euclidean"))


# %%
parser = argparse.ArgumentParser(
    prog="openMM_prepare",
    description="Script to prepare OpenMM system starting from protein and optional ligand files",
)
parser.add_argument(
    "-l",
    "--lig",
    help="Optional SDF file of the ligand",
    required=False,
    default=None,
)
parser.add_argument(
    "-r",
    "--rec",
    help="PDB file of the receptor",
    required=True,
)
parser.add_argument(
    "-L",
    "--log-level",
    help="Choose the logging level to show",
    choices=["debug", "info", "warning", "error", "critical"],
    default="info",
    required=False,
)

args = parser.parse_args()

# Use ligand name as system name if present, otherwise receptor basename
if args.lig:
    sys_name = os.path.splitext(os.path.basename(args.lig))[0]
else:
    sys_name = os.path.splitext(os.path.basename(args.rec))[0]

rec_basename = os.path.splitext(os.path.basename(args.rec))[0]

os.makedirs(sys_name, exist_ok=True)
os.makedirs("logs", exist_ok=True)

logging.basicConfig(
    level=args.log_level.upper(),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(f"logs/{sys_name}.log", mode="w"),
        logging.StreamHandler(),
    ],
)

nb_cutoff = 1.0  # nanometers
hydrogenMass = 4  # amu
timestep = 0.004  # picoseconds
boxShape = "cube"  # cube, dodecahedron
padding = 2  # nanometers
ionicStrength = 0.15  # molar concentration
ligand_ff = "espaloma"  # espaloma, SMIRNOFF, GAFF


# %%
def prepare_protein(
    pdb_file,
    ignore_missing_residues=True,
    ignore_terminal_missing_residues=True,
    ph=7.4,
):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    # fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer
    fixer.findMissingResidues()

    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = list(fixer.missingResidues.keys())
        for key in keys:
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    return fixer


# %%
def prepare_ligand(lig_sdf, allow_undefined_stereo=True):
    """
    Load ligand from SDF and convert to OpenFF Molecule.
    """
    supplier = SDMolSupplier(lig_sdf, removeHs=False)
    if len(supplier) == 0 or supplier[0] is None:
        raise ValueError(f"Could not read ligand from SDF: {lig_sdf}")

    rdkit_mol = supplier[0]

    ligand = Molecule.from_rdkit(
        rdkit_mol,
        allow_undefined_stereo=allow_undefined_stereo,
    )
    return ligand


# %%
def parametrize_ligand(ligand, forcefield, lig_ff="espaloma"):
    """
    Parameterize ligand and register template generator with the forcefield.
    """
    lig_ff_lower = lig_ff.lower()

    if lig_ff_lower == "espaloma":
        from openmmforcefields.generators import EspalomaTemplateGenerator

        template_generator = EspalomaTemplateGenerator(
            molecules=ligand, forcefield="espaloma-0.3.1"
        )
    elif lig_ff_lower == "smirnoff":
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator

        template_generator = SMIRNOFFTemplateGenerator(
            molecules=ligand, forcefield="openff-1.2.0"
        )
    elif lig_ff_lower == "gaff":
        from openmmforcefields.generators import GAFFTemplateGenerator

        template_generator = GAFFTemplateGenerator(
            molecules=ligand, forcefield="gaff-2.11"
        )
    else:
        raise ValueError(
            "Ligand forcefield must be one of: espaloma, SMIRNOFF, or GAFF"
        )

    forcefield.registerTemplateGenerator(template_generator.generator)

    ligand_off_topology = offTopology.from_molecules(molecules=[ligand])
    ligand_omm_topology = ligand_off_topology.to_openmm()
    ligand_positions = offquantity_to_openmm(ligand.conformers[0])

    return ligand_omm_topology, ligand_positions


def select_platform(platform_name=None):
    """
    Select OpenMM platform and set useful defaults.
    """
    if platform_name is None or platform_name == "fastest":
        platform_name = get_fastest_platform().getName()
        logging.info(f"{platform_name} is the fastest platform available.")

    platform = Platform.getPlatformByName(platform_name)

    if platform_name == "OpenCL":
        platform.setPropertyDefaultValue("Precision", "mixed")
        platform.setPropertyDefaultValue("DeviceIndex", "0")

    if platform_name == "CUDA":
        platform.setPropertyDefaultValue("DeterministicForces", "true")
        platform.setPropertyDefaultValue("CudaPrecision", "mixed")
        platform.setPropertyDefaultValue("CudaDeviceIndex", "0")

    return platform


# %%
def export_all_files(system, simulation, sys_name, suffix, modeller, forcefield):
    """
    Export solvated coordinates, serialized system, checkpoint, and Amber files.
    """
    logging.info(f"Exporting files for {suffix}...")
    final_positions = simulation.context.getState(getPositions=True).getPositions()

    # Save solvated PDB using final positions from context
    with open(f"{sys_name}/{sys_name}_solvated.pdb", "w") as f:
        PDBFile.writeFile(modeller.topology, final_positions, f, keepIds=True)

    # Save serialized system
    with open(f"{sys_name}/system.xml", "w") as output:
        output.write(XmlSerializer.serialize(system))

    simulation.saveCheckpoint(f"{sys_name}/{sys_name}_{suffix}.chk")

    # Rebuild for ParmEd export
    new_system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=nb_cutoff * nanometers,
        removeCMMotion=False,
        rigidWater=False,
        hydrogenMass=4 * amu,
    )
    parmed_sys = parmed.openmm.load_topology(
        modeller.getTopology(), new_system, final_positions
    )
    # save amber parmameters
    parmed_sys.save(f"{sys_name}/{sys_name}_solvated.prmtop", overwrite=True)
    parmed_sys.save(
        f"{sys_name}/{sys_name}_solvated.rst7",
        format="rst7",
        overwrite=True,
    )
    # save gromacs parameters
    parmed_sys.save(f"{sys_name}/{sys_name}_solvated.gro", format='gro', overwrite=True)
    parmed_sys.save(f"{sys_name}/{sys_name}_solvated.top", format='gromacs', overwrite=True)


# %%
def main():
    logging.info("Preparing the receptor...")
    pdb_fixed = prepare_protein(
        args.rec,
        ignore_missing_residues=False,
        ignore_terminal_missing_residues=True,
        ph=7.4,
    )

    with open(f"{sys_name}/{rec_basename}_fixed.pdb", "w") as f:
        PDBFile.writeFile(pdb_fixed.topology, pdb_fixed.positions, f, keepIds=True)

    logging.info("Modelling the system...")

    # Create an OpenMM ForceField object with AMBER ff14SB and TIP3P
    forcefield = ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/tip3pfb.xml",
        "amber/tip3p_HFE_multivalent.xml",
    )

    # Make an OpenMM Modeller object with the protein
    modeller = Modeller(pdb_fixed.topology, pdb_fixed.positions)

    # Optional ligand
    if args.lig:
        logging.info(f"Ligand provided: {args.lig}")
        ligand = prepare_ligand(lig_sdf=args.lig, allow_undefined_stereo=True)

        ligand_topology, ligand_positions = parametrize_ligand(
            ligand,
            forcefield=forcefield,
            lig_ff=ligand_ff,
        )

        logging.info("Adding ligand to the modeller...")
        modeller.add(ligand_topology, ligand_positions)
    else:
        logging.info("No ligand provided. Running protein-only setup.")

    logging.info("Adding solvent and ions...")
    modeller.addSolvent(
        forcefield,
        ionicStrength=ionicStrength * molar,
        neutralize=True,
        boxShape=boxShape,
        padding=padding * nanometer,
    )

    logging.info("Selecting MD platform...")
    platform = select_platform("fastest")

    logging.info("Setting up the system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometers,
        removeCMMotion=False,
        rigidWater=True,
        hydrogenMass=hydrogenMass * amu,
        constraints=HBonds,
    )

    logging.info("Setting up the integrator...")
    integrator = LangevinMiddleIntegrator(
        300 * kelvin,
        1 / picoseconds,
        timestep * picoseconds,
    )

    logging.info("Setting up the simulation...")
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    export_all_files(
        system=system,
        simulation=simulation,
        sys_name=sys_name,
        suffix="solvated",
        modeller=modeller,
        forcefield=forcefield,
    )

    logging.info("All done!")


# %%
if __name__ == "__main__":
    main()
