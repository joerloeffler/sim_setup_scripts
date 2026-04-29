#!/usr/bin/env python
# general imports
import os
import sys
import logging

# %% Hacks to get rid of annoying warning messages from imported modules...

os.environ["JAX_ENABLE_X64"] = "True" # get rid of annoying JAX warning
logging.getLogger("pymbar").setLevel(logging.ERROR)

# %% other imports

# import numpy as np
from typing import NamedTuple
from simprepper.argument_parsing import parser
from simprepper.utils import select_platform, get_sysname, prep_filetree, export_all_files
from simprepper.structure_prep import prepare_ligand, prepare_protein, parametrize_ligand

# OpenMM imports
import openmm
from openmm import app as mm_apps
from openmm import unit as mm_units

# %% CONSTANTS
LOG_PATH = 'prot_prep_logs'

# %% Pairwise distances...
#NOTE: PAQ: pw_dist is not actually used anywhere in this module!
#from simprepper.pairwise_distance import pw_dist

# %% setting up process

args = parser.parse_args()
sys_name, rec_basename = get_sysname(args)
prep_filetree(sys_name, log_path=LOG_PATH)
if args.debug:
    print(args)
    print("Log-level: {}".format(args.log_level.upper()))

logging.basicConfig(
    level=args.log_level.upper(),
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y/%m/%d %H:%M:%S",
    handlers=[
        logging.FileHandler(f"{LOG_PATH}/{sys_name}.log", mode="w"),
        logging.StreamHandler(sys.stdout),
    ],
    force=True  # NOTE: paq: otherwise doesn't print to stdout on all systems...
)

# TODO: store these setup parameters in an actual class?
class SimSetup(NamedTuple):
    sys_name      = sys_name
    nb_cutoff     = 1.0 * mm_units.nanometers    # TODO: Make this an input argument?
    hydrogenMass  = args.Hmass * mm_units.amu  # default =4
    timestep      = 0.004 * mm_units.picoseconds   # picoseconds
    temperature   = args.temperature * mm_units.kelvin
    boxShape      = args.box_shape # cube, dodecahedron
    padding       = args.box_padding * mm_units.nanometer 
    ionicStrength = 0.15 * mm_units.molar  # TODO: Make this an input argument?
    ph            = 7.4  # TODO: Make this an input argument?
    ligand_ff     = args.ligand_ff # default= "espaloma"  # espaloma, SMIRNOFF, GAFF
setup = SimSetup()


# %%
def main():

    logging.info("Preparing the receptor...")
    pdb_fixed = prepare_protein(args.rec,
                                ignore_missing_residues=False,
                                ignore_terminal_missing_residues=True,
                                ph=setup.ph,
                                )

    with open(f"{sys_name}/{rec_basename}_fixed.pdb", "w") as f:
        mm_apps.PDBFile.writeFile(pdb_fixed.topology, pdb_fixed.positions, f, keepIds=True)

    logging.info("Modelling the system...")

    # Create an OpenMM ForceField object with AMBER ff14SB and TIP3P
    forcefield = mm_apps.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/tip3pfb.xml",
        "amber/tip3p_HFE_multivalent.xml",
    )

    # Make an OpenMM Modeller object with the protein
    modeller = mm_apps.Modeller(pdb_fixed.topology, pdb_fixed.positions)

    # Optional ligand
    if args.lig:
        logging.info(f"Ligand provided: {args.lig}")
        ligand = prepare_ligand(lig_sdf=args.lig, allow_undefined_stereo=True)

        ligand_topology, ligand_positions = parametrize_ligand(ligand,
                                                               forcefield=forcefield,
                                                               lig_ff=setup.ligand_ff,
                                                               )

        logging.info("Adding ligand to the modeller...")
        modeller.add(ligand_topology, ligand_positions)
    else:
        logging.info("No ligand provided. Running protein-only setup.")

    logging.info("Adding solvent and ions...")
    modeller.addSolvent(
        forcefield,
        ionicStrength=setup.ionicStrength,
        neutralize=True,
        boxShape=setup.boxShape,
        padding=setup.padding,
    )

    logging.info("Selecting MD platform...")
    platform = select_platform("fastest")

    logging.info("Setting up the system...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=mm_apps.PME,
        nonbondedCutoff=setup.nb_cutoff,
        removeCMMotion=False,
        rigidWater=True,
        hydrogenMass=setup.hydrogenMass,
        constraints=mm_apps.HBonds,
    )

    logging.info("Setting up the integrator...")
    integrator = openmm.LangevinMiddleIntegrator(setup.temperature,
                                                 1 / mm_units.picoseconds,
                                                 setup.timestep,
                                                 )

    logging.info("Setting up the simulation...")
    simulation = mm_apps.Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    export_all_files(
        system=system,
        simulation=simulation,
        setup=setup,  # contains all kind of information, e.g. `sys_name`
        suffix="solvated",
        modeller=modeller,
        forcefield=forcefield,
    )

    logging.info("All done!")


# %%
if __name__ == "__main__":
    
    main()
