#!/usr/bin/env python
# general imports
import os
import sys
import logging
import warnings

# %% Hacks to get rid of annoying warning messages from imported modules...

os.environ["JAX_ENABLE_X64"] = "True" # get rid of annoying JAX warning
logging.getLogger("pymbar").setLevel(logging.ERROR)

# %% other imports

# import numpy as np
from simprepper.argument_parsing import parser
from simprepper.utils import select_platform, get_sysname, prep_filetree, export_all_files, sanity_check_pdb_for_TERs
from simprepper.structure_prep import prepare_ligand, prepare_protein, parametrize_ligand
from simprepper.sim_setup import SimSetup

# OpenMM imports
import openmm
from openmm import app as mm_apps
from openmm import unit as mm_units

# %% CONSTANTS
#TODO: make LOG_PATH a parsable argument?
LOG_PATH = 'prot_prep_logs'
SHOULD_SAVE_SETUP = True

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

# Construct this class-instance from the parsed arguments.
setup = SimSetup.from_args(sys_name, args)


# %% define main()
def main():

    logging.info("Preparing the receptor...")
    pdb_fixed = prepare_protein(setup.rec_fname,
                                ignore_missing_residues=False,
                                ignore_terminal_missing_residues=True,
                                ph=setup.ph,
                                )

    
    mm_apps.PDBFile.writeFile(pdb_fixed.topology, 
                              pdb_fixed.positions, 
                              f"{sys_name}/{rec_basename}_fixed.pdb", 
                              keepIds=True
                              )

    logging.info("Modelling the system...")

    # Create an OpenMM ForceField object with AMBER ff14SB and TIP3P
    # TODO: make forcefield a parsable argument?
    forcefield = mm_apps.ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/tip3pfb.xml",
        "amber/tip3p_HFE_multivalent.xml",
    )

    # Make an OpenMM Modeller object with the protein
    sys_modeller = mm_apps.Modeller(pdb_fixed.topology, 
                                    pdb_fixed.positions)

    # Optional ligand
    if setup.lig_fname:
        logging.info(f"Ligand provided: {setup.lig_fname}")

        has_TERS = sanity_check_pdb_for_TERs(setup.rec_fname, verbose=args.verbose)
        if not has_TERS:
            warnings.warn("\n".join(["Did not find any >TER< entry in your pdb-file.",
                                     "This *may* cause problems (depending on how you cap your protein chain).",
                                     "We recommend to add >TER< entries between all chains, and in particular between receptor and ligand"]),
                                     UserWarning)

        ligand = prepare_ligand(lig_sdf=setup.lig_fname, 
                                allow_undefined_stereo=True)

        ligand_topology, ligand_positions, ligand_template_gen = parametrize_ligand(ligand,
                                                                                    lig_ff=setup.ligand_ff
                                                                                    )

        logging.info("Adding ligand to the modeller...")
        forcefield.registerTemplateGenerator(ligand_template_gen)
        
        # if you use the following line: the ligand is added a second time.
        # I will keep this line here, because in alternative preparation pipelines, this is a very useful line...
        # sys_modeller.add(ligand_topology, ligand_positions)
    else:
        logging.info("No ligand provided. Running protein-only setup.")

    logging.info("Adding solvent and ions...")
    sys_modeller.addSolvent(forcefield,
                            ionicStrength=setup.ionicStrength,
                            neutralize=True,
                            boxShape=setup.boxShape,
                            padding=setup.padding,
                            )

    logging.info("Selecting MD platform...")
    platform = select_platform("fastest")

    logging.info("Setting up the system...")
    system = forcefield.createSystem(sys_modeller.topology,
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
    simulation = mm_apps.Simulation(sys_modeller.topology, 
                                    system, 
                                    integrator, 
                                    platform)
    simulation.context.setPositions(sys_modeller.positions)

    export_all_files(system=system,
                     simulation=simulation,
                     setup=setup,  # contains all kind of information, e.g. `sys_name`
                     suffix="solvated",
                     modeller=sys_modeller,
                     forcefield=forcefield,
                     )
    
    if SHOULD_SAVE_SETUP:
        out_fname = os.path.join(LOG_PATH, "simprepper.out.ini")
        logging.info(f"Writing simulation setup to file {out_fname}.")
        setup.to_ini(out_fname)

    logging.info("All done!")


# %% actually run main()
if __name__ == "__main__":
    
    main()
