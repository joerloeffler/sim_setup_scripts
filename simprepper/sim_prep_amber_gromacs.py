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

import numpy as np
from simprepper.argument_parsing import parser
from simprepper.utils.utils import select_platform, get_sysname, prep_filetree, export_all_files, sanity_check_pdb_for_TERs, sanity_check_ligand_extension, check_initial_box_dimensions, calculate_protein_dimensions
from simprepper.utils.utils import ExportPathManager
from simprepper.structure_prep import prepare_ligand, prepare_protein, parametrize_ligand
from simprepper.simprepper_setup import SimSetup
from simprepper.simprepper_setup import SimSetup
from simprepper.config_checker import ConfigChecker

# OpenMM imports
import openmm
from openmm import app as mm_apps
from openmm import unit as mm_units

# %% CONSTANTS
#TODO: make LOG_PATH a parsable argument?
LOG_PATH = 'simprepper_logs'
#TODO: Will we ever not want to save the setup?
SHOULD_SAVE_CONFIG = True

# %% Pairwise distances...
#NOTE: PAQ: pw_dist is not actually used anywhere in this module!
#from simprepper.pairwise_distance import pw_dist

# %% Process ligand
#TODO: Move this function somewhere else. Maybe utils...
def process_ligand(setup, verbose):
    '''Some sanity checks,
    prepare ligand structure, prepare forcefield_gerator
    '''
    # Check if a given ligand is in the correct format
    sanity_check_ligand_extension(setup.lig_fname)

    logging.info(f"Ligand provided: {setup.lig_fname}")

    # will throw a warning, if there are no TERs, but will continue in any case
    sanity_check_pdb_for_TERs(setup.rec_fname, verbose=verbose)

    ligand = prepare_ligand(lig_sdf=setup.lig_fname, 
                            allow_undefined_stereo=True)

    ligand_topology, ligand_positions, ligand_template_gen = parametrize_ligand(ligand,
                                                                                lig_ff=setup.config.ligand_ff
                                                                                )
    return ligand, ligand_template_gen


# %% setting up process

args = parser.parse_args()
sys_name, rec_basename = get_sysname(args)

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

# Construct this class-instance from ini file or from command line arguments
if args.ini is not None:
    setup = SimSetup.from_ini(sys_name, 
                              rec_fname=args.rec,
                              lig_fname=args.lig,
                              ini_fname=args.ini)
else:
    msg = ("Since you did not provide any .ini file, " + \
            "we are using default values for all settings. " + \
            "Consider executing `simprepper-example-config` to see all available settings.")
    logging.info(msg)
    setup = SimSetup.from_defaults(sys_name, 
                                   rec_fname=args.rec,
                                   lig_fname=args.lig)


export_path_manager = ExportPathManager.from_args(sys_name, args)
prep_filetree(subdirs=(sys_name, LOG_PATH),
              subsubdirs=export_path_manager.get_subsubdirectories())


# %% define main()
def main():
    # check for logical errors in the config, before anything else happens.
    config_checker = ConfigChecker(setup.config, logging, verbose=args.verbose)
    config_checker.run_all_checks() 
    config_checker.get_report()

    # the setup contains sys_name and file_names, but also the `SimPrepperConfig`.
    config = setup.config  # We mostly need the SimPrepperConfig fields, which are direct attributes
    logging.info("Preparing the receptor...")
    pdb_fixed = prepare_protein(setup.rec_fname,
                                ignore_missing_residues=False,
                                ignore_terminal_missing_residues=True,
                                ph=config.ph,
                                )

    
    mm_apps.PDBFile.writeFile(pdb_fixed.topology, 
                              pdb_fixed.positions, 
                              f"{sys_name}/{rec_basename}_fixed.pdb", 
                              keepIds=True
                              )

    logging.info("Modelling the system...")
    #BUG: If we pass `None` as an argument to mm_apps.ForceField it will look for a File called `None`
    # Therefore, None-forcefields, should not be passed on...
    # Create an OpenMM ForceField object with parsed forcefields
    if config.lipid_ff is None:
        forcefield = mm_apps.ForceField(
            config.protein_ff,
            config.water_ff,
            config.ion_ff
        )
        logging.info(f"Using {config.protein_ff}, {config.water_ff}, and {config.ion_ff} forcefields.")
    else:
        forcefield = mm_apps.ForceField(
            config.protein_ff,
            config.water_ff,
            config.ion_ff,
            config.lipid_ff,
        )
        logging.info(f"Using {config.protein_ff}, {config.water_ff}, {config.ion_ff}, and {config.lipid_ff} forcefields.")

    # Optional ligand
    if setup.lig_fname:
        #NOTE: `ligand` is currently not acually used anywhere anymore
        ligand, ligand_template_generator = process_ligand(setup, args.verbose)

        logging.info("Adding ligand to the modeller...")
        forcefield.registerTemplateGenerator(ligand_template_generator)
        
        # if you use the following line: the ligand is added a second time.
        # sys_modeller.add(ligand_topology, ligand_positions)  
        # -> ligand_topology, and ligand_positions are now hidden in `process_ligand`
        # I kept this line here, because in alternative preparation pipelines, this might become useful again...
    else:
        logging.info("No ligand provided. Running protein-only setup.")


    # Make an OpenMM Modeller object with the protein
    sys_modeller = mm_apps.Modeller(pdb_fixed.topology, 
                                    pdb_fixed.positions)
    
    #TODO: incorporate this in the SetupChecker?
    #   -> Maybe this is more of a problem for a new module,like a `ModelChecker`?
    check_initial_box_dimensions(sys_modeller)

    if config.membrane_protein:
        logging.info("Preparing membrane protein system. Adding membrane, solvent, and ions...")
        # TODO: find a way to define box size more precisely
        sys_modeller.addMembrane(forcefield,
                                 lipidType=config.lipid_type,
                                 membraneCenterZ=config.membrane_center_z,
                                 minimumPadding=config.minimum_padding,
                                 ionicStrength=config.ionic_strength,
                                 neutralize=True,
                                 )
    else:
        logging.info("Preparing soluble protein system. Adding solvent and ions...")
        # TODO: bring back an option of providing box_length instead of box_padding (was implemented with 
        # argument parsing, but got lost with the new config file approach)
        sys_modeller.addSolvent(forcefield,
                                ionicStrength=config.ionic_strength,
                                neutralize=True,
                                boxShape=config.box_shape,
                                padding=config.padding,
                                )
    logging.info("Selecting MD platform...")
    platform = select_platform("fastest")

    logging.info("Setting up the system...")
    system = forcefield.createSystem(sys_modeller.topology,
                                     nonbondedMethod=mm_apps.PME,
                                     nonbondedCutoff=config.nb_cutoff,
                                     removeCMMotion=False,
                                     rigidWater=True,
                                     hydrogenMass=config.hydrogen_mass,
                                     constraints=mm_apps.HBonds,
                                     )

    logging.info("Setting up the integrator...")
    integrator = openmm.LangevinMiddleIntegrator(config.temperature,
                                                 1 / mm_units.picoseconds,
                                                 config.timestep,
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
                     modeller=sys_modeller,
                     forcefield_obj=forcefield,
                     export_path_manager=export_path_manager,
                     )
    
    if SHOULD_SAVE_CONFIG:
        out_fname = os.path.join(LOG_PATH, "simprepper.out.ini")
        logging.info(f"Writing simulation setup to file {out_fname}.")
        setup.to_ini(out_fname, setup.config.membrane_protein)

    logging.info("All done!")
    return 1  # end of main()


# %% actually run main()
if __name__ == "__main__":
    
    main()
