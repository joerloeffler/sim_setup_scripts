import openmm
import os
import parmed
import logging

from openmmtools.utils import get_fastest_platform
from openmm import app as mm_apps


def select_platform(platform_name=None):
    """
    Select OpenMM platform and set useful defaults.
    """
    if platform_name is None or platform_name == "fastest":
        platform_name = get_fastest_platform().getName()
        logging.info(f"{platform_name} is the fastest platform available.")

    platform = openmm.Platform.getPlatformByName(platform_name)

    if platform_name == "OpenCL":
        platform.setPropertyDefaultValue("Precision", "mixed")
        platform.setPropertyDefaultValue("DeviceIndex", "0")

    if platform_name == "CUDA":
        platform.setPropertyDefaultValue("DeterministicForces", "true")
        platform.setPropertyDefaultValue("CudaPrecision", "mixed")
        platform.setPropertyDefaultValue("CudaDeviceIndex", "0")

    return platform


def get_sysname(args):
    # Use ligand name as system name if present, otherwise receptor basename
    if args.lig:
        sys_name = os.path.splitext(os.path.basename(args.lig))[0]
    else:
        sys_name = os.path.splitext(os.path.basename(args.rec))[0]

    rec_basename = os.path.splitext(os.path.basename(args.rec))[0]
    return sys_name, rec_basename


def prep_filetree(sys_name, log_path):
    os.makedirs(sys_name, exist_ok=True)
    os.makedirs(log_path, exist_ok=True)
    return None


def save_parmed(parmed_sys, fname_trunc, should_save_amber=True, should_save_gmx=True):
    # save amber parmameters
    # fname_trunc = f"{sys_name}/{sys_name}"
    if should_save_amber:
        parmed_sys.save(f"{fname_trunc}_solvated.prmtop", overwrite=True)
        parmed_sys.save(f"{fname_trunc}_solvated.rst7", overwrite=True,
                        format="rst7",
                        )
    # save gromacs parameters
    if should_save_gmx:
        parmed_sys.save(f"{fname_trunc}_solvated.gro", overwrite=True,
                        format='gro',)
        parmed_sys.save(f"{fname_trunc}_solvated.top", overwrite=True,
                        format='gromacs',)
    return None



# %%
def export_all_files(system, simulation, setup, suffix, modeller, forcefield):
    """
    Export solvated coordinates, serialized system, checkpoint, and Amber files.
    """
    logging.info(f"Exporting files for {suffix}...")
    final_positions = simulation.context.getState(getPositions=True).getPositions()
    # most files follow this naming convention:
    fname_trunc = f"{setup.sys_name}/{setup.sys_name}"
    # Save solvated PDB using final positions from context
    #with open(f"{fname_trunc}_solvated.pdb", "w") as fhandle:
    mm_apps.PDBFile.writeFile(modeller.topology, final_positions, 
                                f"{fname_trunc}_solvated.pdb", 
                                keepIds=True)

    # Save serialized system
    with open(f"{setup.sys_name}/system.xml", "w") as output:
        output.write(openmm.XmlSerializer.serialize(system))

    simulation.saveCheckpoint(f"{fname_trunc}_{suffix}.chk")

    # Rebuild for ParmEd export
    new_system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=mm_apps.PME,
        nonbondedCutoff=setup.nb_cutoff,
        removeCMMotion=False,
        rigidWater=False,
        hydrogenMass=setup.hydrogenMass,
    )
    parmed_sys = parmed.openmm.load_topology(
        modeller.getTopology(), new_system, final_positions
    )
    save_parmed(parmed_sys, fname_trunc)
    return None
