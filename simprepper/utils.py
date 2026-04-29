import openmm
import os
import logging
from openmmtools.utils import get_fastest_platform

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


