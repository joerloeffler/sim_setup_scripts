import openmm
import os
import parmed
import logging
from pathlib import Path

from openmmtools.utils import get_fastest_platform
from openmm import app as mm_apps
from openmm.app import forcefield


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


def get_basename(filename):
    """strip away trailing foldernames and filename extension
    """
    return os.path.splitext(os.path.basename(filename))[0]    


def get_sysname(args):
    # Use ligand name as system name if present, otherwise receptor basename
    rec_basename = get_basename(args.rec)
    sys_name = rec_basename

    if args.lig:
        lig_basename = get_basename(args.lig)
        sys_name = f"{rec_basename}-{lig_basename}"

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


def sanity_check_pdb_for_TERs(pdb_filename, verbose=False):
    """
    Checks for the presence of TER records in a PDB file 
    and prints the context of where they appear.
    """
    ter_count = 0
    last_atom_line = ""
    if verbose:
        print(f"Scanning '{pdb_filename}' for TER records...\n")
    
    with open(pdb_filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Keep track of the last ATOM/HETATM line seen before a TER
            if line.startswith(("ATOM  ", "HETATM")):
                last_atom_line = line.strip()
                
            elif line.startswith("TER"):
                ter_count += 1
                if verbose:
                    print(f" Found TER record at line {line_num}!")
                if last_atom_line:
                    # Parse out the residue info from the preceding atom line for context
                    res_name = last_atom_line[17:20].strip()
                    chain_id = last_atom_line[21].strip()
                    res_num = last_atom_line[22:26].strip()
                    if verbose:
                        print(f"   -> Placed after: {res_name} (Chain {chain_id}, Res #{res_num})")
                else:
                    if verbose:
                        print("   -> Placed at the very beginning of the file.")
    if verbose:            
        print("-" * 50)
        
    return ter_count > 0


def print_tree(path: Path, prefix: str = ""):
    """
    Recursively print a tree view of a directory.
    """
    if not path.exists():
        return

    # keep only XML files + directories that may contain XMLs
    entries = [
        p for p in path.iterdir()
        if p.is_dir() or p.suffix == ".xml"
    ]
    entries = sorted(entries, key=lambda p: (not p.is_dir(), p.name.lower()))

    for i, entry in enumerate(entries):
        connector = "`-- " if i == len(entries) - 1 else "|-- "
        print(prefix + connector + entry.name)

        if entry.is_dir():
            extension = "    " if i == len(entries) - 1 else "|   "
            print_tree(entry, prefix + extension)


def find_forcefields():
    """
    Lists available OpenMM forcefields using the _getDataDirectories() method.
    """
    data_dirs = forcefield._getDataDirectories()

    print("\nAvailable OpenMM forcefields:\n")

    for d in data_dirs:
        d = Path(d)
        print(f"{d.name}/")
        print_tree(d)
        print()