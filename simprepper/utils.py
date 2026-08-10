import openmm
import os
import parmed
import logging
import argparse
import warnings

from pathlib import Path
from dataclasses import dataclass, field
from openmmtools.utils import get_fastest_platform
from openmm import app as mm_apps
from openmm.app import forcefield
from openmm import unit as mm_units
from simprepper.sim_setup import SimSetup


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


def prep_filetree(subdirs, subsubdirs):
    # sys_name,logpath = sdirs
    for sdir in subdirs:
        os.makedirs(sdir, exist_ok=True)
    for ssdir in subsubdirs:
        os.makedirs(ssdir, exist_ok=True)
    return None


def save_parmed(modeller, setup, forcefield_of_this_system, 
                final_positions, export_path_manager):
    # just a reminder:
    # fname_trunc = f"{sys_name}/{sys_name}"

    # Rebuild for ParmEd export
    new_system = forcefield_of_this_system.createSystem(
                                modeller.topology,
                                nonbondedMethod=mm_apps.PME,
                                nonbondedCutoff=setup.nb_cutoff,
                                removeCMMotion=False,
                                rigidWater=False,
                                hydrogenMass=setup.hydrogen_mass,
                                )
    parmed_sys = parmed.openmm.load_topology(modeller.getTopology(), 
                                             new_system, 
                                             final_positions
                                             )

    # save amber parmameters
    if export_path_manager.should_export_amber:
        parmed_sys.save(export_path_manager.output_fnames["amber"]["prmtop"], 
                        overwrite=True, format="amber")
        parmed_sys.save(export_path_manager.output_fnames["amber"]["coord_file"], 
                        overwrite=True, format="rst7")
        
    # save gromacs parameters
    if export_path_manager.should_export_gmx:
        #TODO: add gmx subfolder here
        parmed_sys.save(export_path_manager.output_fnames["gmx"]["gro"], 
                        overwrite=True, format='gro')
        
        parmed_sys.save(export_path_manager.output_fnames["gmx"]["top"], 
                        overwrite=True, format='gromacs')
    return None


# %%
def export_all_files(system, simulation, setup, modeller, forcefield_obj, 
                     export_path_manager):
    """
    Export solvated coordinates, serialized system, checkpoint, and Amber files.
    """
    
    logging.info(f"Exporting files for {setup.sys_name}, for the following engines:\n{export_path_manager.printout_flags}")
    # most files follow this naming convention:
    fname_trunc = f"{setup.sys_name}"

    # These positions will be used in all output
    final_positions = simulation.context.getState(getPositions=True).getPositions()

    # Save solvated PDB using final positions from context
    mm_apps.PDBFile.writeFile(modeller.topology, final_positions, 
                              f"{setup.sys_name}/{fname_trunc}_solvated.pdb", 
                              keepIds=True)

    if export_path_manager.should_export_openmm:
        #TODO: add openmm subfolder here
        # Save serialized system
        with open(export_path_manager.output_fnames["openmm"]["xml"], "w") as output:
            output.write(openmm.XmlSerializer.serialize(system))

        simulation.saveCheckpoint(export_path_manager.output_fnames["openmm"]["chk"])

    save_parmed(modeller, setup, forcefield_obj, 
                final_positions, export_path_manager)
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
        
    has_TERS =  ter_count > 0
    if not has_TERS:
            warnings.warn("\n".join(["Did not find any >TER< entry in your pdb-file.",
                                     "This *may* cause problems (depending on how you cap your protein chain).",
                                     "We recommend to add >TER< entries between all chains, and in particular between receptor and ligand"]),
                                     UserWarning)
    return None


def sanity_check_ligand_extension(lig_filename):
    """
    Checks if the ligand file has a supported extension.
    """
    # List of currently supported ligand file extensions
    supported_extensions = [".sdf"]

    lig_path = Path(lig_filename)

    # Check if the ligand file exisits
    if not lig_path.exists():
        raise FileNotFoundError(f"Ligand file '{lig_path}' does not exist.")
    
    # Get the file extension and check if it's in the list of supported extensions
    lig_suffix = lig_path.suffix.lower()
    if lig_suffix not in supported_extensions:
        raise ValueError(f"Unsupported ligand file extension '{lig_suffix}'. Supported extensions are: {', '.join(supported_extensions)}")
    

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


def write_example_ini():
    """
    Writes an example ini file with all default values.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--output",
        default="default_config.ini",
        help="Output .ini filename"
    )
    parser.add_argument(
        "--membrane",
        action="store_true",
        help="Generate a configuration template for a membrane protein system"
    )
    args = parser.parse_args()

    SimSetup().to_ini(args.output, args.membrane)

    print(f"Wrote example config to {args.output}")


def calculate_protein_dimensions(positions):
    """
    Calculate the bounding dimensions (X, Y, Z) of the protein 
    based on its atomic coordinates array (in nanometers).
    """
    min_coords = positions.min(axis=0)
    max_coords = positions.max(axis=0)
    return max_coords - min_coords


def check_initial_box_dimensions(protein_dimensions, box_vectors):
    """
    Check if the initial periodic box dimensions are sufficient to accommodate the protein.
    """
    axis_names = ['X', 'Y', 'Z']
    
    # Extract periodic box diagonal lengths in nanometers
    box_lengths = [
        box_vectors[i][i].value_in_unit(mm_units.nanometers) 
        if hasattr(box_vectors[i][i], 'value_in_unit') 
        else box_vectors[i][i] 
        for i in range(3)
    ]
    
    # Check each axis against protein extent
    for i in range(3):
        if protein_dimensions[i] > box_lengths[i]:
            logging.warning(
                f"Initial box dimension {axis_names[i]} ({box_lengths[i]:.2f} nm) "
                f"is smaller than protein span ({protein_dimensions[i]:.2f} nm). "
                f"OpenMM will and expand the unit cell during membrane/solvent addition. " 
                f"Make sure to inspect if the final box dimensions are sufficient to accommodate the protein and any added solvent/ions."
            )
    print(f"Wrote example config to {args.output}")


@dataclass
class ExportPathManager:
    """
    Tells the simprepper, which output should be generated:
    Gromacs? Yes/No; Amber? Yes/No...
    Small class that care of parsing some arguments, and conveneiently
    """
    sys_name             : str         = "DEFAULT_SYS_NAME"
    should_export_gmx    : bool        = False
    should_export_amber  : bool        = False
    should_export_openmm : bool        = False
    # the subsubdirectories-pathnames will be None by default
    gmx_ssdir            : Path | None = None
    amber_ssdir          : Path | None = None
    openmm_ssdir         : Path | None = None
    output_fnames        : dict        = field(default_factory=dict)

    printout_flags       : str         = ""

    #TODO: Add another constructor, to be used from jupyter notebooks.

    @classmethod
    def from_args(cls, 
                  sys_name: str, args: argparse.Namespace  # argparse arguments
                  ):
        """This is a constructor function, to generate an instance of ExportPathManager() from argparse-arguments
        """
        props = dict(sys_name             = sys_name,
                     should_export_gmx    = args.should_export_gmx,
                     should_export_amber  = args.should_export_amber,
                     should_export_openmm = args.should_export_openmm
                     )
        # if requested, add subsubdirectories:
        output_filenames = {}
        printout_flags = []
        if args.should_export_amber:
            amber_ssdir = Path(sys_name) / Path("amber")
            props["amber_ssdir"] = amber_ssdir
            output_filenames["amber"] = {"prmtop"     : f"{amber_ssdir}/{sys_name}_solvated.prmtop",
                                         "coord_file" : f"{amber_ssdir}/{sys_name}_solvated.rst7"}
            printout_flags.append("Amber")
            
        if args.should_export_gmx:
            gmx_ssdir = Path(sys_name) / Path("gmx")
            props["gmx_ssdir"] = gmx_ssdir
            output_filenames["gmx"]  = {"top" : f"{gmx_ssdir}/{sys_name}_solvated.top",
                                        "gro" : f"{gmx_ssdir}/{sys_name}_solvated.gro"}
            printout_flags.append("Gromacs")
            
        if args.should_export_openmm:
            openmm_ssdir = Path(sys_name) / Path("openmm")
            props["openmm_ssdir"] = openmm_ssdir
            output_filenames["openmm"]  = {"xml" : f"{openmm_ssdir}/{sys_name}_system.xml",
                                           "chk" : f"{openmm_ssdir}/{sys_name}_solvated.chk"}
            printout_flags.append("OpenMM")

        props["output_fnames"] = output_filenames
        props["printout_flags"] = ", ".join(printout_flags)
        return cls(**props)
    
    
    def get_subsubdirectories(self):
        ssdirs = []
        if self.should_export_gmx:
            ssdirs.append(self.gmx_ssdir)
        if self.should_export_amber:
            ssdirs.append(self.amber_ssdir)
        if self.should_export_openmm:
            ssdirs.append(self.openmm_ssdir)
        return ssdirs
