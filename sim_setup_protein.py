# %%
# general imports
import argparse
import os
import pdbfixer
import parmed
import logging
import numpy as np

# OpenMM imports
from openmm import *
from openmm.app import *
from openmm.unit import *

from openmmtools.utils import get_fastest_platform

# cython-optimized pairwise distance function
try:
    # add to PYTHONPATH current workdir and script dir
    sys.path.insert(0, os.curdir)
    sys.path.insert(0, os.path.dirname(__file__))

    from _pwdistance import pw_dist
except ImportError as e:
    logging.warning('Using numpy (slower) pwdist routine for simulation setup.')
    from scipy.spatial.distance import pdist

    def pw_dist(xyz_array):
        return np.amax(pdist(xyz_array, 'euclidean'))
# %%
parser = argparse.ArgumentParser(prog='openMM_prepare', description='Script to prepare OpenMM system starting from protein PDB file')
parser.add_argument('-r','--rec', help='PDB file of the receptor', required=True)
parser.add_argument('-L','--log-level', help='Choose the logging level to show', choices=['debug', 'info', 'warning', 'error', 'critical'], default='info', required=False)

args = parser.parse_args()

sys_name = os.path.splitext(os.path.basename(args.rec))[0]

os.makedirs(sys_name, exist_ok=True)
os.makedirs('logs', exist_ok=True)

logging.basicConfig(
level=args.log_level.upper(),
format="%(asctime)s [%(levelname)s] %(message)s",
handlers=[
    logging.FileHandler(f"logs/{sys_name}.log", mode='w'),
    logging.StreamHandler()
])

nb_cutoff = 1.0 # nanometers
hydrogenMass = 1 # amu
timestep=0.004 # pico
boxShape = 'cube' # cube, dodecahedron
padding=2 # nanometers
ionicStrength=0.15 # molar concentration

#%%
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
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens

    return fixer

# %%
def select_platform(platform_name: str=None):

    if platform_name is None or platform_name == 'fastest':
        platform_name = get_fastest_platform().getName()
        logging.info(f'{platform_name} is the fastest platform available.')

    platform = Platform.getPlatformByName(platform_name)
    
    if platform_name in ['OpenCL']:
        platform.setPropertyDefaultValue('Precision', 'mixed')
        platform.setPropertyDefaultValue('DeviceIndex', '0')
    if platform_name in ['CUDA']:
        platform.setPropertyDefaultValue('DeterministicForces', 'true')
        platform.setPropertyDefaultValue('CudaPrecision', 'mixed')
        platform.setPropertyDefaultValue('CudaDeviceIndex', '0')

    return platform

# %%
def export_all_files(system, simulation, sys_name, suffix):

    logging.info(f"Exporting files for {suffix}...")
    final_positions = simulation.context.getState(getPositions=True).getPositions()

    # Save solvated PDB
    PDBFile.writeFile(modeller.topology, modeller.positions, open(f"{sys_name}/{sys_name}_solvated.pdb", 'w'))

    # Save a serialized version of the system. This stores the forcefield parameters.
    with open(f'{sys_name}/system.xml', 'w') as output:
        output.write(XmlSerializer.serialize(system))

    simulation.saveCheckpoint(f"{sys_name}/{sys_name}_{suffix}.chk")

    new_system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=nb_cutoff*nanometers,
                                    removeCMMotion=False, rigidWater=False, hydrogenMass=1*amu) # adjustable usage of HMR
    parmed_sys = parmed.openmm.load_topology(modeller.getTopology(), new_system, final_positions)
    # save amber parameters
    parmed_sys.save(f"{sys_name}/{sys_name}_solvated.prmtop", overwrite=True)
    parmed_sys.save(f"{sys_name}/{sys_name}_solvated.rst7", format='rst7', overwrite=True)
    # save gromacs parameters
    parmed_sys.save(f"{sys_name}/{sys_name}_solvated.gro", format='gro', overwrite=True)
    parmed_sys.save(f"{sys_name}/{sys_name}_solvated.top", format='gromacs', overwrite=True)

    return

# %%
# prepare the receptor and save it
logging.info('Preparing the receptor..')
pdb_fixed = prepare_protein(args.rec, ignore_missing_residues=False, ignore_terminal_missing_residues=True, ph=7.4) # adjustable pH in case this becomes relevant
PDBFile.writeFile(pdb_fixed.topology, pdb_fixed.positions, open(f'{sys_name}/{sys_name}_fixed.pdb', 'w'))

logging.info('Modelling the system..')

# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P
forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml', 'amber/tip3p_HFE_multivalent.xml') # current water and forcefield model used in the setup

# make an OpenMM Modeller object with the protein
modeller = Modeller(pdb_fixed.topology, pdb_fixed.positions)

# Add solvent to the system
modeller.addSolvent(forcefield, ionicStrength=ionicStrength*molar, neutralize=True, 
                    boxShape=boxShape, padding=padding*nanometer)

logging.info('Selecting MD platform..')
platform = select_platform('fastest')

logging.info('Setting up the system..')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers,
                                 removeCMMotion=False, rigidWater=True, hydrogenMass=hydrogenMass*amu,
                                 constraints=HBonds) # 1.5 is OpenMM default and should work with 4fs, Peter Eastman said

logging.info('Setting up the integrator..')
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picoseconds, timestep*picoseconds) # Set up the integrator

logging.info('Setting up the simulation..')
simulation = Simulation(modeller.topology, system, integrator, platform) # Set up the simulation
simulation.context.setPositions(modeller.positions) # set the positions

export_all_files(system=system, simulation=simulation, sys_name=sys_name, suffix='solvated')

logging.info('All done!')
