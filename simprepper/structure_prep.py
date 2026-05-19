import pdbfixer
# OpenFF-toolkit imports
from openff.units.openmm import to_openmm as off_to_openmm
from openff.toolkit import Molecule
from openff.toolkit import Topology as offTopology

# RDKit imports
from rdkit.Chem import SDMolSupplier

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

    ligand = Molecule.from_rdkit(rdkit_mol,
                                 allow_undefined_stereo=allow_undefined_stereo,
                                 )
    return ligand


# %%
def parametrize_ligand(ligand, lig_ff="espaloma"):
    """
    Parameterize ligand and register template generator with the forcefield.
    """
    lig_ff_lower = lig_ff.lower()

    if lig_ff_lower == "espaloma":
        from openmmforcefields.generators import EspalomaTemplateGenerator

        ff_template = EspalomaTemplateGenerator(molecules=ligand, 
                                                forcefield="espaloma-0.3.1"
                                                )
    elif lig_ff_lower == "smirnoff":
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator

        ff_template = SMIRNOFFTemplateGenerator(molecules=ligand, 
                                                forcefield="openff-1.2.0"
                                                )
    elif lig_ff_lower == "gaff":
        from openmmforcefields.generators import GAFFTemplateGenerator

        ff_template = GAFFTemplateGenerator(molecules=ligand, 
                                            forcefield="gaff-2.11"
                                            )
    else:
        raise ValueError(
            "Ligand forcefield must be one of: espaloma, SMIRNOFF, or GAFF"
        )

    ligand_off_topology = offTopology.from_molecules(molecules=[ligand])
    ligand_omm_topology = ligand_off_topology.to_openmm()
    ligand_positions = off_to_openmm(ligand.conformers[0])

    return ligand_omm_topology, ligand_positions, ff_template.generator
