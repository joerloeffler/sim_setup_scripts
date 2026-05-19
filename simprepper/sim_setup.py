from openmm import unit as mm_units
from openmm.unit.quantity import Quantity as mm_quantity
from dataclasses import dataclass, asdict
import argparse


@dataclass
class SimSetup:
    """For developers: 
    if you want to add a property to this class:
    1. Add it here, to be a field that can be initialized
    2. make sure that it is passed on, in the constructor class `.from_args()`

    Note: you can define default values here, but they can be overwritten during construction (by other default values)
    So defining default values is not actually necessary. 
    HOWEVER, if you do not define a default value, it will be considered a non-default argument, which cannot follow a default argument...
    """
    sys_name      : str         = "DEFAULT_SYS_NAME"
    rec_fname     : str         = "DEFAULT_REC_NAME"
    lig_fname     : str | None  = None
    nb_cutoff     : mm_quantity = 1.0   * mm_units.nanometers
    hydrogenMass  : mm_quantity = 4     * mm_units.amu  # default =4
    timestep      : mm_quantity = 0.004 * mm_units.picoseconds   # picoseconds
    temperature   : mm_quantity = 300.0 * mm_units.kelvin
    boxShape      : str         = "cube" # cube, dodecahedron
    padding       : mm_quantity = 3.0   * mm_units.nanometer 
    ionicStrength : mm_quantity = 0.15  * mm_units.molar  # TODO: Make this an input argument?
    ph            : float       = 7.4   # TODO: Make this an input argument?
    ligand_ff     : str         = "espaloma" # default= "espaloma"  # espaloma, SMIRNOFF, GAFF
    protein_ff    : str         = "amber14/protein.ff14SB.xml" # NOTE: this is currently not used in the tool, but only printed to the out-file
    water_ff      : str         = "amber14/tip3pfb.xml" # NOTE: this is currently not used in the tool, but only printed to the out-file

    @classmethod
    def from_args(cls,
                 sys_name: str,
                 args: argparse.Namespace  # argparse arguments
                 ) -> None:
        """This is a constructor function, to generate an instance of SimSetup() from args
        """
        # what happens here?
        # we generate a dictionary from the `args` which is used to generate an instance of SimSetup
        # in the process of generating these dictionaries, we also add the corresponding units 
        props = dict(
                sys_name      = sys_name,
                rec_fname     = args.rec,
                lig_fname     = args.lig,
                hydrogenMass  = args.Hmass * mm_units.amu,  # default =4
                temperature   = args.temperature * mm_units.kelvin,
                boxShape      = args.box_shape, # cube, dodecahedron
                padding       = args.box_padding * mm_units.nanometer, 
                ligand_ff     = args.ligand_ff, # default= "espaloma"  # espaloma, SMIRNOFF, GAFF
        )
        #TODO: the following properties are currently fixed, because the argument parser doesn't know them
        # nb_cutoff
        # timestep
        # ionicStrength
        # ph

        return cls(**props)
    
    @classmethod
    def from_ini(cls, ini_fname):
        print("Reading the simulation config from an ini-file is currently not implemented, yet.")
        pass
    
    def to_ini(self, out_fname):
        props = asdict(self)
        key_length = get_longest_key(props.keys())

        with open(out_fname, "w") as f:
            for key,val in props.items():
                f.write(f"{key:<{key_length}} = {val}\n")
        return None


def get_longest_key(keys):
    longest_string = max([len(key) for key in keys])
    return longest_string+1
