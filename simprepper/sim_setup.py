from openmm import unit as mm_units
from openmm.unit.quantity import Quantity as mm_quantity
from dataclasses import dataclass, asdict
import argparse


# Units for different properties
QUANTITY_FIELDS = {
    "nb_cutoff": (mm_units.nanometer, "nm"),
    "hydrogenMass": (mm_units.amu, "amu"),
    "timestep": (mm_units.picoseconds, "ps"),
    "temperature": (mm_units.kelvin, "K"),
    "padding": (mm_units.nanometer, "nm"),
    "ionicStrength": (mm_units.molar, "M"),
}


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
    ligand_ff     : str         = "espaloma" # default = "espaloma"  # espaloma, SMIRNOFF, GAFF
    protein_ff    : str         = "amber14/protein.ff14SB.xml" # default = "amber14/protein.ff14SB.xml"
    water_ff      : str         = "amber14/tip3pfb.xml" # default = "amber14/tip3pfb.xml"
    ion_ff        : str         = "amber/tip3p_HFE_multivalent.xml" # default = "amber/tip3p_HFE_multivalent.xml"
    lipid_ff      : str         = "amber14/lipid17.xml" # default = "amber14/lipid17.xml"

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
                ligand_ff     = args.ligand_ff, # default = "espaloma"  # espaloma, SMIRNOFF, GAFF
                protein_ff    = args.protein_ff, # default = "amber14/protein.ff14SB.xml"
                water_ff      = args.water_ff, # default = "amber14/tip3pfb.xml"
                ion_ff        = args.ion_ff, # default = "amber/tip3p_HFE_multivalent.xml"
                lipid_ff      = args.lipid_ff, # default = "amber14/lipid17.xml"
        )
        #TODO: the following properties are currently fixed, because the argument parser doesn't know them
        # nb_cutoff
        # timestep
        # ionicStrength
        # ph
        # NOTE: not sure if we keep this method since we have the option to construct this class from an ini 
        # file, which is more flexible and easier to maintain when we have many parameters.

        return cls(**props)
    
    @classmethod
    def from_ini(cls, ini_fname):
        """
        This is a constructor function, to generate an instance of SimSetup() from an ini file.
        """
        raw = {}
        with open(ini_fname) as f:
            for line in f:
                # Remove inline comments
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
                key, value = [x.strip() for x in line.split("=", 1)]
                raw[key] = value

        props = {}
        for key, value in raw.items():
            if key in QUANTITY_FIELDS:
                unit_obj, _ = QUANTITY_FIELDS[key]
                props[key] = float(value) * unit_obj
            elif key == "ph":
                props[key] = float(value)
            elif key == "lig_fname":
                props[key] = None if value == "None" else value
            else:
                props[key] = value

        return cls(**props)
    
    def to_ini(self, out_fname):
        """
        Write the simulation setup to an ini file. This can be used to reload the same setup later, or to use it as a template for other setups.
        """
        props = asdict(self)
        key_length = get_longest_key(props.keys())

        with open(out_fname, "w") as f:
            for key, value in props.items():
                if key in QUANTITY_FIELDS:
                    # Add the unit of the quantity as an inline comment
                    unit_obj, unit_str = QUANTITY_FIELDS[key]
                    value = value.value_in_unit(unit_obj)
                    f.write(f"{key:<{key_length}} = {value} # {unit_str}\n")
                else:
                    f.write(f"{key:<{key_length}} = {value}\n")
        return None


def get_longest_key(keys):
    longest_string = max([len(key) for key in keys])
    return longest_string+1
