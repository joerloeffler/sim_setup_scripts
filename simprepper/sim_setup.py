from openmm import unit as mm_units
from openmm.unit.quantity import Quantity as mm_quantity
from dataclasses import dataclass, asdict
from configparser import ConfigParser
import argparse
import warnings

# %% global constants
# Units for different properties
QUANTITY_FIELDS = {
    "nb_cutoff": (mm_units.nanometer, "nm"),
    "hydrogenMass": (mm_units.amu, "amu"),
    "timestep": (mm_units.picoseconds, "ps"),
    "temperature": (mm_units.kelvin, "K"),
    "padding": (mm_units.nanometer, "nm"),
    "ionicStrength": (mm_units.molar, "M"),
}

# Sections for the config file
FIELD_SECTIONS = {
    "boxShape": "simulation box",
    "padding": "simulation box",

    "nb_cutoff": "simulation",
    "hydrogenMass": "simulation",
    "timestep": "simulation",
    "temperature": "simulation",
    "ionicStrength": "simulation",
    "ph": "simulation",

    "ligand_ff": "forcefields",
    "protein_ff": "forcefields",
    "water_ff": "forcefields",
    "ion_ff": "forcefields",
    "lipid_ff": "forcefields",
}

# %% SimSetup class
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
    ionicStrength : mm_quantity = 0.15  * mm_units.molar
    ph            : float       = 7.4
    ligand_ff     : str         = "GAFF" # default = "GAFF"  # espaloma, SMIRNOFF, GAFF
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
                ligand_ff     = args.ligand_ff, # default = "GAFF"  # espaloma, SMIRNOFF, GAFF
                protein_ff    = args.protein_ff, # default = "amber14/protein.ff14SB.xml"
                water_ff      = args.water_ff, # default = "amber14/tip3pfb.xml"
                ion_ff        = args.ion_ff, # default = "amber/tip3p_HFE_multivalent.xml"
                lipid_ff      = args.lipid_ff, # default = "amber14/lipid17.xml"
        )
        # NOTE: This method is still functional but it is a legacy feature, so it will not be futher developed.
        # Using the .ini file is the recommended way to go, as it allows for more flexibility and a wider range of options.

        return cls(**props)
    
    @classmethod
    def from_ini(cls, sys_name, args):
        """
        This is a constructor function, to generate an instance of SimSetup() from an ini file.
        """
        config = ConfigParser(inline_comment_prefixes=("#",))
        config.read(args.ini)

        # Print a warning if there are unknown sections in the ini file, to help users identify typos or misplaced options
        expected_sections = set(FIELD_SECTIONS.values())
        actual_sections = set(config.sections())
        unknown_sections = actual_sections - expected_sections
        if unknown_sections:
            warnings.warn(
                f"Unknown section(s) ignored: {', '.join(unknown_sections)}",
                UserWarning,
                stacklevel=2
            )

        defaults = cls()  # default fallback values from dataclass

        raw = {}
        for key, section in FIELD_SECTIONS.items():
            section_data = config[section] if section in config else {}
            value = section_data.get(key, None)
            # fallback to default
            if value is None:
                value = getattr(defaults, key)
            raw[key] = value

        # Set the receptor and ligand file names from the command line arguments
        props = dict(
            sys_name=sys_name,
            rec_fname=args.rec,
            lig_fname=args.lig,
        )
        # Add the other properties from the ini file
        for key, value in raw.items():
            if key in QUANTITY_FIELDS:
                unit_obj, _ = QUANTITY_FIELDS[key]
                if isinstance(value, mm_units.Quantity):
                    props[key] = value
                else:
                    props[key] = float(value) * unit_obj
            elif key == "ph":
                props[key] = float(value)
            else:
                props[key] = value

        return cls(**props)
    
    def to_ini(self, out_fname):
        """
        Write the simulation setup to an ini file. This can be used to reload the same setup later, or to use it as a template for other setups.
        """
        props = asdict(self)
        config = ConfigParser()

        for key, value in props.items():
            if key not in FIELD_SECTIONS:
                continue
            section = FIELD_SECTIONS[key]
            if section not in config:
                config[section] = {}
            if key in QUANTITY_FIELDS:
                unit_obj, unit_str = QUANTITY_FIELDS[key]
                value = value.value_in_unit(unit_obj)
                config[section][key] = f"{value}  # {unit_str}"
            else:
                config[section][key] = str(value)

        with open(out_fname, "w") as f:
            config.write(f)
        return None


# NOTE: this function is not actually used anywhere in the code anymore
def get_longest_key(keys):
    longest_string = max([len(key) for key in keys])
    return longest_string+1
