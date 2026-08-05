from openmm import unit as mm_units
from openmm.unit.quantity import Quantity as mm_quantity
from openff.units.openmm import openmm_unit_to_string
from dataclasses import dataclass, field, asdict
from collections import namedtuple
from configparser import ConfigParser
import argparse
import warnings
import os

# %% global constants
# Sections for the config file
# NOTE: At the moment, when adding a new parameter, you need to add it to the CONFIG_SECTIONS dictionary, 
# and also add it to the SimSetup dataclass. Ideally, we would like to have a single source of truth for 
# the parameters, but this can be solved later.
SimSetupField = namedtuple("SimSetupField", ["name", "unit", "default", "data_type"])
CONFIG_SECTIONS = {
    "system": [
        SimSetupField("membrane", None, False, bool),
    ],
    "membrane": [     # membrane parameters present only if membrane=True
        SimSetupField("lipid_type", None, "POPC", str),
        SimSetupField("membrane_center_z", None, 0.0, float),
        SimSetupField("minimum_padding", mm_units.nanometer, 1.0, float),
    ],
    "simulation": [
        SimSetupField("nb_cutoff", mm_units.nanometer, 1.0, float),
        SimSetupField("hydrogen_mass", mm_units.amu, 4.0, float),
        SimSetupField("timestep", mm_units.picoseconds, 0.004, float),
        SimSetupField("temperature", mm_units.kelvin, 300.0, float),
        SimSetupField("ionic_strength", mm_units.molar, 0.15, float),
        SimSetupField("ph", None, 7.4, float),
    ],
    "simulation box": [
        SimSetupField("box_shape", None, "cube", str),
        SimSetupField("padding", mm_units.nanometer, 3.0, float),
    ],
    "forcefields": [
        SimSetupField("ligand_ff", None, "GAFF", str),
        SimSetupField("protein_ff", None, "amber14/protein.ff14SB.xml", str),
        SimSetupField("water_ff", None, "amber14/tip3pfb.xml", str),
        SimSetupField("ion_ff", None, "amber/tip3p_HFE_multivalent.xml", str),
        SimSetupField("lipid_ff", None, None, str),
    ]
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
    sys_name         : str         = "DEFAULT_SYS_NAME"
    rec_fname        : str         = "DEFAULT_REC_NAME"
    lig_fname        : str | None  = None
    membrane         : bool        = field(default_factory=lambda: get_field_default("membrane", CONFIG_SECTIONS))
    lipid_type       : str         = field(default_factory=lambda: get_field_default("lipid_type", CONFIG_SECTIONS))
    membrane_center_z: float       = field(default_factory=lambda: get_field_default("membrane_center_z", CONFIG_SECTIONS))
    minimum_padding  : mm_quantity = field(default_factory=lambda: get_field_default("minimum_padding", CONFIG_SECTIONS))
    nb_cutoff        : mm_quantity = field(default_factory=lambda: get_field_default("nb_cutoff", CONFIG_SECTIONS))
    hydrogen_mass    : mm_quantity = field(default_factory=lambda: get_field_default("hydrogen_mass", CONFIG_SECTIONS))
    timestep         : mm_quantity = field(default_factory=lambda: get_field_default("timestep", CONFIG_SECTIONS))
    temperature      : mm_quantity = field(default_factory=lambda: get_field_default("temperature", CONFIG_SECTIONS))
    box_shape        : str         = field(default_factory=lambda: get_field_default("box_shape", CONFIG_SECTIONS))
    padding          : mm_quantity = field(default_factory=lambda: get_field_default("padding", CONFIG_SECTIONS))
    ionic_strength   : mm_quantity = field(default_factory=lambda: get_field_default("ionic_strength", CONFIG_SECTIONS))
    ph               : float       = field(default_factory=lambda: get_field_default("ph", CONFIG_SECTIONS))
    ligand_ff        : str         = field(default_factory=lambda: get_field_default("ligand_ff", CONFIG_SECTIONS))
    protein_ff       : str         = field(default_factory=lambda: get_field_default("protein_ff", CONFIG_SECTIONS))
    water_ff         : str         = field(default_factory=lambda: get_field_default("water_ff", CONFIG_SECTIONS))
    ion_ff           : str         = field(default_factory=lambda: get_field_default("ion_ff", CONFIG_SECTIONS))
    lipid_ff         : str | None  = field(default_factory=lambda: get_field_default("lipid_ff", CONFIG_SECTIONS))

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
                hydrogen_mass = args.Hmass * mm_units.amu,  # default =4
                temperature   = args.temperature * mm_units.kelvin,
                box_shape     = args.box_shape, # cube, dodecahedron
                padding       = args.box_padding * mm_units.nanometer, 
                ligand_ff     = args.ligand_ff, # default = "GAFF"  # espaloma, SMIRNOFF, GAFF
                protein_ff    = args.protein_ff, # default = "amber14/protein.ff14SB.xml"
                water_ff      = args.water_ff, # default = "amber14/tip3pfb.xml"
                ion_ff        = args.ion_ff, # default = "amber/tip3p_HFE_multivalent.xml"
                lipid_ff      = args.lipid_ff, # default = None
        )
        # NOTE: This method is still functional but it is a legacy feature, so it will not be futher developed.
        # Using the .ini file is the recommended way to go, as it allows for more flexibility and a wider range of options.

        return cls(**props)
    
    @classmethod
    def from_ini(cls, sys_name, args):
        """
        This is a constructor function, to generate an instance of SimSetup() from an ini file.
        """
        # Check if the ini file exists
        sanity_check_ini_file_if_exists(args.ini, CONFIG_SECTIONS)

        # Read the config
        config = ConfigParser(inline_comment_prefixes=("#",))
        config.read(args.ini)

        # Check if the ini file is valid
        sanity_check_ini_file_content(config, CONFIG_SECTIONS)

        defaults = cls()  # default fallback values from dataclass

        # Set the receptor and ligand file names from the command line arguments
        props = dict(
            sys_name=sys_name,
            rec_fname=args.rec,
            lig_fname=args.lig,
        )

        for section, fields in CONFIG_SECTIONS.items():
            section_data = config[section] if section in config else {}

            for field in fields:
                value = section_data.get(field.name, None)

                # Fall back to the dataclass default if the field is missing
                if value is None:
                    props[field.name] = getattr(defaults, field.name)
                    continue

                # Allow explicit None in the config file (e.g. lipid_ff = None)
                if value.strip().lower() == "none":
                    props[field.name] = None
                    continue

                # Convert the value to the appropriate type
                value = field.data_type(value)

                # Attach the OpenMM unit if present
                if field.unit is not None:
                    value *= field.unit

                props[field.name] = value

        return cls(**props)
    
    def to_ini(self, out_fname, membrane_flag=False):
        """
        Write the simulation setup to an ini file. This can be used to reload the same setup later, or to use it as a template for other setups.
        """
        props = asdict(self)
        config = ConfigParser()
        
        for key, value in props.items():
            # retrieve the SimSetupField object for the current key
            field = get_field(key, CONFIG_SECTIONS)
            if field is None:
                continue
            section = get_section_of_field(key, CONFIG_SECTIONS) # get the name of the section
            if section not in config:
                config[section] = {}
            # add a comment with the unit string if the field has a unit
            if field.unit is not None:
                value = value.value_in_unit(field.unit)
                unit_str = openmm_unit_to_string(field.unit)
                config[section][key] = f"{value}  # {unit_str}"
            else:
                config[section][key] = str(value)

        # NOTE: This is a temporary solution to handle the membrane_flag, I'm not fully content with this solution. 
        # Ideally, we should have a more elegant way to handle this, but for now, this will suffice.
        if membrane_flag:
            # If membrane_flag is True, include the membrane section
            sections_to_include = ["system", "membrane", "simulation", "forcefields"]
            config["system"]["membrane"] = "True"
            # NOTE: The forcefield section will soon change so this fragment will also need to be adapted.
            config["forcefields"]["lipid_ff"] = "amber14/lipid17.xml" 
        else:
            # If membrane_flag is False, exclude the membrane section
            sections_to_include = ["system", "simulation", "simulation box", "forcefields"]

        # Create a filtered config & filter sections
        filtered_config = ConfigParser()
        for section in sections_to_include:
            if section in config:
                filtered_config[section] = dict(config[section])

        with open(out_fname, "w") as f:
            filtered_config.write(f)
        return None


# NOTE: this function is not actually used anywhere in the code anymore
def get_longest_key(keys):
    longest_string = max([len(key) for key in keys])
    return longest_string+1


def sanity_check_ini_file_if_exists(ini_filename, field_sections):
    """
    Checks if the ini file exists. Raises a FileNotFoundError if the file does not exist.
    """
    if not os.path.isfile(ini_filename):
        raise FileNotFoundError(
            f"Provided config .ini file not found: {ini_filename}\n"
            "Use simprepper-example-config to generate a template configuration."
        )

def sanity_check_ini_file_content(config_content, field_sections):
    """
    Checks if the ini file is valid. Raises a ValueError if there are unknown 
    sections or fields in the ini file.
    """
    # Raise an error if there are unknown sections in the ini file, to help users identify typos or misplaced options
    expected_sections = set(field_sections.keys())
    actual_sections = set(config_content.sections())
    unknown_sections = actual_sections - expected_sections
    if unknown_sections:
        raise ValueError(
            f"Unknown section(s) in config file: {', '.join(sorted(unknown_sections))}\n"
            f"Expected sections: {', '.join(sorted(expected_sections))}"
        )
    
    # Check for unknown fields in each section
    actual_fields = {section: set(config_content[section].keys()) for section in actual_sections}
    unknown_fields = {}
    for section in actual_sections:
        expected_fields = {field.name.lower(): field.name for field in field_sections.get(section, [])}
        unknown = {field for field in actual_fields[section] if field.lower() not in expected_fields}
        if unknown:
            unknown_fields[section] = {
                "unknown": unknown,
                "expected": sorted(expected_fields.values())
            }
    # Raise an error
    if unknown_fields:
        msg = ["Unknown field(s) found in config file:"]
        for section, info in unknown_fields.items():
            msg.append(f"Found in [{section}]: {', '.join(sorted(info['unknown']))}")
            msg.append(f"  Expected: {', '.join(info['expected'])}")
        raise ValueError("\n".join(msg))


def get_section_of_field(field_name, field_sections):
    """
    Returns the section name of a given field.
    """
    for section, fields in field_sections.items():
        for field in fields:
            if field.name == field_name:
                return section
    return None  # if the field doesn't belong to any section


def get_field(field_name, field_sections):
    """
    Returns the SimSetupField object for a given field name.
    """
    for fields in field_sections.values():
        for field in fields:
            if field.name == field_name:
                return field
    return None

def get_field_default(field_name, field_sections):
    """
    Get the default value of a field and add OpenMM unit if present.
    """
    for fields in field_sections.values():
        for field in fields:
            if field.name == field_name:
                if field.unit is not None:
                    return field.default * field.unit
                return field.default
    raise KeyError(f"Field '{field_name}' not found.")