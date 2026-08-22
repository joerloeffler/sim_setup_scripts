from openmm import unit as mm_units
from openmm.unit.quantity import Quantity as mm_quantity
from openff.units.openmm import openmm_unit_to_string
from dataclasses import dataclass, field, asdict
from collections import namedtuple
from configparser import ConfigParser
import argparse
# import warnings  # not used here
import os

os.environ["JAX_ENABLE_X64"] = "True" # get rid of annoying JAX warning

# %% global constants
# Sections for the config file
# NOTE: At the moment, when adding a new parameter, you need to add it to the CONFIG_SECTIONS dictionary, 
# and also add it to the SimSetup dataclass. Ideally, we would like to have a single source of truth for 
# the parameters, but this can be solved later.
SimSetupField = namedtuple("SimSetupField", ["name", 
                                             "unit", 
                                             "default_value", # this will always be passed on
                                             "given_value", # this will be None, if not initialized differently
                                             "data_type"])


def initialize_field_from_kwarg(default_field_behind_key, val):
    if type(val) != default_field_behind_key.data_type:
        #TODO: actually handle this type comparison...
        print("Oh Boy! This should have a different type!")
        
    # initialize the field properly, with units and all that, not just the value:
    field_from_user_input = SimSetupField(default_field_behind_key.name,
                                          default_field_behind_key.unit,
                                          default_field_behind_key.default_value, 
                                          val,
                                          default_field_behind_key.data_type)
    return field_from_user_input


# %% Template class for Config-sections
@dataclass
class AbstractConfigSection:
    section_name         : str   = ""
    default_fields       : dict  = field(default_factory=dict)
    fields               : dict  = field(default_factory=dict)

    def __init__(self, input_fields):
        self.fields = input_fields

    @classmethod
    def from_defaults(cls):
        """Constructs an instance using the class's default_fields dict."""
        initial_fields = {
            name: setup_field  for name, setup_field in cls.default_fields.items()
            }
        
        return cls(section_name=cls.section_name,
                   fields=initial_fields)

    @classmethod
    def from_kwargs(cls, **kwargs):
        """Constructs an instance using input arguments. 
        For any field that is missing, use default instead..."""
        initial_fields = {
                    name: setup_field  for name, setup_field in cls.default_fields.items()
                    }
        # check the kwargs, if they exists in the default definition. 
        # Otherwise arbitrary properties could be added
        filtered_kwargs = {}
        for key,val in kwargs.items():
            if key not in cls.default_fields.keys():
                print("Ingnoring unkown field {}".format(key))  # This should never happen, from ini-files,
                # because in ini files, we already sanity check the fields before initialization.
                continue
            default_field_behind_key = cls.default_fields[key]
            filtered_kwargs[key] = initialize_field_from_kwarg(default_field_behind_key,val)
        # Update default values with any explicit kwargs that have been given
        initial_fields.update(filtered_kwargs)
    
        return cls(section_name=cls.section_name,
                   fields=initial_fields)

@dataclass
class SystemSection(AbstractConfigSection):
    section_name       = "system"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("membrane_protein", None, False, None, bool)
                          ]
    default_fields     = {f.name: f for f in _list_of_defaults}


@dataclass
class WaterBoxSection(AbstractConfigSection):
    section_name       = "water_box"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("box_shape", None, "cube", None, str),
                          SimSetupField("padding", mm_units.nanometer, 3.0, None, float),
                          ]
    default_fields     = {f.name: f for f in _list_of_defaults}

@dataclass
class SimulationSection(AbstractConfigSection):
    section_name       = "simulation"
    # class-level definition of the default values
    _list_of_defaults  = [
        SimSetupField("nb_cutoff", mm_units.nanometer, 1.0, None, float),
        SimSetupField("hydrogen_mass", mm_units.amu, 4.0, None, float),
        SimSetupField("timestep", mm_units.picoseconds, 0.004, None, float),
        SimSetupField("temperature", mm_units.kelvin, 300.0, None, float),
        SimSetupField("ionic_strength", mm_units.molar, 0.15, None, float),
        SimSetupField("ph", None, 7.4, None, float),
    ]
    default_fields    = {f.name: f for f in _list_of_defaults}


@dataclass
class ForceFieldsSectionAmber(AbstractConfigSection):
    section_name       = "forcefields"
    # class-level definition of the default values
    _list_of_defaults  = [
        SimSetupField("ligand_ff", None, "GAFF", None, str),
        SimSetupField("protein_ff", None, "amber14/protein.ff14SB.xml", None, str),
        SimSetupField("water_ff", None, "amber14/tip3pfb.xml", None, str),
        SimSetupField("ion_ff", None, "amber/tip3p_HFE_multivalent.xml", None, str),
        SimSetupField("lipid_ff", None, "amber14/lipid17.xml", None, str),
    ]
    default_fields     = {f.name: f for f in _list_of_defaults}


# This would be the new version, using the Section class:
_CONFIG_SECTIONS_vanilla = [SystemSection.from_kwargs(membrane_protein=False),
                            WaterBoxSection, 
                            SimulationSection, 
                            ForceFieldsSectionAmber]
_CONFIG_SECTIONS_vanilla = {section.section_name:section for section in _CONFIG_SECTIONS_vanilla}
#TODO: implement _CONFIG_SECTIONS_membrane!
_CONFIG_SECTIONS_membrane = [SystemSection.from_kwargs(membrane_protein=False),
                            #MembraneBoxSection,  # does not exist, yet
                            SimulationSection, 
                            ForceFieldsSectionAmber]
_CONFIG_SECTIONS_membrane = {section.section_name:section for section in _CONFIG_SECTIONS_membrane}


SECTION_LABELS_vanilla     = list(_CONFIG_SECTIONS_vanilla.keys())
#SECTION_LABELS_vanilla     = ["system", "water_box",    "simulation", "forcefields"]

#MEMBRANE_PROTEIN_SECTIONS  = ["system", "membrane_box", "simulation", "forcefields"]
SECTION_LABELS_membrane  = ["system", "membrane_box", "simulation", "forcefields"]
# SECTION_LABELS_membrane     = list(_CONFIG_SECTIONS_membrane.keys())


# %% DEVELOPMENT: This is supposed to substitute SimSetup (still defined below)
@dataclass
class SimSetup_fromAbstractSections:
    """For developers: 
    if you want to add a property to this class:
    1. Add it here, to be a field that can be initialized
    2. make sure that it is passed on, in the constructor class `.from_args()`

    Note: you can define default values here, but they can be overwritten during construction (by other default values)
    So defining default values is not actually necessary. 
    HOWEVER, if you do not define a default value, it will be considered a non-default argument, which cannot follow a default argument...
    """
    sys_name        : str         = "DEFAULT_SYS_NAME"
    rec_fname       : str         = "DEFAULT_REC_NAME"
    lig_fname       : str | None  = None
    sections        : dict        = field(default_factory=dict)
    fields          : dict        = field(default_factory=dict)

    @classmethod
    def from_defaults(cls, membrane_flag=False):
        """
        This is a constructor function, to generate an instance of SimSetup() from the defaults defined above.
        """

        # Set the receptor and ligand file names from the command line arguments
        props = dict(
            sys_name=cls.sys_name,
            rec_fname=cls.rec_fname,
            lig_fname=cls.lig_fname,
        )
        
        sections=_CONFIG_SECTIONS_vanilla
        if membrane_flag:
            sections=_CONFIG_SECTIONS_membrane

        # What happens below? We initialize the sections from the default values as defined in the section classes above
        sections_out = {}
        for section_key, section in sections.items():
            default_section = section.from_defaults()
            sections_out[section_key] = default_section
        return cls(**props, 
                   sections=sections_out)  # call __init__ 

    def __init__(self, sys_name, rec_fname, lig_fname, sections):
        self.sys_name     = sys_name
        self.rec_fname    = rec_fname
        self.lig_fname    = lig_fname
        self.sections     = sections

        # What happens below? We flatten all fields to a single dictionary and remove the hierarchical layer of the sections
        fields ={}
        # for section_key,section in sections.items():  # we do not actually care about the keys here...
        for section in sections.values():
            for field_key,field in section.fields.items():
                fields[field_key] = field

        self.fields       = fields
        return None

    @classmethod
    def from_ini(cls, sys_name, args):
        """
        This is a constructor function, to generate an instance of SimSetup() from an ini file.
        """
        #TODO: implement read values to given, values, otherwise use default values
        pass

    def to_ini(self, out_fname, membrane_flag=False):
        """
        Write the simulation setup to an ini file. This can be used to reload the same setup later, or to use it as a template for other setups.
        """
        #TODO: implement, use given_values if not None, otherwise default values
        pass


# <---- DEVELOP

# this is the original version (until AUG/2026)
CONFIG_SECTIONS = {
    "system": [
        SimSetupField("membrane_protein", None, False, None, bool),
    ],
    "membrane_box": [     # membrane parameters present only if membrane=True
        SimSetupField("lipid_type", None, "POPC", None, str),
        SimSetupField("membrane_center_z", None, 0.0, None, float),
        SimSetupField("minimum_padding", mm_units.nanometer, 1.0, None, float),
    ],
    "simulation": [
        SimSetupField("nb_cutoff", mm_units.nanometer, 1.0, None, float),
        SimSetupField("hydrogen_mass", mm_units.amu, 4.0, None, float),
        SimSetupField("timestep", mm_units.picoseconds, 0.004, None, float),
        SimSetupField("temperature", mm_units.kelvin, 300.0, None, float),
        SimSetupField("ionic_strength", mm_units.molar, 0.15, None, float),
        SimSetupField("ph", None, 7.4, None, float),
    ],
    "water_box": [
        SimSetupField("box_shape", None, "cube", None, str),
        SimSetupField("padding", mm_units.nanometer, 3.0, None, float),
    ],
    "forcefields": [
        SimSetupField("ligand_ff", None, "GAFF", None, str),
        SimSetupField("protein_ff", None, "amber14/protein.ff14SB.xml", None, str),
        SimSetupField("water_ff", None, "amber14/tip3pfb.xml", None, str),
        SimSetupField("ion_ff", None, "amber/tip3p_HFE_multivalent.xml", None, str),
        SimSetupField("lipid_ff", None, None, None, str),
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
    membrane_protein : bool        = field(default_factory=lambda: get_field_default("membrane_protein", CONFIG_SECTIONS))
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
        sanity_check_ini_file_if_exists(args.ini)

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
                if field.data_type is bool:
                    # boolean values are a special case (otherwise any string would be interpreted as True)
                    value = section_data.getboolean(field.name)
                else:
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
            sections_to_include = ["system", "membrane_box", "simulation", "forcefields"]
            config["system"]["membrane_protein"] = "True"
            # NOTE: The forcefield section will soon change so this fragment will also need to be adapted.
            config["forcefields"]["lipid_ff"] = "amber14/lipid17.xml" 
        else:
            # If membrane_flag is False, exclude the membrane section
            sections_to_include = ["system", "simulation", "water_box", "forcefields"]

        # Create a filtered config & filter sections
        filtered_config = ConfigParser()
        for section in sections_to_include:
            if section in config:
                filtered_config[section] = dict(config[section])

        with open(out_fname, "w") as f:
            filtered_config.write(f)
        return None


# NOTE: this function is not actually used anywhere in the code anymore
#def get_longest_key(keys):
#    longest_string = max([len(key) for key in keys])
#    return longest_string+1


def sanity_check_ini_file_if_exists(ini_filename):
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
    #QUESTION: Paq: Should the following ever happen?
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
                    return field.default_value * field.unit
                return field.default_value
    raise KeyError(f"Field '{field_name}' not found.")
