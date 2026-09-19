from openmm import unit as mm_units
from collections import namedtuple
from openff.units.openmm import openmm_unit_to_string
from configparser import ConfigParser
import argparse
import warnings
import os
from simprepper.utils.utils import quit_with_error
from simprepper.utils.config_utils import SystemSection, WaterBoxSection 
from simprepper.utils.config_utils import MembraneBoxSection, SimulationSection, ForceFieldsSectionAmber
from simprepper.utils.config_utils import sanity_check_ini_file_if_exists

os.environ["JAX_ENABLE_X64"] = "True" # get rid of annoying JAX warning

#TODO: separate Setup from Config
# Currently, most of the Config related stuff, is still in the constructor of the Setup Class...

# %% global constants

# NOTE: If you want to add a new parameter, you need to add it to one specific configsections.
# also NOTE: you cannot have the same name for parameters, even if they are in different sections!

# This would be the new version, using the Section class:
_CONFIG_SECTIONS_vanilla = [SystemSection.from_kwargs(membrane_protein=False),
                            WaterBoxSection, 
                            SimulationSection, 
                            ForceFieldsSectionAmber]
_CONFIG_SECTIONS_vanilla = {section.section_name:section for section in _CONFIG_SECTIONS_vanilla}
SECTION_LABELS_vanilla   = list(_CONFIG_SECTIONS_vanilla.keys())
#                        = ["system", "water_box",    "simulation", "forcefields"]


_CONFIG_SECTIONS_membrane = [SystemSection.from_kwargs(membrane_protein=True),
                             MembraneBoxSection,  # does not exist, yet
                             SimulationSection, 
                             ForceFieldsSectionAmber]
_CONFIG_SECTIONS_membrane = {section.section_name:section for section in _CONFIG_SECTIONS_membrane}
SECTION_LABELS_membrane   = list(_CONFIG_SECTIONS_membrane.keys())
#                         = ["system", "membrane_box", "simulation", "forcefields"]


def check_for_unknown_sections(configfile_content, expected_sections):
    expected_section_keys = set(expected_sections.keys())
    actual_section_keys = set(configfile_content.sections())
    unknown_sections = actual_section_keys - expected_section_keys
    if unknown_sections:
        raise ValueError(
            f"Unknown section(s) in config file: {', '.join(sorted(unknown_sections))}\n"
            f"Expected sections: {', '.join(sorted(expected_section_keys))}"
        )
    return actual_section_keys


def sanity_check_ini_file_content(configfile_content, expected_sections):
    """
    Checks if the ini file is valid. Raises a ValueError if there are unknown 
    sections or fields in the ini file.
    """
    # Raise an error if there are unknown sections in the ini file, to help users identify typos or misplaced options
    sections_in_file = check_for_unknown_sections(configfile_content, expected_sections)
    
    # Check for unknown fields in each section
    # TODO: reimplement this check!
    # actual_fields = {section: set(config_content[section].keys()) for section in actually_given_section_keys}
    # unknown_fields = {}
    # for section_key in actually_given_section_keys:
    #     section_of_interest = expected_sections.get(section_key)
    #     print([field for field in section_of_interest.default_fields])
    #     expected_fields = {field.name.lower(): field.name for field in section_of_interest.default_fields}
    #     unknown = {field for field in actual_fields[section_key] if field.lower() not in expected_fields.keys()}
    #     if unknown:
    #         unknown_fields[section_key] = {
    #             "unknown": unknown,
    #             "expected": sorted(expected_fields.values())
    #         }
    # # Raise an error
    # if unknown_fields:
    #     msg = ["Unknown field(s) found in config file:"]
    #     for section, info in unknown_fields.items():
    #         msg.append(f"Found in [{section}]: {', '.join(sorted(info['unknown']))}")
    #         msg.append(f"  Expected: {', '.join(info['expected'])}")
    #     raise ValueError("\n".join(msg))


def check_for_membrane_protein(config):
    """Do a config preflight, to check if we have a membrane protein, or not
    """
    # hack: if capitalization in config is wrong:
    # map: any capitalization fo the lowcaps version
    lowcaps_keys = {key.lower():key for key in config.keys()}
    
    if "system" in lowcaps_keys.keys():
        #TODO : potential bug "System"...
        section_of_interest = config[lowcaps_keys["system"]]
        if section_of_interest is None:
            return False
        lowcaps_kkeys = {key.lower():key for key in section_of_interest.keys()}
        if "membrane_protein" in lowcaps_kkeys.keys():
            field_of_interest = section_of_interest[lowcaps_kkeys["membrane_protein"]]
            return str2bool(field_of_interest)
    return False


# NOTE: Essentially the same function exists in argument_parsing as well, *BUT* it raises a different error!!
def str2bool(v) -> bool:
    """This function enables the translation of user input to actual boolean values
    otherwise bool("false") would return `True`,
    which is really not what the user wants...
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise TypeError('Boolean value expected.')


def parse_value(value, field):
    # None: Parsing of "none" would be disturbed by the type-cast below, if we do not directly return it...
    if value.strip().lower() == "none":
        return None

    # Convert the value to the appropriate type
    if field.data_type is bool:
        # boolean values are a special case (otherwise any string would be interpreted as True)
        value = str2bool(value)
    else:
        value = field.data_type(value)
    return value


class SimPrepperConfig:
    """This is a need little class to collect all the fields from all the sections.
    It has similarities to a @dataclass, but it is way less strict...
    """
    _fields_dict:   dict       = {}
    fields:         namedtuple = None

    def __init__(self, fields_dict):
        self._fields_dict       = fields_dict
        # this namedtuple-hack is really neat
        fields = namedtuple("ConfigFields", list(fields_dict.keys()))(**fields_dict)
        self.fields = fields
        return None

    def __getattr__(self, name):
        """This allows directly accessing fields as an attribute of `SimPrepperConfig` 
        """
        if name in self._fields_dict.keys():
            # by default return the value property, without a unit...
            return self._fields_dict[name].with_or_without_u
        print("Unknown setup attribute: {}".format(name))
        return None

    def has_field(self, field_key):
        return field_key in self._fields_dict.keys()

    def print_fields(self):
        """Mostly for debugging...
        """
        for key,val in self._fields_dict.items():
            print(key, val)
        return None



# %% DEVELOPMENT: This is supposed to substitute SimSetup (still defined below)
class SimSetup:
    """This class contains a file filenames, and the sys_name.
    All settings, are collected in the `SimPrepperConfig`, which is constructed from an ini file, 
    and/or default values.
    """
    sys_name        : str                       = "DEFAULT_SYS_NAME"
    rec_fname       : str                       = "DEFAULT_REC_NAME"
    lig_fname       : str | None                = None
    config_sections : dict                      = {}
    config:           SimPrepperConfig | None   = None  

    def __init__(self, sys_name, rec_fname, lig_fname, sections):
        self.sys_name           = sys_name
        self.rec_fname          = rec_fname
        self.lig_fname          = lig_fname
        self.config_sections    = sections

        # What happens below? We flatten all fields to a single dictionary and remove the hierarchical layer of the sections
        fields_dict ={}
        for section in sections.values():   # we do not actually care about the section names any more...
            for field_key,field in section.fields.items():
                fields_dict[field_key] = field

        self.config = SimPrepperConfig(fields_dict)
        return None

    #TODO: Eventually remove this deprecated function
    @classmethod
    def from_args(cls) -> None:
        """This was a constructor function, to generate an instance of SimSetup() from argparse-arguments
        """
        warnings.warn("SimSetup.from_args() is deprecated!", DeprecationWarning)
        quit_with_error("SimSetup.from_args() is deprecated!")
        return None

    @classmethod
    def from_defaults(cls, sys_name=None, rec_fname=None, lig_fname=None, membrane_flag=False):
        """
        This is a constructor function, to generate an instance of SimSetup() from the defaults defined above.
        """

        # Set the receptor and ligand file names from the command line arguments
        props = dict(
            sys_name=(sys_name or cls.sys_name),
            rec_fname=(rec_fname or cls.rec_fname),
            lig_fname=(lig_fname or cls.lig_fname),
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


    @classmethod
    def from_ini(cls, sys_name, ini_fname, rec_fname, lig_fname):
        """
        This is a constructor function, to generate an instance of SimSetup() from an ini file.
        """
        # Check if the ini file exists
        sanity_check_ini_file_if_exists(ini_fname)

        # Read the config
        config = ConfigParser(inline_comment_prefixes=("#",))
        config.read(ini_fname)

        membrane_flag = check_for_membrane_protein(config)

        sections=_CONFIG_SECTIONS_vanilla
        if membrane_flag:
            print("Detected a membrane protein config!!")
            sections=_CONFIG_SECTIONS_membrane
        # Check if the ini file is valid
        sanity_check_ini_file_content(config, sections)

        # Set the receptor and ligand file names from the command line arguments
        kwargs_for_init = dict(sys_name=sys_name,
                               rec_fname=rec_fname,
                               lig_fname=lig_fname,
                               )
        sections_out = {}
        for section_key, section in sections.items():
            section_loop = section.from_defaults()
            section_in_ini = config[section_key]
            for field_key, given_value in section_in_ini.items():
                parsed_value = parse_value(given_value, section_loop.fields[field_key])                
                section_loop.fields[field_key].set_value(parsed_value)

            sections_out[section_key] = section_loop
        kwargs_for_init["sections"] = sections_out

        return cls(**kwargs_for_init)

    def to_ini(self, out_fname, membrane_flag=False):
        """
        Write the simulation setup to an ini file. This can be used to reload the same setup later, or to use it as a template for other setups.
        """
        config = ConfigParser()
        for section_key, section in self.config_sections.items():
            # retrieve the SimSetupField object for the current key
            config[section_key] = {}
            for field_key, field in section.fields.items():
                value = str(field.value)
                # print(field_key, value)  # debug
                if field.unit is not None:  # add a comment with the unit string if the field has a unit
                    unit_str = openmm_unit_to_string(field.unit)
                    value = f"{value}  # {unit_str}"
                config[section_key][field_key] = f"{value} "

        # TODO: Check if this is even necessary anymore. I think if it is properly initialized, this does not need to be filrered anymore..
        # NOTE: This is a temporary solution to handle the membrane_flag, I'm not fully content with this solution. 
        # Ideally, we should have a more elegant way to handle this, but for now, this will suffice.
        filtered_config = get_filtered_config(config, membrane_flag)

        with open(out_fname, "w") as f:
            filtered_config.write(f)
        return None


def get_filtered_config(config, membrane_flag):
    sections_to_include = SECTION_LABELS_vanilla
    if membrane_flag:
        # If membrane_flag is True, include the membrane section
        sections_to_include = SECTION_LABELS_membrane
        config["system"]["membrane_protein"] = "True"
        # NOTE: The forcefield section will soon change so this fragment will also need to be adapted.
        config["forcefields"]["lipid_ff"] = "amber14/lipid17.xml" 
        
    # Create a filtered config & filter sections
    filtered_config = ConfigParser()
    for section in sections_to_include:
        if section in config:
            filtered_config[section] = dict(config[section])
    return filtered_config


# <---- DEVELOP

# this is the original version (until AUG/2026)
# CONFIG_SECTIONS_legacy = {
#     "system": [
#         SimSetupField("membrane_protein", None, False, bool),
#     ],
#     "membrane_box": [     # membrane parameters present only if membrane=True
#         SimSetupField("lipid_type", None, "POPC", str),
#         SimSetupField("membrane_center_z", None, 0.0, float),
#         SimSetupField("minimum_padding", mm_units.nanometer, 1.0, float),
#     ],
#     "simulation": [
#         SimSetupField("nb_cutoff", mm_units.nanometer, 1.0, float),
#         SimSetupField("hydrogen_mass", mm_units.amu, 4.0, float),
#         SimSetupField("timestep", mm_units.picoseconds, 0.004, float),
#         SimSetupField("temperature", mm_units.kelvin, 300.0, float),
#         SimSetupField("ionic_strength", mm_units.molar, 0.15, float),
#         SimSetupField("ph", None, 7.4, float),
#     ],
#     "water_box": [
#         SimSetupField("box_shape", None, "cube", str),
#         SimSetupField("padding", mm_units.nanometer, 3.0, float),
#     ],
#     "forcefields": [
#         SimSetupField("ligand_ff", None, "GAFF", str),
#         SimSetupField("protein_ff", None, "amber14/protein.ff14SB.xml", str),
#         SimSetupField("water_ff", None, "amber14/tip3pfb.xml", str),
#         SimSetupField("ion_ff", None, "amber/tip3p_HFE_multivalent.xml", str),
#         SimSetupField("lipid_ff", None, None, str),
#     ]
# }


# def get_section_of_field(field_name, field_sections):
#     """
#     Returns the section name of a given field.
#     """
#     for section, fields in field_sections.items():
#         for field in fields:
#             if field.name == field_name:
#                 return section
#     #QUESTION: Paq: Should the following ever happen?
#     return None  # if the field doesn't belong to any section


# def get_field(field_name, field_sections):
#     """
#     Returns the SimSetupField object for a given field name.
#     """
#     for fields in field_sections.values():
#         for field in fields:
#             if field.name == field_name:
#                 return field
#     return None


# def get_field_default(field_name, field_sections):
#     """
#     Get the default value of a field and add OpenMM unit if present.
#     """
#     for fields in field_sections.values():
#         for field in fields:
#             if field.name == field_name:
#                 if field.unit is not None:
#                     return field.default_value * field.unit
#                 return field.default_value
#     raise KeyError(f"Field '{field_name}' not found.")


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

    SimSetup.from_defaults(membrane_flag=args.membrane).to_ini(args.output, args.membrane)

    print(f"Wrote example config to {args.output}")
