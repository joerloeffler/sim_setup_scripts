from openff.units.openmm import openmm_unit_to_string
from configparser import ConfigParser
import argparse
import warnings
import os
from simprepper.utils.utils import quit_with_error
from simprepper.simprepper_config import SimPrepperConfig, _CONFIG_SECTIONS_membrane, _CONFIG_SECTIONS_vanilla,SECTION_LABELS_vanilla, SECTION_LABELS_membrane
from simprepper.simprepper_config import check_for_membrane_protein, sanity_check_ini_file_content, parse_value
from simprepper.utils.config_utils import check_if_ini_file_exists


os.environ["JAX_ENABLE_X64"] = "True" # get rid of annoying JAX warning

#TODO: separate Setup from Config
# Currently, most of the Config related stuff, is still in the constructor of the Setup Class...

def sections_from_config_file(template_sections, config_file_content):
    sections_out = {}
    for section_key, section in template_sections.items():
        section_loop = section.from_defaults()
        section_in_ini = config_file_content[section_key]
        for field_key, given_value in section_in_ini.items():
            parsed_value = parse_value(given_value, section_loop.fields[field_key])                
            section_loop.fields[field_key].set_value(parsed_value)

        sections_out[section_key] = section_loop
    return sections_out


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
        kwargs_for_init = dict(sys_name=(sys_name or cls.sys_name),
                               rec_fname=(rec_fname or cls.rec_fname),
                               lig_fname=(lig_fname or cls.lig_fname),
                               )

        #TODOL this should be SimPrepperConfig.from_defaults()
        template_sections=_CONFIG_SECTIONS_vanilla
        if membrane_flag:
            template_sections=_CONFIG_SECTIONS_membrane

        # What happens below? We initialize the sections from the default values as defined in the section classes above
        sections_out = {}
        for section_key, section in template_sections.items():
            default_section = section.from_defaults()
            sections_out[section_key] = default_section

        return cls(**kwargs_for_init, 
                   sections=sections_out)  # call __init__ 


    @classmethod
    def from_ini(cls, sys_name, ini_fname, rec_fname, lig_fname):
        """
        This is a constructor function, to generate an instance of SimSetup() from an ini file.
        """
        # Set the receptor and ligand file names from the command line arguments
        kwargs_for_init = dict(sys_name=sys_name,
                                rec_fname=rec_fname,
                                lig_fname=lig_fname,
                                )

        # Check if the ini file exists
        check_if_ini_file_exists(ini_fname)

        # Read the config
        config_file = ConfigParser(inline_comment_prefixes=("#",))
        config_file.read(ini_fname)

        membrane_flag = check_for_membrane_protein(config_file)

        template_sections = _CONFIG_SECTIONS_vanilla
        if membrane_flag:
            #TODO: use logging, not just print
            print("Detected a membrane protein config!!")
            template_sections = _CONFIG_SECTIONS_membrane
        # Check if the ini file is valid
        sanity_check_ini_file_content(config_file, template_sections)
        
        kwargs_for_init["sections"] = sections_from_config_file(template_sections=template_sections, 
                                                                config_file_content=config_file)

        return cls(**kwargs_for_init)

    #TODO: this should SimprepperConfig.to_ini()
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
