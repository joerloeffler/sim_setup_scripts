from collections import namedtuple
from simprepper.utils.config_utils import SystemSection, WaterBoxSection 
from simprepper.utils.config_utils import MembraneBoxSection, SimulationSection, ForceFieldsSectionAmber

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


def check_for_unknown_sections(configfile_content, expected_sections):
    expected_section_keys = set(expected_sections.keys())
    given_section_keys = set(configfile_content.sections())
    unknown_sections = given_section_keys - expected_section_keys
    if unknown_sections:
        raise ValueError(
            f"Unknown section(s) in config file: {', '.join(sorted(unknown_sections))}\n"
            f"Expected sections: {', '.join(sorted(expected_section_keys))}"
        )
    return given_section_keys


def sanity_check_ini_file_content(configfile_content, expected_sections):
    """
    Checks if the ini file is valid. Raises a ValueError if there are unknown 
    sections or fields in the ini file.
    """
    # Raise an error if there are unknown sections in the ini file, to help users identify typos or misplaced options
    given_section_keys_in_file = check_for_unknown_sections(configfile_content, expected_sections)
    
    # Check for unknown fields in each section
    given_fields_in_sections = {section: set(configfile_content[section].keys()) for section in given_section_keys_in_file}
    unknown_fields = {}
    for section_key in given_section_keys_in_file:
        section_of_interest = expected_sections.get(section_key)

        # print(section_key, [field for field in section_of_interest.default_fields])  # debug
        expected_fields = list(section_of_interest.default_fields.keys())
        unknown = [field for field in given_fields_in_sections[section_key] if field.lower() not in expected_fields]
        if unknown:
            unknown_fields[section_key] = {"unknown": unknown,
                                           "expected": sorted(expected_fields)
                                           }
    # Raise an error
    if unknown_fields:
        msg = ["\nUnknown field(s) found in config file!"]
        for section, info in unknown_fields.items():
            msg.append(f"Fields found in section [{section}]: {', '.join(sorted(info['unknown']))}")
            msg.append(f"  Expected: {', '.join(info['expected'])}")
        raise ValueError("\n".join(msg))
    return None


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


