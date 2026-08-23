
import os
from collections import namedtuple
from dataclasses import dataclass, field
from openmm import unit as mm_units


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


# %% Sections for the config file
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
class MembraneBoxSection(AbstractConfigSection):
    section_name       = "membrane_box"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("lipid_type", None, "POPC", None, str),
                          SimSetupField("membrane_center_z", None, 0.0, None, float),
                          SimSetupField("minimum_padding", mm_units.nanometer, 1.0, None, float)
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



def sanity_check_ini_file_if_exists(ini_filename):
    """
    Checks if the ini file exists. Raises a FileNotFoundError if the file does not exist.
    """
    if not os.path.isfile(ini_filename):
        raise FileNotFoundError(
            f"Provided config .ini file not found: {ini_filename}\n"
            "Use simprepper-example-config to generate a template configuration."
        )

