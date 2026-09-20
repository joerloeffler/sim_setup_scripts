
import os
# from collections import namedtuple
# from dataclasses import dataclass, field
from openmm import unit as mm_units

#TODO: remove this eventually
# SimSetupField = namedtuple("SimSetupField", ["name", 
#                                              "unit", 
#                                              "default_value", # this will always be passed on
#                                              "given_value", # this will be None, if not initialized differently
#                                              "data_type"])

class SimSetupField:
    name:           str | None  = None
    unit                        = None    
    default_value               = None  # this will always be passed on
    given_value                 = None  # this will be None, if not initialized differently
    data_type                   = None
    value_was_set:  bool        = False

    def __init__(self, name, unit, default_value, data_type):
        self.name = name 
        self.unit = unit 
        self.default_value = default_value 
        self.data_type = data_type 
        # Note: Given value is not set via init!
        return None
        
    def set_value(self, value):
        self.given_value = value
        self.value_was_set = True
        return self  # this is very convienient, because you can initialize and set_value in one line (alled mutator)

    @property
    def value(self):
        """Only give the value, regardless of whether this field has a unit or not...
        """
        if self.value_was_set:
            return self.given_value
        return self.default_value

    @property
    def value_with_unit(self):
        """This property, will probably never be used.
        """
        if self.unit is None:
            print("SetupField {} does not have a unit...".format(self.name))
            return self.value
        return self.value * self.unit

    def value_in_unit(self, unit):
        """Construct an openmm `quantity`, and immediately use the builtin function
        """
        if self.unit is None:  # in this case, we cannot construct an openmm `quantity``
            return self.value
        return self.value_with_unit.value_in_unit(unit)

    @property
    def with_or_without_u(self):
        """Great song by U2!
        ... but in this case it means _unit
        -> give the value of the fields, if this field has a unit then add the unit, too.
        """
        if self.unit is None:
            return self.value
        return self.value_with_unit

    def __repr__(self):
        all_vars_except_name = {key:val for (key,val) in vars(self).items() if not key=="name"}
        return "{}({})".format(self.name, all_vars_except_name)


def initialize_field_from_kwarg(default_field_behind_key, val):
    if type(val) != default_field_behind_key.data_type:
        #TODO: actually handle this type comparison...
        print("Oh Boy! This should have a different type!")
        
    # initialize the field properly, with units and all that, not just the value:
    field_from_user_input = SimSetupField(default_field_behind_key.name,
                                          default_field_behind_key.unit,
                                          default_field_behind_key.default_value, 
                                          default_field_behind_key.data_type)
    field_from_user_input.set_value(val)
    return field_from_user_input


# %% Template class for Config-sections
class AbstractConfigSection:
    section_name         : str   = ""
    default_fields       : dict
    fields               : dict

    def __init__(self, input_fields, section_name=None):
        self.section_name = section_name or self.section_name
        self.fields = input_fields
        return None

    @classmethod
    def from_defaults(cls):
        """Constructs an instance using the class's default_fields dict."""
        initial_fields = {
            name: setup_field  for name, setup_field in cls.default_fields.items()
            }
        
        return cls(section_name=cls.section_name,
                   input_fields=initial_fields)

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
                   input_fields=initial_fields)


# %% Sections for the config file
class SystemSection(AbstractConfigSection):
    section_name       = "system"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("membrane_protein", None, False, bool)
                          ]
    default_fields     = {f.name: f for f in _list_of_defaults}


class WaterBoxSection(AbstractConfigSection):
    section_name       = "water_box"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("box_shape", None, "cube", str),
                          SimSetupField("padding", mm_units.nanometer, 3.0, float),
                          ]
    default_fields     = {f.name: f for f in _list_of_defaults}


class MembraneBoxSection(AbstractConfigSection):
    section_name       = "membrane_box"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("lipid_type", None, "POPC", str),
                          SimSetupField("membrane_center_z", mm_units.nanometer, 0.0, float),
                          SimSetupField("minimum_padding", mm_units.nanometer, 1.0, float)
                          ]
    default_fields     = {f.name: f for f in _list_of_defaults}


class SimulationSection(AbstractConfigSection):
    section_name       = "simulation"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("nb_cutoff", mm_units.nanometer, 1.0, float),
                          SimSetupField("hydrogen_mass", mm_units.amu, 4.0, float),
                          SimSetupField("timestep", mm_units.picoseconds, 0.004, float),
                          SimSetupField("temperature", mm_units.kelvin, 300.0, float),
                          SimSetupField("ionic_strength", mm_units.molar, 0.15, float),
                          SimSetupField("ph", None, 7.4, float),
                          ]
    default_fields    = {f.name: f for f in _list_of_defaults}


class ForceFieldsSectionAmber(AbstractConfigSection):
    section_name       = "forcefields"
    # class-level definition of the default values
    _list_of_defaults  = [SimSetupField("ligand_ff", None, "GAFF", str),
                          SimSetupField("protein_ff", None, "amber14/protein.ff14SB.xml", str),
                          SimSetupField("water_ff", None, "amber14/tip3pfb.xml", str),
                          SimSetupField("ion_ff", None, "amber/tip3p_HFE_multivalent.xml", str),
                          SimSetupField("lipid_ff", None, "amber14/lipid17.xml", str),
                          ]
    default_fields     = {f.name: f for f in _list_of_defaults}


def check_if_ini_file_exists(ini_filename):
    """
    Checks if the ini file exists. Raises a FileNotFoundError if the file does not exist.
    """
    if not os.path.isfile(ini_filename):
        raise FileNotFoundError(
            f"Provided config .ini file not found: {ini_filename}\n"
            "Use simprepper-example-config to generate a template configuration file."
        )

