import logging
from openmm import unit as mm_units

from simprepper.sim_setup import SimSetup

# Maybe this function should be somewhere else...
def quit_with_error(msg="Aborting, because of erroneous setup for preparation. Please fine the report above."):
    logging.error(msg)
    raise RuntimeError
    

# NOTE: This class is the important feature here
class CheckItem:
    """ The expression of the test condition should return a boolean.
    Best expressed in terms of a lambda function, which accepts a SimSetup-object.
    Theoretically, you could also pass on an actual function as a condition. 
    Maybe this might be done for some particularly complex logical conditions...
    """
    name        : str  = ""
    explanation : str  = ""
    should_quit : bool = False
    condition          = None
    skip_condition     = None
    has_failed         = None
    was_skipped : bool = False
    result      : str  = ""
    def __init__(self, name, condition, explanation, should_quit, skip_condition=None):
        self.name = name
        self.explanation = explanation
        self.should_quit = should_quit
        self.condition = condition
        self.skip_condition = skip_condition
        return None

    def perform_check(self, sim_setup):
        if self.skip_condition is not None:
            if self.skip_condition(sim_setup):
                self.was_skipped = True
                self.result="Skipped"
                return None
        test_result = self.condition(sim_setup)
        self.has_failed = test_result
        result = {True:"Failed", False:"Passed"}.get(test_result)
        self.result=result
        return result


# Example for devs:
#TODO: Currently the SimSetup class does not include sections, but only the fields from the sections, so we cannot actually check for sections...
template_item = CheckItem(name="Mutually exclusive sections: [propA] and [propB]",
                          condition=lambda sim_setup: 
                          (sim_setup["propA"]==True and 
                           sim_setup["propB"]==True),
                          explanation="Dear User, propA and propB are mutually exclusive...",
                          should_quit=True  # if this will definitely lead to an erroneous setup, instruct termination of simprepper
                          )
#NOTE: Some checks will not be performed, if some other condition
template_item = CheckItem(name="sign of [propA], but propA is only given if propB is false.",
                          condition=lambda sim_setup: 
                          (sim_setup["propA"]<0),
                          explanation="Dear User, propA should be larger than 0...",
                          skip_condition=lambda sim_setup: sim_setup["propB"]==True,
                          should_quit=False  # if the anticipated error happens after our RunTime, we let it slide
                          )

# Actual list of checks for runtime
LIST_OF_CHECKS = [CheckItem(name="Mutually exclusive fields: membrane_protein and box_shape",
                            condition=lambda sim_setup: 
                            (sim_setup.membrane_protein==True and 
                             sim_setup.box_shape is not None),
                            explanation="If you prepare a membrane protein, you cannot use the [water box] section in the config...",
                            should_quit=False  # TODO: This test currently fails with default configs, because if nothing was given, there is always a default for box_shape...
                            ),
                  CheckItem(name="Sign of ionic strength",
                            condition=lambda sim_setup: 
                            (sim_setup.ionic_strength._value < 0.0),  # note we don't use .value_in_unit, because for the sign, the unit doesnt matter...
                            explanation="Ionic strength shouldn't be negative",
                            should_quit=True  # modelling of solvent box will fail
                            ),
                  CheckItem(name="Sign of timestep",
                            condition=lambda sim_setup: 
                            (sim_setup.timestep._value < 0.0),  # note we don't use .value_in_unit, because for the sign, the unit doesnt matter...
                             explanation="Timestep shouldn't be negative",
                             should_quit=False  # This error will only happen after our RunTime, so we let it slide...
                             ),
                  CheckItem(name="Sign of temperature",
                            condition=lambda sim_setup: 
                            (sim_setup.temperature.value_in_unit(mm_units.kelvin) < 0.0),
                            explanation="Temperature (in Kelvin) shouldn't be negative",
                            should_quit=False  # This error will only happen after our RunTime, so we let it slide...
                            ),
                  CheckItem(name="Sign of box_padding",
                            condition=lambda sim_setup: 
                            (sim_setup.padding._value < 0.0),
                            explanation="Padding shouldn't be negative",
                            skip_condition=lambda sim_setup: sim_setup.membrane_protein==True,
                            should_quit=False  # This error will only happen after our RunTime, so we let it slide...
                            ),  
                  CheckItem(name="Sign of minimum padding",
                            condition=lambda sim_setup: 
                            (sim_setup.minimum_padding._value < 0.0),
                            explanation="Minimum padding shouldn't be negative",
                            skip_condition=lambda sim_setup: sim_setup.membrane_protein==False,
                            should_quit=False  # This error will only happen after our RunTime, so we let it slide...
                            ),  
                  # ...
                  ]


# This is the workhorse here, but in the best base, this never needs to be touched
class SetupChecker:
    sim_setup          : SimSetup | None = None
    list_of_checks     : list            = []
    verbose            : bool            = False
    _results_of_checks : list            = []
    _has_finished      : bool            = False
    logging = None

    def __init__(self, 
                 sim_setup, 
                 logging=None,
                 list_of_checks=LIST_OF_CHECKS,
                 verbose=False):
        self.sim_setup = sim_setup
        self.list_of_checks = list_of_checks
        self.verbose = verbose
        self.logging = logging
        return None

    def run_all_checks(self):
        if self.verbose:
            self.logging.info("Checking Setup for problems...")
        for check in self.list_of_checks:
            #DEBUG:
            # self.logging.info(f"Running check: {check.name}")
            check.perform_check(self.sim_setup)
        self._has_finished = True
        return None

    def get_report(self):
        if not self._has_finished:
            self.run_all_checks()
        if self.verbose:  # print all checks, independent of results
            for i,check in enumerate(self.list_of_checks):    
                N_checks = len(self.list_of_checks)
                self.logging.info(f"SetupChecker sanity check {i+1:>2}/{N_checks} - >{check.name}< has: {check.result}")

        for check in self.list_of_checks:
            if check.has_failed:
                self.logging.warning(f"SetupCheck: >{check.name}< has failed!")
                self.logging.warning(check.explanation)
                if check.should_quit:
                    quit_with_error()
        return None
