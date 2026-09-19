from openmm import unit as mm_units

from simprepper.sim_setup import SimSetup
from simprepper.utils.utils import quit_with_error    

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
#TODO: Currently the SimSetup class does not include sections, but only the fields from the sections, 
#   -> so we cannot actually check for sections...
template_item = CheckItem(name="Mutually exclusive sections: [propA] and [propB]",
                          condition=lambda config: 
                          (config.propA == True and 
                           config.propB ==True ),
                          explanation="Dear User, propA and propB are mutually exclusive...",
                          should_quit=True  # if this will definitely lead to an erroneous setup, instruct termination of simprepper
                          )
#NOTE: Some checks can be skipped, if some other condition is met (skip_condition)
template_item = CheckItem(name="sign of [propA]",
                          condition=lambda config: 
                          (config.propA < 0),
                          explanation="Dear User, propA should be larger than 0...",
                          skip_condition=lambda config: config.propB == True,
                          should_quit=False  # if the anticipated error happens after our RunTime, we let it slide
                          )

# Actual list of checks for runtime
LIST_OF_CHECKS = [CheckItem(name="Mutually exclusive fields: membrane_protein and box_shape",
                            condition=lambda config: 
                            (config.membrane_protein==True and 
                             config.has_field("box_shape") ),
                            explanation="If you prepare a membrane protein, you cannot use the [water box] section in the config...",
                            should_quit=False  # TODO: I think this test is actually obsolete, \
                            #                       -> because this is already captured during sanity checking of the config file...
                            ),
                
                  CheckItem(name="Sign of ionic strength",
                            condition=lambda config: 
                            (config.fields.ionic_strength.value < 0.0),  # note we ignore the unit, because for the sign, it doesnt matter...
                            explanation="Ionic strength shouldn't be negative",
                            should_quit=True  # modelling of solvent box will fail
                            ),
                  CheckItem(name="Sign of timestep",
                            condition=lambda config: 
                            (config.fields.timestep.value < 0.0),  # note we don't use .value_in_unit, because for the sign, the unit doesnt matter...
                             explanation="Timestep shouldn't be negative",
                             should_quit=False  # This will only lead to errors after the Runtime of `SimPrepper`, so we let it slide...
                             ),
                  CheckItem(name="Sign of temperature",
                            condition=lambda config: 
                            (config.fields.temperature.value_in_unit(mm_units.kelvin) < 0.0),
                            explanation="Temperature (in Kelvin) shouldn't be negative",
                            should_quit=False  # This will only lead to errors after the Runtime of `SimPrepper`, so we let it slide...
                            ),
                  CheckItem(name="Sign of box_padding",
                            condition=lambda config: 
                            (config.fields.padding.value < 0.0),
                            explanation="Padding shouldn't be negative",
                            skip_condition=lambda config: config.membrane_protein==True,
                            should_quit=False  # This will only lead to errors after the Runtime of `SimPrepper`, so we let it slide...
                            ),  
                  CheckItem(name="Sign of minimum padding",
                            condition=lambda config: 
                            (config.fields.minimum_padding.value < 0.0),
                            explanation="Minimum padding shouldn't be negative",
                            skip_condition=lambda config: config.membrane_protein==False,
                            should_quit=False  # This will only lead to errors after the Runtime of `SimPrepper`, so we let it slide...
                            ),  
                  # ...
                  ]


#TODO: Rename to ConfigChecker
# This is the workhorse here, but in the best base, this never needs to be touched
class SetupChecker:
    sim_setup          : SimSetup | None = None
    list_of_checks     : list[CheckItem] = []
    #TODO: instead of verbose flag, parse logging level...?
    verbose            : bool            = False
    _results_of_checks : list            = []
    _has_finished      : bool            = False
    logging                              = None

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
