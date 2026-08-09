import logging
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
    name        : str = ""
    explanation : str = ""
    should_quit : bool = False
    condition = None
    has_passed = None
    result      : str = ""
    def __init__(self, name, condition, explanation, should_quit):
        self.name = name
        self.explanation = explanation
        self.should_quit = should_quit
        self.condition = condition
        return None

    def perform_check(self, sim_setup):
        test_result = self.condition(sim_setup)
        self.has_passed = test_result
        result = {True:"Passed", False:"Failed"}.get(test_result)
        self.result=result
        return result

# Example for devs:
template_item = CheckItem(name="Mutually exclusive sections: [propA] and [propB]",
                          condition=lambda sim_setup: 
                          not all((sim_setup["propA"]==True, 
                                   sim_setup["propB"]==True)),
                          explanation="Dear User, propA and propB are mutually exclusive...",
                          should_quit=True  # if this will definitely lead to an erroneous setup, instruct termination of simprepper
                          )
        

# Actual list of checks for runtime
LIST_OF_CHECKS = [CheckItem(name="Mutually exclusive sections: [Membrane box] and [Water box]",
                            condition=lambda sim_setup: 
                            not all((sim_setup["water_box"]==True, 
                                     sim_setup["membrane_box"]==True)),
                            explanation="This is what the user should see",
                            should_quit=True  # if this will definitely lead to an erroneous setup, instruct termination of simprepper
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

    def __init__(self, 
                 sim_setup, 
                 list_of_checks=LIST_OF_CHECKS,
                 verbose=False):
        self.sim_setup = sim_setup
        self.list_of_checks = list_of_checks
        self.verbose = verbose
        return None

    def run_all_checks(self):
        if self.verbose:
            logging.info("Checking Setup for problems...")
        for check in self.list_of_checks:
            check.perform_check(self.sim_setup)
        self._has_finished = True
        return None

    def get_report(self):
        if not self._has_finished:
            self.run_all_checks()
        if self.verbose:  # print all checks, independent of results
            for check in self.list_of_checks:    
                logging.info(f"Test result: >{check.name}< is: {check.result}")

        for check in self.list_of_checks:
            # TODO: check.explanation should be printed here, somehwere...!
            if not(check.has_passed):
                logging.warning(f"{check.name} has not passed!")
                if check.should_quit:
                    quit_with_error()
        return None
