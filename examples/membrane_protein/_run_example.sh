# generate a config file for a membrane protein system
simprepper-example-config --membrane # By default writes the config to default_config.ini

simprepper -r 1afo_opm_no_dum.pdb -L info -i default_config.ini
