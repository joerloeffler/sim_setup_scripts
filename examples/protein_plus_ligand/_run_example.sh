simprepper -r 181L_clean.pdb -l BNZ.sdf -L info --ligand_ff GAFF

# example with the ini file
simprepper-to-ini # By default writes the config to default_config.ini
simprepper -r 181L_clean.pdb -l BNZ.sdf -L info -i default_config.ini

# this will cause an error:
# simprepper -r 181L_clean_noTER.pdb -l BNZ.sdf -L info --ligand_ff GAFF

# test for ligand extension
# simprepper -r 181L_clean.pdb -l BNZ.mol2 -L info --ligand_ff GAFF
