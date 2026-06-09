simprepper -r 181L_clean.pdb -l BNZ.sdf -L info --ligand_ff GAFF

# this will cause an error:
# simprepper -r 181L_clean_noTER.pdb -l BNZ.sdf -L info --ligand_ff GAFF

# test for ligand extension
# simprepper -r 181L_clean.pdb -l BNZ.mol2 -L info --ligand_ff GAFF
