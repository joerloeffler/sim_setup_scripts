simprepper -r 181L_only_prot.pdb -L info

# example with the ini file
simprepper-example-config # By default writes the config to default_config.ini
simprepper -r 181L_only_prot.pdb -L info -i default_config.ini

# Example to play with forcefields
# simprepper -r 181L_only_prot.pdb -L info --protein_ff amber14/protein.ff14SB.xml --water_ff  amber14/tip3pfb.xml --ion_ff amber/tip3p_HFE_multivalent.xml --lipid_ff amber14/lipid17.xml

