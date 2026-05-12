import argparse

parser = argparse.ArgumentParser(
                prog="simprepper",
                description="\n".join(["Module to prepare OpenMM system starting from protein and optional ligand files.",
                                       "Exports the prepared systems to amber and gromacs formats."
                                       ]),
                formatter_class=argparse.ArgumentDefaultsHelpFormatter  # show defaults in `--help`
            )

parser.add_argument("--debug",
                    help="Print debug messenges (mostly for developers)",
                    action="store_true", default=False)

parser.add_argument("-l", "--lig",
                    help="Optional SDF file of the ligand",
                    required=False, default=None,
                    )
parser.add_argument("-r", "--rec",
                    help="PDB file of the receptor",
                    required=True,
                    )
parser.add_argument("-L", "--log-level",
                    help="Choose the logging level to show",
                    choices=["debug", "info", "warning", "error", "critical"],
                    default="info",
                    required=False,
                    )
# %% various

parser.add_argument("--ligand_ff",
                    help="Choose force field for ligand",
                    choices=["espaloma", "SMIRNOFF", "GAFF"],
                    default="espaloma",
                    required=False
                    )
parser.add_argument('--temperature', 
                    help='Simulation temperature', 
                    type=float, default=300.0, required=False)

parser.add_argument('--Hmass', 
                    help='Mass of Hydrogens (for hydrogen-mass repartitioning)', 
                    type=float, default=4.0)

# %% simulation box
box_args = parser.add_argument_group(description="\nEither set the box length explicitly, or let it be determined automatically via box padding...")
box_args_mut_ex = box_args.add_mutually_exclusive_group(required=False)
box_args_mut_ex.add_argument('--box_padding', help='Padding for water box (in nm)', default=2.0, required=False, type=float)
box_args_mut_ex.add_argument('--box_length', help='Length of cubic box (in nm), e.g., 5.0', default=None, type=float)
parser.add_argument('--box_shape', 
                    help='Shape of box', 
                    choices=['cube', 'dodecahedron'], 
                    default='cube', required=False)


