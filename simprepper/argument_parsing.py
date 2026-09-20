import argparse

def str2bool(v) -> bool:
    """This function enables the translation of user input to actual boolean values
    otherwise bool("false") would return `True`,
    which is really not what the user wants...
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


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

parser.add_argument("--verbose",
                    help="Print some additional info during runtime",
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
parser.add_argument("-i", "--ini",
                    help="ini-file with simulation parameters",
                    required=False,
                    )
# %% different output

parser.add_argument("--gmx_output",
                    help="Flag: Generate files for gromacs",
                    type=str2bool, default=True, dest="should_export_gmx")

parser.add_argument("--amber_output",
                    help="Flag: Generate files for amber",
                    type=str2bool, default=True, dest="should_export_amber")

parser.add_argument("--openmm_output",
                    help="Flag: Generate files for openmm",
                    type=str2bool, default=True, dest="should_export_openmm")

