import argparse
from pathlib import Path
from argparse import RawTextHelpFormatter

def get_args():
    """
    Function defines all possible arguments for the script
    """

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)


    parser.add_argument(
        "--i",
        nargs=1,
        metavar=("file1"),
        help="Path to input XYZ file containing the molecular structure"
    )

    parser.add_argument(
        "--test",
        # User can choose multiple tests from the list
        nargs="+",
        choices=["method_basis_combinations"],
        help="""Choose tests to perform:\n
              method_basis_combinations: Determine viable method and basis set combinations for the molecule"""
    )

    args = parser.parse_args()
    return args