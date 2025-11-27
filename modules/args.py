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


    args = parser.parse_args()
    return args