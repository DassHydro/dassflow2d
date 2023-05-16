#!/usr/bin/env python3
"""
Functions to run a notebook
"""

# Wee need nbformat to run validation of the notebooks
import sys
from argparse import ArgumentParser


def run_notebook(filename, timeout, update_nb=False):
    """
    Run a notebook

    @param filename (string) Name of the notebook to run
    @param timeout (float) Time after which a notebook is killed
    @param update_nb (boolean) Overwrite file with results notebook
    """
    try:
        import nbformat
        from nbconvert.preprocessors import ExecutePreprocessor
    except ImportError:
        print("    Doing nothing nbformat is not avaialable")
        return

    # Reading the file
    with open(filename) as f:
        nbf = nbformat.read(f, as_version=4)

    # Configuring execution noe
    epr = ExecutePreprocessor(timeout=timeout, kernel_name="python3")

    # Running it

    epr.preprocess(nbf, {'metadata': {'path': './'}})

    if update_nb:
        # Saving ouput
        with open(filename, 'wt') as f:
            nbformat.write(nbf, f)

def main():
    """ Main function """
    parser = ArgumentParser()
    parser.add_argument(
        "filename", default="",
        help="name of the notebook to run")
    # time out for the execution
    parser.add_argument(
        "-t", "--timeout",
        dest="timeout", type=int, default=600,
        help="name of the output file also defines the output format")
    parser.add_argument(
        "--update-notebook",
        dest="store_true", default=False,
        help="Update notebook file with the run one")

    args = parser.parse_args()
    run_notebook(args.filename, args.timeout)

if __name__ == "__main__":
    main()
