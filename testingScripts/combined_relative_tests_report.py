#!/usr/bin/env python


"""Create the combined MAGE relative test report.

This script creates the combined MAGE relative test report for the serial and
MPI versions of the test.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import copy
import datetime
import glob
import os
import platform
# import re
# import subprocess
import sys

# Import 3rd-party modules.
# from astropy.time import Time
# import matplotlib as mpl
# import matplotlib.dates as mdates
# import matplotlib.pyplot as plt

# Import project modules.
import common
# import kaipy.kaiH5 as kh5
# import kaipy.kaiViz as kv


# Program constants

# Program description.
DESCRIPTION = (
    "Create a combined report for the serial and MPI relative tests."
)


# Default values for command-line arguments.
DEFAULT_ARGUMENTS = {
    "debug": False,
    "loud": False,
    "slack_on_fail": False,
    "test": False,
    "verbose": False,
    "mpi1": "",
    "mpi2": "",
    "serial1": "",
    "serial2": "",
}

# # Root of directory tree for all tests.
# MAGE_TEST_ROOT = os.environ["MAGE_TEST_ROOT"]

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ["MAGE_TEST_SET_ROOT"]

# Directory for relative tests
TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, "compTest")

# Glob pattern for individual relative test directories.
COMP_TESTS_DIRECTORY_GLOB_PATTERN = "compTest_*"

# # Regular expression for git hash read from weekly dash output log.
# GIT_HASH_PATTERN = "Git hash   = (.{8})"

# Name of subdirectory containing binaries and test results.
BIN_DIR = "bin"

# # Name of file containg PBS job IDs.
# JOB_LIST_FILE = "jobs.txt"

# String naming branch or commit used in this test.
BRANCH_OR_COMMIT = os.environ["BRANCH_OR_COMMIT"]


def create_command_line_parser():
    """Create the command-line parser.

    Create the parser for the command line.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line parser for this script.

    Raises
    ------
    None
    """
    parser = common.create_command_line_parser(DESCRIPTION)
    parser.add_argument(
        "mpi1", type=str, default=DEFAULT_ARGUMENTS["mpi1"],
        help="Folder for first MPI case."
    )
    parser.add_argument(
        "mpi2", type=str, default=DEFAULT_ARGUMENTS["mpi2"],
        help="Folder for second MPI case."
    )
    parser.add_argument(
        "serial1", type=str, default=DEFAULT_ARGUMENTS["serial1"],
        help="Folder for first serial case."
    )
    parser.add_argument(
        "serial2", type=str, default=DEFAULT_ARGUMENTS["serial2"],
        help="Folder for second serial case."
    )
    return parser


def combined_relative_tests_report(args: dict):
    """Create a combined test report for relativeTests.py.

    Create a combined test report for relativeTests.py.

    Parameters
    ----------
    args: dict
        Dictionary of command-line options.

    Returns
    -------
    int 0 on success

    Raises
    ------
    AssertionError
        If an invalid argument is provided.
    """
    # Set defaults for command-line options, then update with values passed
    # from the caller.
    local_args = copy.deepcopy(DEFAULT_ARGUMENTS)
    local_args.update(args)
    args = local_args

    # Local convenience variables.
    debug = args["debug"]
    loud = args["loud"]
    slack_on_fail = args["slack_on_fail"]
    test = args["test"]
    verbose = args["verbose"]
    mpi1 = args["mpi1"]
    mpi2 = args["mpi2"]
    serial1 = args["serial1"]
    serial2 = args["serial2"]

    # Validate arguments.
    assert len(serial1) > 0
    assert len(serial2) > 0
    assert len(mpi1) > 0
    assert len(mpi2) > 0

    # ------------------------------------------------------------------------

#     # Add additional arguments
#     parser.add_argument(
#         "-d1", type=str, help="Folder for the first case"
#     )
#     parser.add_argument(
#         "-d2", type=str, help="Folder for the second case"
#     )
#     parser.add_argument(
#         "-cn", type=str, help="Name of this post process case"
#     )

#     # Parse the command-line arguments.
#     args = parser.parse_args()
#     if args.debug:
#         print(f"args = {args}")
#     debug = args.debug
#     be_loud = args.loud
#     # slack_on_fail = args.slack_on_fail
#     is_test = args.test
#     verbose = args.verbose
#     d1 = args.d1
#     d2 = args.d2
#     cn = args.cn

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Move to the top-level weekly dash directory.
    os.chdir(TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Get list of comparison test directories.
    test_directories = glob.glob(COMP_TESTS_DIRECTORY_GLOB_PATTERN)
    if debug:
        print(f"comp_test_directories = {test_directories}")

    # <HACK>
    # Use only the first directory for now.
    test_directory = test_directories[0]
    # </HACK>

    # ------------------------------------------------------------------------

    # Read results from the latest run.
    if verbose:
        print(f"Reading results for latest run in {test_directory}.")

    # Go to test folder
    os.chdir(test_directory)

    # Move down to the build directory
    os.chdir(BIN_DIR)

    # -----------------------------------------------------------------------

    # Detail the test results.
    test_report_details_string = ""
    test_report_details_string += (
        f"Test results are on `derecho` in `{test_directory}`.\n"
    )

    # Summarize the test results.
    test_report_summary_string = (
        "Comparative test result plots complete on branch "
        f"`{BRANCH_OR_COMMIT}`."
    )

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If loud mode is on, post results to Slack.
    if loud:
        slack_client = common.slack_create_client()
        if debug:
            print(f"slack_client = {slack_client}")
        slack_response = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=test
        )
        if slack_response["ok"]:
            parent_ts = slack_response["ts"]
            message = (
                "All results are from a Double Resolution Run using"
                " various Gamera and Voltron settings.\n"
            )
            slack_response = common.slack_send_message(
                slack_client, message, thread_ts=parent_ts, is_test=test
            )
            message = f"This result compared `{mpi1}` and `{mpi2}`.\n"
            slack_response = common.slack_send_image(
                slack_client, os.path.join("vidData", "MpiRestartComp.mp4"),
                initial_comment=message,
                thread_ts=parent_ts, is_test=test
            )
            message = f"This result compared `{serial1}` and `{serial2}`.\n"
            slack_response = common.slack_send_image(
                slack_client, os.path.join("vidData", "SerialMpiComp.mp4"),
                initial_comment=message,
                thread_ts=parent_ts, is_test=test
            )
        else:
            print("Failed to post parent message and images to Slack.")

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")

    # Return normally.
    return 0


def main():
    """Driver for command-line version of code."""
    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")

    # Convert the arguments from Namespace to dict.
    args = vars(args)

    # Call the main program code.
    return_code = combined_relative_tests_report(args)
    sys.exit(return_code)


if __name__ == "__main__":
    main()
