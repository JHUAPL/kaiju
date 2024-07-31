#!/usr/bin/env python

"""Null MAGE test.

This script does not run any tests, and only sends a message too Slacak if
requested.

NOTE: These tests are performed on a load-balance-assigned login node on
derecho. No PBS job is submitted.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
# import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
import common


# Program constants

# Program description.
DESCRIPTION = 'Script for null MAGE test'

# # Root of directory tree for this set of tests.
# MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# # Directory for build tests
# BUILD_TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'buildTest')

# # Path to directory to use for building executable list
# EXECUTABLE_LIST_BUILD_DIRECTORY = os.path.join(BUILD_TEST_DIRECTORY,
#                                                'build_executable_list')

# # Home directory of kaiju installation
# KAIJUHOME = os.environ['KAIJUHOME']

# # Path to directory containing the test scripts
# TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# # Path to directory containing module lists
# MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
#                                      'mage_build_test_modules')

# # Path to module list file to use when generating the list of executables
# # Use a module set without MKL.
# EXECUTABLE_LIST_MODULE_LIST = os.path.join(MODULE_LIST_DIRECTORY, '04.lst')

# # Path to file containing list of module sets to use for build tests
# BUILD_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'build_test.lst')

# # Prefix for naming build test directories
# BUILD_TEST_DIRECTORY_PREFIX = 'buildTest_'

# # Name of subdirectory of current build subdirectory containing binaries
# BUILD_BIN_DIR = 'bin'

# # Branch or commit (or tag) used for testing.
# BRANCH_OR_COMMIT = os.environ['BRANCH_OR_COMMIT']


def main():
    """Begin main program.

    This is the main program code.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    subprocess.CalledProcessError
        If an exception occurs in subprocess.run()
    """
    # Set up the command-line parser.
    parser = common.create_command_line_parser(DESCRIPTION)

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    be_loud = args.loud
    is_test = args.test
    slack_on_fail = args.slack_on_fail
    verbose = args.verbose

    # -------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # -------------------------------------------------------------------------

    # # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        'Test results are on `derecho` in '
        f"`{os.environ['MAGE_TEST_SET_ROOT']}`.\n"
    )

    # Summarize the test results.
    test_report_summary_string = 'This was a null test.'

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If loud mode is on, post report to Slack.
    if be_loud:
        slack_client = common.slack_create_client()
        if debug:
            print(f"slack_client = {slack_client}")
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=is_test
        )
        if debug:
            print(f"slack_response_summary = {slack_response_summary}")
        thread_ts = slack_response_summary['ts']
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_details_string, thread_ts=thread_ts,
            is_test=is_test
        )
        if debug:
            print(f"slack_response_summary = {slack_response_summary}")

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    main()
