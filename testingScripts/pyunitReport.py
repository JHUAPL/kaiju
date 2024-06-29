#!/usr/bin/env python

"""Create a report for the results of the MAGE python unit tests.

Create a report for the results of the MAGE python unit tests.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Generate a report for MAGE python unit test results.'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Path to directory for running python unit tests
PYTHON_UNIT_TEST_DIRECTORY = os.path.join(KAIJUHOME, 'pytests')

# Name of python unit test log file
PYTHON_UNIT_TEST_LOG_FILE = 'kaiju-pyunit.txt'


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
    None
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
    verbose = args.verbose

    # -------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # -------------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # -------------------------------------------------------------------------

    # Move to the python unit tests directory.
    os.chdir(PYTHON_UNIT_TEST_DIRECTORY)

    # -------------------------------------------------------------------------

    # Detail the test results
    test_details_message = ''
    test_report_details_string += (
        f"Test results are in {os.getcwd()}.\n"
    )

    # Check that the python unit tests completed.
    has_fail = False
    has_error = False
    has_pass = False
    if os.path.exists(PYTHON_UNIT_TEST_LOG_FILE):
        test_details_message += 'Python unit tests ran to completion.\n'
        with open(PYTHON_UNIT_TEST_LOG_FILE, encoding='utf-8') as f:
            lines = f.readlines()
        last_line = lines[-1]
        if 'fail' in last_line:
            has_fail = True
        if 'error' in last_line:
            has_error = True
        if 'passed' in last_line:
            has_pass = True
        if has_error:
            test_details_message += 'Python unit tests error detected.\n'
        if has_fail:
            test_details_message += 'Python unit tests: *FAILED*\n'
        if has_pass and not has_error and not has_fail:
            test_details_message += 'Python unit tests: *PASSED*\n'
        if not has_pass and not has_error and not has_fail:
            test_details_message += (
                'Unexpected error occured during python unit tests.\n'
            )
    else:
        test_details_message += 'Python unit tests did not run to completion.\n'
        has_pass = False

    # Summarize the test results.
    test_summary_message = (
        'Summary of python unit test results from `pyunitReport.py`: '
    )
    if not has_pass:
        test_summary_message += '*FAILED*\n'
    else:
        test_summary_message += '*ALL PASSED*\n'

    # If loud mode is on, post report to Slack.
    if be_loud:
        test_summary_message += 'Details in thread for this messsage.\n'
        slack_response_summary = common.slack_send_message(
            slack_client, test_summary_message, is_test=is_test
        )
        if slack_response_summary['ok']:
            thread_ts = slack_response_summary['ts']
            slack_response_details = common.slack_send_message(
                slack_client, test_details_message, thread_ts=thread_ts,
                is_test=is_test
            )
        else:
            print('*ERROR* Unable to post test result summary to Slack.')

    # -------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
