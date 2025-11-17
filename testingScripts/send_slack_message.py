#!/usr/bin/env python

"""Send a message to Slack.

Send a message to Slack.

Authors
-------
Eric Winter
"""


# Import standard modules.
import datetime
# import glob
import os
import sys

# # Import 3rd-party modules.

# Import project modules.
import common


# Program constants

# Program description.
DESCRIPTION = "Send a message to Slack."

# # Root of directory tree for this set of tests.
# MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# # Directory for unit tests
# UNIT_TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'unitTest')

# # glob pattern for naming unit test directories
# UNIT_TEST_DIRECTORY_GLOB_PATTERN = 'unitTest_*'

# # Name of build subdirectory containing binaries
# BUILD_BIN_DIR = 'bin'

# # Name of file containing job IDs for each unit test directory.
# JOB_ID_LIST_FILE = 'jobs.txt'


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.

    Raises
    ------
    None
    """
    parser = common.create_command_line_parser(DESCRIPTION)
    parser.add_argument(
        "message",
        default="",
        help="Message to send to Slack (default: %(default)s)"
    )
    return parser


def send_slack_message(args: dict = None):
    """Send a message to Slack.

    Send a message to Slack.

    Parameters
    ----------
    args : dict
        Dictionary of program options.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Local convenience variables.
    debug = args["debug"]
    be_loud = args["loud"]
    slack_on_fail = args["slack_on_fail"]
    is_test = args["test"]
    verbose = args["verbose"]
    message = args["message"]

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Create the Slack client.
    slack_client = common.slack_create_client()
    slack_response_summary = common.slack_send_message(
        slack_client, message, is_test=is_test
    )

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


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
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")

    # ------------------------------------------------------------------------

    # Call the main program logic. Note that the Namespace object (args)
    # returned from the option parser is converted to a dict using vars().
    send_slack_message(vars(args))


if __name__ == "__main__":
    main()
