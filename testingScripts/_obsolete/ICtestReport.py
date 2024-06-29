#!/usr/bin/env python

"""Generate reports for MAGE initial condition tests.

This script generates reports for the MAGE initial condition tests.

Authors
-------
Jeff Garretson
Eric Winter
"""


# Import standard modules.
import argparse
import datetime
import glob
import os
import re
import subprocess
import sys

# Import 3rd-party modules.
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = 'Generate reports for MAGE initial condition testing'


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
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--account', default=None,
        help='PBS account to use for testing (default: %(default)s)'
    )
    parser.add_argument(
        '--all', '-a', action='store_true',
        help='Run all tests (default: %(default)s).'
    )
    parser.add_argument(
        '--debug', '-d', action='store_true',
        help='Print debugging output (default: %(default)s).'
    )
    parser.add_argument(
        '--force', '-f', action='store_true',
        help='Force all tests to run (default: %(default)s).'
    )
    parser.add_argument(
        '--loud', '-l', action='store_true',
        help='Enable loud mode (post results to Slack) (default: %(default)s).'
    )
    parser.add_argument(
        '--test', '-t', action='store_true',
        help='Enable testing mode (default: %(default)s).'
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Print verbose output (default: %(default)s).'
    )
    return parser


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
    account = args.account
    doAll = args.all
    debug = args.debug
    forceRun = args.force
    beLoud = args.loud
    isTest = args.test
    verbose = args.verbose

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")

    #--------------------------------------------------------------------------

    # Set up for communication with Slack.

    # Get the Slack API token
    slack_token = os.environ["SLACK_BOT_TOKEN"]
    if debug:
        print(f"slack_token = {slack_token}")

    # Create the Slack client.
    slack_client = WebClient(token=slack_token)
    if debug:
        print(f"slack_client = {slack_client}")

    #--------------------------------------------------------------------------

    # Determine the path to the MAGE installation to use for testing.

    # Fetch the path to this running script.
    called_from = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print(f"called_from = {called_from}")

    # Assume this script is in a subdirectory of the kaiju directory.
    os.chdir(called_from)
    os.chdir('..')

    # Use this directory as the home directory for testing.
    home = os.getcwd()
    if debug:
        print(f"home = {home}")
    if verbose:
        print('I am the initial condition check report script. '
              'This is my current home directory:')
        print(home)

    #--------------------------------------------------------------------------

    # Determine the current git branch.
    cmd = 'git symbolic-ref --short HEAD'
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    git_branch = cproc.stdout.rstrip()
    if debug:
        print(f"git_branch = {git_branch}")
    
    #--------------------------------------------------------------------------

    # Go to the folder containing the initial condition builds.
    os.chdir(home)
    if debug:
        print(f"Now in directory {os.getcwd()}.")
    os.chdir('ICBuilds')
    if debug:
        print(f"Now in directory {os.getcwd()}.")

    # Get a list of subdirectories, one per module set/initial
    # condition file.
    initial_condition_build_directories = os.listdir(os.getcwd())
    if debug:
        print("initial_condition_build_directories = "
              f"{initial_condition_build_directories}")

    # Check each subdirectory to ensure gamera.x was built.
    pattern = '^gamera_(.*)_(.+)$'
    ic_re = re.compile(pattern)
    for initial_condition_build_directory in initial_condition_build_directories:
        gamera_path = os.path.join(initial_condition_build_directory,
                            'bin', 'gamera.x')
        gamera_built_ok = os.path.exists(gamera_path)
        if gamera_built_ok:
            continue
        if verbose:
            print('Found a bad one!')
        m = ic_re.match(dir_name)
        initial_condition_name, module_set_name = (m.group(1), m.group(2))
        message = (
            f"Module set {module_set_name} failed to build gamera.x "
            f"for initial condition {initial_condition_name}."
        )

        # If this is a test run, don't post to Slack.
        if isTest:
            pass
        else:
            # If loud, or an error occurred, send Slack message.
            if beLoud:
                try:
                    response = slack_client.chat_postMessage(
                        channel="#kaijudev",
                        text=message,
                    )
                except SlackApiError as e:
                    # You will get a SlackApiError if "ok" is False
                    # str like 'invalid_auth', 'channel_not_found'
                    assert e.response["error"]

        # Send message to stdout.
        print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
