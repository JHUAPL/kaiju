#!/usr/bin/env python


"""Common code for MAGE regression tests

This module provides common code used by scripts in the MAGE
regression test suite.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import argparse
import os
import subprocess
import sys

# Import 3rd-party modules.
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

# Import project modules.


# Module constants

# Name of Slack channels to use as message target
SLACK_CHANNEL_NAME = '#kaijudev'
SLACK_TEST_CHANNEL_NAME = '#kaijudev-testing'


def read_strip_blanks_and_comments(path):
    """Read a text file, stripping comments and blank lines.

    Read the specified text file, and strip out:
    0. Leading and trailing whitespace
    1. Lines containing only whitespace
    2. Lines containing comments (First non-whitespace character is '#')
    3. In-line comments (everything after '#' on a line) and resiual
       whitespace.

    Parameters
    ----------
    path : str
        Path to text file to read

    Returns
    -------
    lines : list of str
        List containing strings for each line in the file, with comments and
        blank lines removed

    Raises
    ------
    None
    """
    # Read the file.
    lines = open(path, encoding='ut-8').readlines()

    # Strip whitespace.
    lines = [line.strip() for line in lines]

    # Strip blank lines.
    lines = [line for line in lines if len(line) > 0]

    # Strip comment lines and residual trailing whitespace.
    lines = [line for line in lines if not line.startswith('#')]

    # Strip in-line comments.
    newlines = []
    for line in lines:
        comment_pos = line.find('#')
        if comment_pos > -1:
            line = line[:comment_pos].rstrip()
        newlines.append(line)
    lines = newlines

    # Return the file contents as a list of strings.
    return lines


def create_command_line_parser(description):
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    description : str
        Description of script

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '--account', default=os.environ['DERECHO_TESTING_ACCOUNT'],
        help='PBS account to use for testing (default: %(default)s)'
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


def run_mage_test_script(test_script, args):
    """Run a single MAGE test script.

    Run a single MAGE test script.

    Parameters
    ----------
    test_script : str
        Path to test script to run
    args : dict
        Options provided on command line

    Returns
    -------
    cproc : subprocess.CompletedProcess
        Object containing results of running the command

    Raises
    ------
    None
    """
    # Assemble command-line flags to pass to the individual test
    # scripts.
    test_script_options = ''
    test_script_options += f" --account {args.account}"
    if args.all:
        test_script_options += ' -a'
    if args.debug:
        test_script_options += ' -d'
    if args.force:
        test_script_options += ' -f'
    if args.loud:
        test_script_options += ' -l'
    if args.test:
        test_script_options += ' -t'
    if args.verbose:
        test_script_options += ' -v'

    # Run the test script.
    cmd = f"python {test_script} {test_script_options}"
    cproc = subprocess.run(cmd, shell=True, check=True)
    return cproc


def read_build_module_list_file(list_file):
    """Read a MAGE build module list file

    Read a MAGE build module list file.

    Parameters
    ----------
    list_file : str
        Path to MAGE build module list/file list file to read

    Returns
    -------
    module_or_file_names : list of str
        List of modules, or list of module list files
    cmake_environment : str
        Environment variables and values to set when invoking cmake with this
        module set
    cmake_options : str
        Command-line options for cmake when building MAGE with this module set

    Raises
    ------
    None
    """
    # Read the file, ignoring blank lines and comments.
    lines = read_strip_blanks_and_comments(list_file)

    # Extract the optional cmake environment variable settings,
    cmake_environment = ''
    label = 'CMAKE_ENV='
    if lines[0].startswith(label):
        cmake_environment = lines[0][len(label):].rstrip()
        lines.pop(0)  # Remove cmake environment line.

    # Extract the optional cmake command-line options,
    cmake_options = ''
    label = 'CMAKE_OPTIONS='
    if lines[0].startswith(label):
        cmake_options = lines[0][len(label):].rstrip()
        lines.pop(0)  # Remove cmake options line.

    # Save the remaining lines as a module or file list.
    module_or_file_names = lines

    # Return the file contents.
    return module_or_file_names, cmake_environment, cmake_options


# -----------------------------------------------------------------------------

# Git utilities


def git_get_branch_name():
    """Fetch the name of the current git branch,

    Fetch the name of the current git branch.

    Parameters
    ----------
    None

    Returns
    -------
    git_branch_name : str
        Name of the current git branch

    Raises
    ------
    None
    """
    cmd = 'git symbolic-ref --short HEAD'
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    git_branch_name = cproc.stdout.rstrip()

    # Return the git branch name,
    return git_branch_name


def git_pull():
    """Pull the current branch from the git repository.

    Pull the current branch from the git repository.

    Parameters
    ----------
    None

    Returns
    -------
    git_pull_output : str
        Output from git pull command

    Raises
    ------
    None
    """
    cmd = 'git pull'
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    git_pull_output = cproc.stdout.rstrip()

    # Return the git pull output.
    return git_pull_output


# -----------------------------------------------------------------------------

# Slack utilities


def slack_create_client():
    """Create a client for Slack communication.

    Create a client for Slack communication.

    Parameters
    ----------
    None

    Returns
    -------
    slack_client : slack_sdk.WebClient
        Client for Slack communication

    Raises
    ------
    None
    """
    # Get the Slack API token
    slack_token = os.environ['SLACK_BOT_TOKEN']

    # Create the Slack client.
    slack_client = WebClient(token=slack_token)

    # Return the Slack client.
    return slack_client


def slack_send_message(slack_client, message, thread_ts=None, is_test=False):
    """Send a message to Slack.

    Send a message to Slack. Errors during message sending are not considered
    to be fatal. Errors are caught, an error message is printed, and the
    program continues normally.

    Parameters
    ----------
    slack_client : slack_sdk.WebClient
        Client for Slack communication
    message : str
        Message to send to Slack
    thread_ts : XXX, default None
        Set to desired Slack thread identifier (timestamp), if any
    is_test : bool
        If True, use the testing channel as the message target.

    Returns
    -------
    response : XXX
        Response from Slack API when posting the message.

    Raises
    ------
    None
    """
    if is_test:
        channel = SLACK_TEST_CHANNEL_NAME
    else:
        channel = SLACK_CHANNEL_NAME
    try:
        response = slack_client.chat_postMessage(
            channel=channel,
            thread_ts=thread_ts,
            text=message,
        )
    except SlackApiError as e:
        print('Sending message to Slack failed.', file=sys.stderr)
        response = e.response
        print(f"response = {response}", file=sys.stderr)
    return response


def slack_send_image(slack_client, image_file_path, initial_comment='',
                     thread_ts=None, is_test=False):
    """Send an image file to Slack.

    Send an image file to Slack.

    Parameters
    ----------
    slack_client : slack_sdk.WebClient
        Client for Slack communication
    image_file_path : str
        Path to image file to send to Slack
    initial_comment : str
        Comment to include with image, default ''
    thread_ts : XXX, default None
        Set to desired Slack thread identifier (timestamp), if any
    is_test : bool
        If True, use the testing channel as the message target.

    Returns
    -------
    response : XXX
        Response from Slack API when posting the image.

    Raises
    ------
    None
    """
    if is_test:
        channel = SLACK_TEST_CHANNEL_NAME
    else:
        channel = SLACK_CHANNEL_NAME
    try:
        response = slack_client.files_upload(
            channels=channel,
            thread_ts=thread_ts,
            file=image_file_path,
            initial_comment=initial_comment,
        )
    except SlackApiError as e:
        print('Sending image to Slack failed.', file=sys.stderr)
        response = e.response
        print(f"response = {response}", file=sys.stderr)
    return response
