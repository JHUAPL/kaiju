#!/usr/bin/env python


"""Master script for automated MAGE regression testing.

This script controls automated MAGE regression testing on
derecho. This script should be run in a cron job using the testing
setup script run_mage_test_00.sh, which must run
kaiju/scripts/setupEnvironment.sh before running this script, to
ensure that the KAIJUHOME environment variable is set.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Master script for MAGE regression testing'

# Subdirectory of KAIJUHOME containing the test scripts.
KAIJU_TEST_SCRIPTS_DIRECTORY = 'testingScripts'


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
    run_all_tests = args.all
    debug = args.debug
    force_tests = args.force
    verbose = args.verbose

    if verbose:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # #--------------------------------------------------------------------------

    # # Move to the MAGE installation directory.
    # kaiju_home = os.environ['KAIJUHOME']
    # os.chdir(kaiju_home)

    # #--------------------------------------------------------------------------

    # # Check to see if any changes have been made in the MAGE
    # # repository.

    # # Find the current branch.
    # if verbose:
    #     print(f"Fetching git branch name for directory {kaiju_home}.")
    # git_branch_name = common.git_get_branch_name()
    # if debug:
    #     print(f"git_branch_name = {git_branch_name}")

    # # Pull any changes from the repository.
    # if verbose:
    #     print(f"Pulling changes from repository on branch {git_branch_name}.")
    # git_pull_output = common.git_pull()
    # if debug:
    #     print(git_pull_output)

    # # If there are no changes to the code, and the forced-test flag is
    # # not set, abort this test.
    # if git_pull_output == 'Already up to date.' and not force_tests:
    #     print(f"No test today. Branch {git_branch_name} is already up to date!")
    #     exit(0)

    # #--------------------------------------------------------------------------

    # # Run the tests.

    # # Run the tests.
    # if run_all_tests:
    #     if verbose:
    #         print('Running all tests.')
    #     test_scripts = [
    #         # 'buildTest.py',
    #         # 'unitTest.py',
    #         # 'intelChecks.py',
    #         # 'ICtest.py',
    #         # 'pyunitTest.py',
    #         # 'weeklyDash.py',
    #         # 'weeklyDashReport.py',
    #     ]
    # else:
    #     if verbose:
    #         print('Running typical tests.')
    #     test_scripts = [
    #         # 'buildTest.py',
    #         # 'unitTest.py',
    #         # 'ICtest.py',
    #         # 'pyunitTest.py',
    #     ]
    # for test_script in test_scripts:
    #     path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
    #                         test_script)
    #     if verbose:
    #         print(f"Running test script {test_script}.")
    #     common.run_mage_test_script(path, args)

    #--------------------------------------------------------------------------

    if verbose:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
