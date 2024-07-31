#!/usr/bin/env python

"""Run MAGE initial condition build tests.

This script runs a series of initial condition build tests of the MAGE
software. Each set of initial conditions provided in the kaiju source
tree is tested with each set of modules used for these tests. The module
sets are listed in files under:

$KAIJUHOME/testingScripts/mage_build_test_modules

This script reads the file initial_condition_build_test.lst from this
directory, and uses the contents as a list of module list files to use for
MAGE iinitial condition build tests.

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
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE initial condition build testing'

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for initial condition build tests
INITIAL_CONDITION_BUILD_TEST_DIRECTORY = os.path.join(
    MAGE_TEST_SET_ROOT, 'ICtest'
)

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'mage_build_test_modules'
)

# Path to file containing names of module lists to use for initial
# condition build tests
# Path to module list file to use when generating the list of executables
INITIAL_CONDITION_BUILD_TEST_LIST_FILE = os.path.join(
    MODULE_LIST_DIRECTORY, 'initial_condition_build_test.lst'
)

# Path to directory containing initial condition source code
INITIAL_CONDITION_SRC_DIRECTORY = os.path.join(
    KAIJUHOME, 'src', 'gamera', 'ICs'
)

# Prefix for subdirectory names for individual builds.
INITIAL_CONDITION_BUILD_DIR_PREFIX = 'ICtest_'

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'

# Branch or commit (or tag) used for testing.
BRANCH_OR_COMMIT = os.environ['BRANCH_OR_COMMIT']


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

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the initial condition build tests.
    print(f"Creating ${INITIAL_CONDITION_BUILD_TEST_DIRECTORY}.")
    os.mkdir(INITIAL_CONDITION_BUILD_TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    with open(INITIAL_CONDITION_BUILD_TEST_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Get a list of initial conditions to try, ignoring files in the
    # "deprecated" folder.
    initial_condition_paths = []

    for root, _, filenames in os.walk(
        INITIAL_CONDITION_SRC_DIRECTORY
    ):
        if 'deprecated' not in root and 'underdev' not in root:
            for filename in filenames:
                initial_condition_paths.append(os.path.join(root, filename))
    if debug:
        print(f"initial_condition_paths = {initial_condition_paths}")

    # ------------------------------------------------------------------------

    # Initalize test results for all module sets and initial conditions to
    # False (failed).
    test_passed = []
    for _ in module_list_files:
        test_passed.append([False]*len(initial_condition_paths))

    # Define the make command for each build.
    make_cmd = 'make gamera.x'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Build with each initial condition with each set of modules.
    for (i_test, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing initial condition build tests with module set '
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_set_name = {module_set_name}.")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
        if debug:
            print(f"path = {path}")
        module_names, cmake_environment, cmake_options = (
            common.read_build_module_list_file(path)
        )
        if debug:
            print(f"module_names = {module_names}")
            print(f"cmake_environment = {cmake_environment}")
            print(f"cmake_options = {cmake_options}")

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Build with each initial condition.
        for (j_ic, initial_condition_path) in enumerate(
            initial_condition_paths
        ):

            # Extract the name of the initial condition.
            initial_condition_name = os.path.split(
                os.path.splitext(initial_condition_path)[0]
            )[-1]
            if debug:
                print(f"initial_condition_name={initial_condition_name}")

            if verbose:
                print(f"Building with module set {module_set_name} and "
                      f"initial condition {initial_condition_name}.")

            # Make a directory for this build, and go there.
            build_directory = os.path.join(
                INITIAL_CONDITION_BUILD_TEST_DIRECTORY,
                f"{INITIAL_CONDITION_BUILD_DIR_PREFIX}{initial_condition_name}"
                f"_{module_set_name}"
            )
            if debug:
                print(f"build_directory = {build_directory}")
            os.mkdir(build_directory)
            os.chdir(build_directory)

            # Add extra cmake options for initial condition test builds.
            IC_cmake_options = (
                f"{cmake_options} -DGAMIC:FILEPATH={initial_condition_path}"
            )
            if debug:
                print(f"IC_cmake_options = {IC_cmake_options}")

            # Run cmake to build the Makefile.
            if verbose:
                print(
                    'Running cmake to create Makefile for module set'
                    f" {module_set_name},"
                    f" initial condition {initial_condition_name}."
                )
            cmd = (
                f"{module_cmd}; {cmake_environment} cmake {IC_cmake_options}"
                f" {KAIJUHOME} >& cmake.out"
            )
            if debug:
                print(f"cmd = {cmd}")
            try:
                # NOTE: stdout and stderr goes into cmake.out.
                _ = subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(
                    f"ERROR: cmake for module set {module_set_name}, "
                    f"initial condition {initial_condition_name} failed.\n"
                    f"e.cmd = {e.cmd}\n"
                    f"e.returncode = {e.returncode}\n"
                    f"See {os.path.join(build_directory, 'cmake.out')}"
                    ' for output from cmake.\n'
                    "Skipping remaining steps for module set "
                    f"{module_set_name}, initial condition "
                    "{initial_condition_name}.",
                    file=sys.stderr
                )
                continue

            # Run the build.
            if verbose:
                print(
                    'Running make to build kaiju for module set'
                    f" {module_set_name}, initial condition"
                    f" {initial_condition_name}."
                )
            cmd = f"{module_cmd}; {make_cmd} >& make.out"
            if debug:
                print(f"cmd = {cmd}")
            try:
                # NOTE: stdout and stderr go to makeout.
                _ = subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(
                    f"ERROR: make for module set {module_set_name}, "
                    f"initial condition {initial_condition_name} failed.\n"
                    f"e.cmd = {e.cmd}\n"
                    f"e.returncode = {e.returncode}\n"
                    f"See {os.path.join(build_directory, 'make.out')}"
                    ' for output from make.\n'
                    "Skipping remaining steps for module set "
                    f"{module_set_name}, initial condition "
                    "{initial_condition_name}.",
                    file=sys.stderr
                )
                continue

            # Check for gamera.x.
            executable_list = ['gamera.x']
            missing = []
            for executable in executable_list:
                path = os.path.join(build_directory, BUILD_BIN_DIR, executable)
                if not os.path.isfile(path):
                    missing.append(executable)
            if len(missing) > 0:
                for executable in missing:
                    print(f"ERROR: Did not build {executable}.")
            else:
                test_passed[i_test][j_ic] = True

        # End loop over initial conditions
    # End loop over module sets

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        'Test results are on `derecho` in '
        f"`{INITIAL_CONDITION_BUILD_TEST_DIRECTORY}`.\n"
    )
    for (i_test, module_list_file) in enumerate(module_list_files):
        module_set_name = module_list_file.rstrip('.lst')
        for (j_ic, initial_condition_path) in enumerate(
            initial_condition_paths
        ):
            initial_condition_name = os.path.split(
                os.path.splitext(initial_condition_path)[0]
            )[-1]
            test_report_details_string += (
                f"Module set `{module_set_name}`, initial condition "
                f"`{initial_condition_name}`: "
            )
            if test_passed[i_test][j_ic]:
                test_report_details_string += '*PASSED*\n'
            else:
                test_report_details_string += '*FAILED*\n'

    # Summarize the test results.
    test_report_summary_string = (
        'Summary of initial condition build test results from `ICtest.py`'
        f" for branch or commit or tag `{BRANCH_OR_COMMIT}`: "
    )
    if 'FAILED' in test_report_details_string:
        test_report_summary_string += '*FAILED*\n'
    else:
        test_report_summary_string += '*PASSED*\n'

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If a test failed, or loud mode is on, post report to Slack.
    if (slack_on_fail and 'FAILED' in test_report_details_string) or be_loud:
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
