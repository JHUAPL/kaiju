#!/usr/bin/env python

"""Run MAGE build regression tests.

This script runs a series of builds of the MAGE software using sets of
modules listed in files under:

$KAIJUHOME/testingScripts/mage_build_test_modules

This script reads the file build_test.lst from this directory, and
uses the contents as a list of module list files to use for MAGE build
tests. Each build takes about 10 minutes.

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
DESCRIPTION = 'Script for MAGE build testing'

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for build tests
BUILD_TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'buildTest')

# Path to directory to use for building executable list
EXECUTABLE_LIST_BUILD_DIRECTORY = os.path.join(BUILD_TEST_DIRECTORY,
                                               'build_executable_list')

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Path to module list file to use when generating the list of executables
# Use a module set without MKL.
EXECUTABLE_LIST_MODULE_LIST = os.path.join(MODULE_LIST_DIRECTORY,
                                           'intel_mpich.lst')

# Path to file containing list of module sets to use for build tests
BUILD_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'build_test.lst')

# Prefix for naming build test directories
BUILD_TEST_DIRECTORY_PREFIX = 'buildTest_'

# Name of subdirectory of current build subdirectory containing binaries
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

    # Make a directory to hold all of the build tests.
    if verbose:
        print(f"Creating {BUILD_TEST_DIRECTORY}.")
    os.mkdir(BUILD_TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Do a preliminary cmake run to generate the list of executables.

    if verbose:
        print('Generating list of executables.')

    # Create and move to the preliminary build folder.
    os.mkdir(EXECUTABLE_LIST_BUILD_DIRECTORY)
    os.chdir(EXECUTABLE_LIST_BUILD_DIRECTORY)

    # Read the module list file for building the executable list.
    if verbose:
        print('Reading module list for executable list generation.')
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(EXECUTABLE_LIST_MODULE_LIST)
    )
    if debug:
        print(f"module_names = {module_names}")
        print(f"cmake_environment = {cmake_environment}")
        print(f"cmake_options = {cmake_options}")

    # Assemble the commands to load the listed modules.
    module_cmd = f"module --force purge; module load {' '.join(module_names)}"
    if debug:
        print(f"module_cmd = {module_cmd}")

    # Run cmake to build the Makefile.
    if verbose:
        print('Running cmake for executable list generation.')
    cmd = (
        f"{module_cmd}; {cmake_environment} cmake {cmake_options} {KAIJUHOME}"
        ' >& cmake.out'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: cmake for building executable list failed.\n'
            f"e.cmd = {e.cmd}"
            f"e.returncode = {e.returncode}\n"
            f"See {os.path.join(EXECUTABLE_LIST_BUILD_DIRECTORY, 'cmake.out')}"
            ' for output from cmake.\n'
            'Unable to generate executable list.',
            file=sys.stderr
        )
        raise

    # Run make to build the list of executable targets.
    if verbose:
        print('Running make for executable list generation.')
    pattern = r'\.x'
    cmd = f"{module_cmd}; make help | grep '{pattern}'"
    if debug:
        print(f"cmd = {cmd}")
    try:
        # NOTE: stdout and stderr combined into stdout
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        path = os.path.join(
            EXECUTABLE_LIST_BUILD_DIRECTORY, 'make_to_grep.out'
        )
        with open(path, 'w', encoding='utf-8') as f:
            f.write(e.stdout)
        print(
            'ERROR: make for building executable list failed.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f"See {path} for output from make piped to grep.\n"
            'Unable to generate executable list.',
            file=sys.stderr
        )
        raise
    executable_list = cproc.stdout.splitlines()
    if debug:
        print(f"executable_list = {executable_list}")

    # Remove the first four characters (dots and spaces).
    executable_list = [_[4:] for _ in executable_list]
    if debug:
        print(f"executable_list = {executable_list}")

    # Create the make command to build the executable list.
    make_cmd = f"make {' '.join(executable_list)}"
    if debug:
        print(f"make_cmd = {make_cmd}")

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    module_list_files, _, _ = common.read_build_module_list_file(
        BUILD_TEST_LIST_FILE)
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Initalize test results for all module sets to False (failed).
    test_passed = [False]*len(module_list_files)

    # Do a build with each set of modules.
    for (i_test, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing build test with module list file '
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_set_name = {module_set_name}.")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
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

        # Make a directory for this build, and go there.
        dir_name = f"{BUILD_TEST_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(BUILD_TEST_DIRECTORY, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Run cmake to build the Makefile.
        if verbose:
            print(
                'Running cmake to create Makefile for module set'
                f" {module_set_name}."
            )
        cmd = (
            f"{module_cmd}; {cmake_environment} cmake {cmake_options}"
            f" {KAIJUHOME} >& cmake.out"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            # NOTE: stdout and stderr goes to stdout (into log file)
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: cmake for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'cmake.out')}"
                ' for output from cmake.\n'
                f"Skipping remaining steps for module set {module_set_name}",
                file=sys.stderr
            )
            continue

        # Run the build.
        if verbose:
            print(
                'Running make to build kaiju for module set'
                f" {module_set_name}."
            )
        cmd = f"{module_cmd}; {make_cmd} >& make.out"
        if debug:
            print(f"cmd = {cmd}")
        try:
            # NOTE: stdout and stderr go into make.out.
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: make for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'make.out')}"
                ' for output from make.\n'
                f"Skipping remaining steps for module set {module_set_name}",
                file=sys.stderr
            )
            continue

        # Check for all executables.
        if verbose:
            print(
                f"Checking build results for module set {module_set_name}."
            )
        missing = []
        for executable in executable_list:
            path = os.path.join(build_directory, BUILD_BIN_DIR, executable)
            if not os.path.isfile(path):
                missing.append(executable)
        if len(missing) > 0:
            for executable in missing:
                print(f"ERROR: Did not build {executable}.")
        else:
            test_passed[i_test] = True

    # End of loop over module sets.

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Test results are on `derecho` in `{BUILD_TEST_DIRECTORY}`.\n"
    )
    for (i_test, module_list_file) in enumerate(module_list_files):
        module_set_name = module_list_file.rstrip('.lst')
        test_report_details_string += f"Module set `{module_set_name}`: "
        if test_passed[i_test]:
            test_report_details_string += '*PASSED*\n'
        else:
            test_report_details_string += '*FAILED*\n'

    # Summarize the test results.
    test_report_summary_string = (
        f"Build test results for `{BRANCH_OR_COMMIT}`: "
    )
    if 'FAILED' in test_report_details_string:
        test_report_summary_string += '*FAILED*'
    else:
        test_report_summary_string += '*PASSED*'

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If a test failed, or loud mode is on, post report to Slack.
    if (slack_on_fail and 'FAILED' in test_report_summary_string) or be_loud:
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
