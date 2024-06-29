#!/usr/bin/env python

"""Run MAGE Fortran unit tests.

This script runs a series of unit tests of the MAGE Fortran software. These
tests are run as PBS jobs on derecho. There will be one job which generates
the data for testing, then 1 or more dependent jobs that use the newly-
generated data for unit testing, then a job for the test report.

NOTE: If this script is run as part of a set of tests for run_mage_tests.sh,
this script must be listed *last*, since it makes changes to the kaiju source
code tree that are incompatible with the other tests.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
# import glob
import os
import shutil
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE Fortran unit testing'

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for unit tests
UNIT_TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'unitTest')

# Top-level directory for testing on derecho.
MAGE_TEST_ROOT = os.environ['MAGE_TEST_ROOT']

# Home directory for pFUnit compiled code
PFUNIT_HOME = os.path.join(
    MAGE_TEST_ROOT, 'pfunit', 'pFUnit-4.2.0', 'ifort-23-mpich-derecho'
)

# List of pFUnit directories to copy from PFUNIT_HOME into kaiju/external
PFUNIT_BINARY_DIRECTORIES = [
    'FARGPARSE-1.1',
    'GFTL-1.3',
    'GFTL_SHARED-1.2',
    'PFUNIT-4.2',
]

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Path to kaiju subdirectory for external code
KAIJU_EXTERNAL_DIRECTORY = os.path.join(KAIJUHOME, 'external')

# Prefix for naming unit test directories
UNIT_TEST_DIRECTORY_PREFIX = 'unitTest_'

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for unit tests
UNIT_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'unit_test.lst')

# Path to directory containing unit test files.
TESTS_DIRECTORY = os.path.join(KAIJUHOME, 'tests')

# PBS scripts for unit test jobs.
UNIT_TEST_PBS_SCRIPTS = [
    'genTestData.pbs',
    'runCaseTests.pbs',
    'runNonCaseTests1.pbs',
    # 'runNonCaseTests2.pbs',  # Hangs for 12 hours
    'unitTestReport.pbs',
]

# Input files for unit tests
UNIT_TEST_DATA_INPUT_DIRECTORY = os.path.join(
    os.environ['MAGE_TEST_ROOT'], 'unit_test_inputs'
)
UNIT_TEST_DATA_INPUT_FILES = [
    'bcwind.h5',
    'geo_mpi.xml',
    'lfmD.h5',
    'rcmconfig.h5',
]

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'

# Name of PBS account to use for testing jobs.
DERECHO_TESTING_ACCOUNT = os.environ['DERECHO_TESTING_ACCOUNT']

# Token string for access to Slack.
SLACK_BOT_TOKEN = os.environ['SLACK_BOT_TOKEN']

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
    account = args.account
    debug = args.debug
    be_loud = args.loud
    is_test = args.test
    verbose = args.verbose

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the unit tests.
    if verbose:
        print(f"Creating ${UNIT_TEST_DIRECTORY}.")
    os.mkdir(UNIT_TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Make a copy of the pFUnit code under kaiju/external.
    if verbose:
        print('Copying compiled pFUnit binaries.')
    for directory in PFUNIT_BINARY_DIRECTORIES:
        from_path = os.path.join(PFUNIT_HOME, directory)
        to_path = os.path.join(KAIJU_EXTERNAL_DIRECTORY, directory)
        shutil.copytree(from_path, to_path)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for unit tests.
    with open(UNIT_TEST_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Initalize job ID to None for all module set/PBS script combinations.
    job_ids = []
    for _ in module_list_files:
        job_ids.append([None]*len(UNIT_TEST_PBS_SCRIPTS))

    # Create the common make command for all module sets.
    make_cmd = 'make gamera_mpi voltron_mpi allTests'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Run the unit tests with each set of modules.
    for (i_set, module_list_file) in enumerate(module_list_files):
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

        # <HACK>
        # Extra argument needed for unit test build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=RELWITHDEBINFO'
        # </HACK>

        # Make a directory for this test, and go there.
        dir_name = f"{UNIT_TEST_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(UNIT_TEST_DIRECTORY, dir_name)
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
            # NOTE: stdout and stderr goes cmake.out.
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

        # Go to the bin directory for testing.
        os.chdir(BUILD_BIN_DIR)

        # Copy in the PBS scripts for unit testing.
        for filename in UNIT_TEST_PBS_SCRIPTS:
            from_path = os.path.join(TESTS_DIRECTORY, filename)
            to_path = os.path.join('.', filename)
            shutil.copyfile(from_path, to_path)

        # Copy in inputs for unit test data generation.
        for filename in UNIT_TEST_DATA_INPUT_FILES:
            from_path = os.path.join(
                UNIT_TEST_DATA_INPUT_DIRECTORY, filename
            )
            to_path = os.path.join('.', filename)
            shutil.copyfile(from_path, to_path)

        # Submit the jobs to create the test data and run the unit
        # tests. Note that the unit test jobs will only run if the
        # data generation job completes successfully.
        for (j_pbs, pbs_file) in enumerate(UNIT_TEST_PBS_SCRIPTS):
            job_id = None
            cmd = (
                f"qsub -A {DERECHO_TESTING_ACCOUNT} "
                f"-v MODULE_LIST='{' '.join(module_names)}',"
                f"KAIJUROOTDIR={KAIJUHOME},"
                f"MAGE_TEST_SET_ROOT={MAGE_TEST_SET_ROOT},"
                f"DERECHO_TESTING_ACCOUNT={DERECHO_TESTING_ACCOUNT},"
                f"SLACK_BOT_TOKEN={SLACK_BOT_TOKEN}"
            )
            # <HACK>
            # Assumes data generation job is first.
            if j_pbs > 0:
                cmd += f" -W depend=afterok:{job_ids[i_set][0]}"
            # </HACK>
            # <HACK>
            # Assumes report generation job is last.
            if pbs_file == 'unitTestReport.pbs':
                cmd += f" -W depend=afterok:{':'.join(job_ids[i_set][1:-1])}"
            # </HACK>
            cmd += f" {pbs_file}"
            if debug:
                print(f"cmd = {cmd}")
            try:
                cproc = subprocess.run(cmd, shell=True, check=True,
                                       text=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                path = os.path.join(
                    build_directory, BUILD_BIN_DIR, f"qsub_{j_pbs}.out"
                )
                with open(path, 'w', encoding='utf-8') as f:
                    f.write(e.stdout)
                print(
                    'ERROR: Job submission failed.\n'
                    f"e.cmd = {e.cmd}\n"
                    f"e.returncode = {e.returncode}\n"
                    f"See {path} for output from qsub.\n"
                    'Skipping remaining steps for module set '
                    f"{module_set_name}.",
                    file=sys.stderr
                )
                continue

            # Save the job ID.
            job_id = cproc.stdout.split('.')[0]
            if debug:
                print(f"job_id = {job_id}")
            job_ids[i_set][j_pbs] = job_id
        # End of loop over PBS scripts

        # Record the job IDs in a text file.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            for job_id in job_ids[i_set]:
                f.write(f"{job_id}\n")

    # End of loop over module sets

    # ------------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # ------------------------------------------------------------------------

    # NOTE: Assumes only 1 module set was used.

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Test results are in {os.getcwd()}.\n"
    )
    test_report_details_string += (
        'Fortran unit test PBS job script `genTestData.pbs` submitted as job '
        f"{job_ids[0][0]}.\n"
    )
    test_report_details_string += (
        'Fortran unit test PBS job script `runCaseTests.pbs` submitted as job '
        f"{job_ids[0][1]}.\n"
    )
    test_report_details_string += (
        'Fortran unit test PBS job script `runNonCaseTests1.pbs` submitted as'
        f" job {job_ids[0][2]}.\n"
    )
    test_report_details_string += (
        'Fortran unit test PBS job script `runNonCaseTests2.pbs` skipped'
        ' since it currently hangs on `derecho`.\n'
    )
    test_report_details_string += (
        'Fortran unit test report PBS job script `unitTestReport.pbs` submitted as'
        f" job {job_ids[0][3]}.\n"
    )

    # Summarize the test results
    test_report_summary_string = (
        'Fortran unit tests submitted by `unitTest.py`'
        f" for branch or commit or tag {BRANCH_OR_COMMIT}\n"
    )

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If loud mode is on, post report to Slack.
    if be_loud:
        test_report_summary_string += 'Details in thread for this messsage.\n'
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=is_test
        )
        if slack_response_summary['ok']:
            thread_ts = slack_response_summary['ts']
            slack_response_details = common.slack_send_message(
                slack_client, test_report_details_string, thread_ts=thread_ts,
                is_test=is_test
            )
            if 'ok' not in slack_response_details:
                print('*ERROR* Unable to post test details to Slack.')
        else:
            print('*ERROR* Unable to post test summary to Slack.')

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    main()
