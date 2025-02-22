#!/usr/bin/env python

"""Run MAGE Fortran unit tests.

This script runs a series of unit tests of the MAGE Fortran software. These
tests are run as PBS jobs on derecho. There will be one job which generates
the data for testing, then 3 dependent jobs that use the newly-generated data
for unit testing, then a job for the test report.

There are 5 PBS job scripts used per module set. Each is generated from a
jinja2 template.

1. genTestData.pbs - Data generation. Runs in about 4-5 minutes on 5 derecho
   nodes. Output in PBS job file genTestData.o*, and cmiD_deep_8_genRes.out.

2. runCaseTests.pbs - Runs in about 35 minutes on 1 derecho node. Only runs if
   genTestData.pbs completes successfully. Output in PBS log file
   runCaseTests.o*, caseTests.out, and caseMpiTests.out.

3. runNonCaseTests1.pbs - Runs in about 2 minutes on 1 derecho node. Only runs
   if genTestData.pbs completes successfully. Output in PBS log file
   runNonCaseTests1.o*, gamTests.out, mixTests.out, voltTests.out,
   baseMpiTests.out, gamMpiTests.out. shgrTests.out.
   NOTE: As of 2024-08-22, voltTests.out will contain errors like this when
   run on the development branch:

   ...
   [testebsquish.pf:151]
   Squish Fake Projection Latitude value is wrong. Check Squish Processing and Output.
   AssertEqual failure:
         Expected: <147591.2572518899>
           Actual: <143412.6753716097>
       Difference: <-4178.581880280253> (greater than tolerance of .1000000000000000E-06)
   ...

4. runNonCaseTests2.pbs - Runs in about XX minutes on 2 derecho nodes. Only
   runs if genTestData.pbs completes successfully. Output in PBS log file
   runNonCaseTests2.o*, and voltMpiTests.out.

5. unitTestReport.pbs - Report generation. Runs in about XX minutes on 1
   derecho node. Only runs if jobs 2-4 complete successfully. Output in PBS
   log file unitTestReport.o*, and unitTestReport.out.

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
import os
import shutil
import subprocess
import sys

# Import 3rd-party modules.
from jinja2 import Template

# Import project modules.
import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE Fortran unit testing'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

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

# List of pFUnit directories to copy from PFUNIT_HOME into
# kaiju_private/external
PFUNIT_BINARY_DIRECTORIES = [
    'FARGPARSE-1.1',
    'GFTL-1.3',
    'GFTL_SHARED-1.2',
    'PFUNIT-4.2',
]

# Path to kaiju subdirectory for external code
KAIJU_EXTERNAL_DIRECTORY = os.path.join(KAIJUHOME, 'external')

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for unit tests
UNIT_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'unit_test.lst')

# Path to directory containing the unit test scripts
UNIT_TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'tests')

# Paths to jinja2 template files for PBS scripts.
DATA_GENERATION_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_SCRIPTS_DIRECTORY, 'genTestData-template.pbs'
)
RUN_CASE_TESTS_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_SCRIPTS_DIRECTORY, 'runCaseTests-template.pbs'
)
RUN_NON_CASE_TESTS_1_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_SCRIPTS_DIRECTORY, 'runNonCaseTests1-template.pbs'
)
RUN_NON_CASE_TESTS_2_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_SCRIPTS_DIRECTORY, 'runNonCaseTests2-template.pbs'
)
UNIT_TEST_REPORT_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_SCRIPTS_DIRECTORY, 'unitTestReport-template.pbs'
)

# Prefix for naming unit test directories
UNIT_TEST_DIRECTORY_PREFIX = 'unitTest_'

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'

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

# Names of PBS scripts to create from templates.
DATA_GENERATION_PBS_SCRIPT = 'genTestData.pbs'
RUN_CASE_TESTS_PBS_SCRIPT = 'runCaseTests.pbs'
RUN_NON_CASE_TESTS_1_PBS_SCRIPT = 'runNonCaseTests1.pbs'
RUN_NON_CASE_TESTS_2_PBS_SCRIPT = 'runNonCaseTests2.pbs'
UNIT_TEST_REPORT_PBS_SCRIPT = 'unitTestReport.pbs'

# Branch or commit (or tag) used for testing.
BRANCH_OR_COMMIT = os.environ['BRANCH_OR_COMMIT']

# Name of file to hold job list.
JOB_LIST_FILE = 'jobs.txt'


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
    slack_on_fail = args.slack_on_fail
    is_test = args.test
    verbose = args.verbose

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the Fortran unit tests.
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
        if debug:
            print(f"Copying {from_path} to {to_path}.")
        shutil.copytree(from_path, to_path)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.
    if verbose:
        print(f"Reading module set list from {UNIT_TEST_LIST_FILE}.")

    # Read the list of  module sets to use for unit tests.
    with open(UNIT_TEST_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    if verbose:
        print('Reading templates for PBS scripts.')

    # Read the template for the PBS script used for the test data generation.
    with open(DATA_GENERATION_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    data_generation_pbs_template = Template(template_content)
    if debug:
        print(f"data_generation_pbs_template = {data_generation_pbs_template}")

    # Read the template for the PBS script used for the case tests.
    with open(RUN_CASE_TESTS_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    run_case_tests_pbs_template = Template(template_content)
    if debug:
        print(f"run_case_tests_pbs_template = {run_case_tests_pbs_template}")

    # Read the template for the PBS script used for the 1st non-case tests.
    with open(RUN_NON_CASE_TESTS_1_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    run_non_case_tests_1_pbs_template = Template(template_content)
    if debug:
        print('run_non_case_tests_1_pbs_template = '
              f"{run_non_case_tests_1_pbs_template}")

    # Read the template for the PBS script used for the 2nd non-case tests.
    with open(RUN_NON_CASE_TESTS_2_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    run_non_case_tests_2_pbs_template = Template(template_content)
    if debug:
        print('run_non_case_tests_2_pbs_template = '
              f"{run_non_case_tests_2_pbs_template}")

    # Read the template for the PBS script used for report generation.
    with open(UNIT_TEST_REPORT_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    unit_test_report_pbs_template = Template(template_content)
    if debug:
        print('unit_test_report_pbs_template = '
              f"{unit_test_report_pbs_template}")

    # ------------------------------------------------------------------------

    # Create the common make command for all module sets.
    make_cmd = 'make gamera_mpi voltron_mpi allTests'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Create the list for submit results. Only set to True if all qsub commands
    # for a set are OK.
    submit_ok = [False]*len(module_list_files)
    if debug:
        print(f"submit_ok = {submit_ok}")

    # Create a list of lists for job IDs. There are 5 job IDs per set - one for
    # data generration, case tests, non-case tests 1, non-case tests 2, and the
    # test report.
    job_ids = [[None, None, None, None, None]]*len(module_list_files)
    if debug:
        print(f"job_ids = {job_ids}")

    # Run the unit tests with each set of modules.
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing unit tests tests with module set '
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_set_name = {module_set_name}.")

        # --------------------------------------------------------------------

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

        # <HACK>
        # Extra argument needed for unit test build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=RELWITHDEBINFO'
        if debug:
            print(f"cmake_options = {cmake_options}")
        # </HACK>

        # Assemble the command to load the listed modules.
        module_cmd = (
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

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
            cproc = subprocess.run(cmd, shell=True, check=True)
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

        # --------------------------------------------------------------------

        # Assemble common data to fill in the PBS templates.
        pbs_options = {}
        pbs_options['account'] = os.environ['DERECHO_TESTING_ACCOUNT']
        pbs_options['queue'] = os.environ['DERECHO_TESTING_QUEUE']
        pbs_options['job_priority'] = os.environ['DERECHO_TESTING_PRIORITY']
        pbs_options['modules'] = module_names
        pbs_options['kaijuhome'] = KAIJUHOME
        pbs_options['branch_or_commit'] = BRANCH_OR_COMMIT

        # Go to the bin directory for testing.
        os.chdir(BUILD_BIN_DIR)

        # --------------------------------------------------------------------

        # Copy in inputs for unit test data generation.
        for filename in UNIT_TEST_DATA_INPUT_FILES:
            from_path = os.path.join(
                UNIT_TEST_DATA_INPUT_DIRECTORY, filename
            )
            to_path = os.path.join('.', filename)
            if debug:
                print(f"Copying {from_path} to {to_path}.")
            shutil.copyfile(from_path, to_path)

        # Set options specific to the data generation job, then render the
        # template.
        pbs_options['job_name'] = 'genTestData'
        pbs_options['walltime'] = '00:30:00'
        pbs_content = data_generation_pbs_template.render(pbs_options)
        if verbose:
            print(f"Creating {DATA_GENERATION_PBS_SCRIPT}.")
        with open(DATA_GENERATION_PBS_SCRIPT, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Run the data generation job.
        cmd = f"qsub {DATA_GENERATION_PBS_SCRIPT}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: qsub failed.\n'
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See test log for output.\n'
                  'Skipping remaining steps for module set '
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split('.')[0]
        job_ids[i_module_set][0] = job_id
        if debug:
            print(f"job_id = {job_id}")
            print(f"job_ids = {job_ids}")

        # --------------------------------------------------------------------

        # Set options specific to the case tests job, then render the
        # template.
        pbs_options['job_name'] = 'runCaseTests'
        pbs_options['walltime'] = '00:40:00'
        pbs_content = run_case_tests_pbs_template.render(pbs_options)
        if verbose:
            print(f"Creating {RUN_CASE_TESTS_PBS_SCRIPT}.")
        with open(RUN_CASE_TESTS_PBS_SCRIPT, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Run the case tests job if data was generated.
        cmd = (
            f"qsub -W depend=afterok:{job_ids[i_module_set][0]} "
            f"{RUN_CASE_TESTS_PBS_SCRIPT}"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: qsub failed.\n'
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See test log for output.\n'
                  'Skipping remaining steps for module set '
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split('.')[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids[i_module_set][1] = job_id

        # --------------------------------------------------------------------

        # Set options specific to the 1st non-case tests job, then render the
        # template.
        pbs_options['job_name'] = 'runNonCaseTests1'
        pbs_options['walltime'] = '00:05:00'
        if verbose:
            print(f"Creating {RUN_NON_CASE_TESTS_1_PBS_SCRIPT}.")
        pbs_content = run_non_case_tests_1_pbs_template.render(pbs_options)
        with open(RUN_NON_CASE_TESTS_1_PBS_SCRIPT, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Run the 1st non-case tests job if data was generated.
        cmd = (
            f"qsub -W depend=afterok:{job_ids[i_module_set][0]} "
            f"{RUN_NON_CASE_TESTS_1_PBS_SCRIPT}"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: qsub failed.\n'
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See test log for output.\n'
                  'Skipping remaining steps for module set '
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split('.')[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids[i_module_set][2] = job_id

        # --------------------------------------------------------------------

        # Set options specific to the 2nd non-case tests job, then render the
        # template.
        pbs_options['job_name'] = 'runNonCaseTests2'
        pbs_options['walltime'] = '12:00:00'
        pbs_content = run_non_case_tests_2_pbs_template.render(pbs_options)
        with open(RUN_NON_CASE_TESTS_2_PBS_SCRIPT, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Run the 2nd non-case tests job if data was generated.
        cmd = (
            f"qsub -W depend=afterok:{job_ids[i_module_set][0]} "
            f"{RUN_NON_CASE_TESTS_2_PBS_SCRIPT}"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: qsub failed.\n'
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See test log for output.\n'
                  'Skipping remaining steps for module set '
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split('.')[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids[i_module_set][3] = job_id

        # --------------------------------------------------------------------

        # Set options specific to the report generation job, then render the
        # template.
        pbs_options['job_name'] = 'unitTestReport'
        pbs_options['walltime'] = '00:10:00'
        pbs_options['slack_bot_token'] = os.environ['SLACK_BOT_TOKEN']
        pbs_options['mage_test_root'] = os.environ['MAGE_TEST_ROOT']
        pbs_options['mage_test_set_root'] = os.environ['MAGE_TEST_SET_ROOT']
        pbs_options['report_options'] = ''
        if debug:
            pbs_options['report_options'] += ' -d'
        if be_loud:
            pbs_options['report_options'] += ' -l'
        if slack_on_fail:
            pbs_options['report_options'] += ' -s'
        if is_test:
            pbs_options['report_options'] += ' -t'
        if verbose:
            pbs_options['report_options'] += ' -v'
        pbs_content = unit_test_report_pbs_template.render(pbs_options)
        with open(UNIT_TEST_REPORT_PBS_SCRIPT, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Run the report generation job if all others ran OK.
        cmd = (
            f"qsub -W depend=afterok:{':'.join(job_ids[i_module_set][1:4])} "
            f"{UNIT_TEST_REPORT_PBS_SCRIPT}"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: qsub failed.\n'
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See test log for output.\n'
                  'Skipping remaining steps for module set '
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split('.')[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids[i_module_set][4] = job_id

        # --------------------------------------------------------------------

        # Record the job IDs for this module set in a file.
        if verbose:
            print(f"Saving job IDs for module set {module_set_name} "
                  f"in {JOB_LIST_FILE}.")
        with open(JOB_LIST_FILE, 'w', encoding='utf-8') as f:
            for job_id in job_ids[i_module_set]:
                f.write(f"{job_id}\n")

        # This module set worked.
        submit_ok[i_module_set] = True

    # End of loop over module sets
    if debug:
        print(f"submit_ok = {submit_ok}")
        print(f"job_ids = {job_ids}")

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Test results are on `derecho` in `{UNIT_TEST_DIRECTORY}`.\n"
    )
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if not submit_ok[i_module_set]:
            test_report_details_string += (
                f"Module set `{module_list_file}`: *FAILED*"
            )
            continue
        test_report_details_string += (
            f"`{DATA_GENERATION_PBS_SCRIPT}` for module set "
            f"`{module_list_file}` submitted as PBS job "
            f"{job_ids[i_module_set][0]}.\n"
        )
        test_report_details_string += (
            f"`{RUN_CASE_TESTS_PBS_SCRIPT}` for module set "
            f"`{module_list_file}` submitted as PBS job "
            f"{job_ids[i_module_set][1]}.\n"
        )
        test_report_details_string += (
            f"`{RUN_NON_CASE_TESTS_1_PBS_SCRIPT}` for module set "
            f"`{module_list_file}` submitted as PBS job "
            f"{job_ids[i_module_set][2]}.\n"
        )
        test_report_details_string += (
            f"`{RUN_NON_CASE_TESTS_2_PBS_SCRIPT}` for module set "
            f"`{module_list_file}` submitted as PBS job "
            f"{job_ids[i_module_set][3]}.\n"
        )
        test_report_details_string += (
            f"`{UNIT_TEST_REPORT_PBS_SCRIPT}` for module set "
            f"`{module_list_file}` submitted as PBS job "
            f"{job_ids[i_module_set][4]}.\n"
        )

    # Summarize the test results
    test_report_summary_string = (
        f"Unit test submission for `{os.environ['BRANCH_OR_COMMIT']}`: "
    )
    if 'FAILED' in test_report_details_string:
        test_report_summary_string += '*FAILED*'
    else:
        test_report_summary_string += '*PASSED*'

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
