#!/usr/bin/env python

"""Run MAGE Intel Inspector tests.

This script runs a series of tests of the MAGE software using Intel tools.

The Intel Inspector memory checks run in about an hour on two derecho nodes.

The Intel Inspector thread checks run in about 45 minutes on two derecho nodes.

The report script runs in about 40 minutes on one derecho node.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import sys
import shutil
import subprocess

# Import 3rd-party modules.
from jinja2 import Template

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE checks with Intel Inspector tools'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for Intel Inspector checks
INTEL_CHECKS_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'intelChecks')

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for Intel checks
INTEL_CHECKS_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'intelChecks.lst')

# Prefix for naming Intel Inspector checks directories
INTEL_CHECKS_DIRECTORY_PREFIX = 'intelChecks_'

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'

# Name of PBS account to use for testing jobs.
DERECHO_TESTING_ACCOUNT = os.environ['DERECHO_TESTING_ACCOUNT']

# Path to jinja2 template file for PBS script for the memory tests.
MEM_CHECK_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'intelCheckSubmitMem-template.pbs'
)

# Path to jinja2 template file for PBS script for the thread tests.
THREAD_CHECK_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'intelCheckSubmitThread-template.pbs'
)

# Path to jinja2 template file for PBS script for the reporting script.
REPORT_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'intelCheckSubmitReport-template.pbs'
)

# Name of PBS file for memory checks.
MEM_CHECK_PBS_FILENAME = 'intelCheckSubmitMem.pbs'

# Name of PBS file for thread checks.
THREAD_CHECK_PBS_FILENAME = 'intelCheckSubmitThread.pbs'

# Name of PBS file for report.
REPORT_PBS_FILENAME = 'intelCheckSubmitReport.pbs'


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
    # account = args.account
    debug = args.debug
    be_loud = args.loud
    is_test = args.test
    verbose = args.verbose

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the Intel Inspector tests.
    if verbose:
        print(f"Creating {INTEL_CHECKS_DIRECTORY}.")
    os.mkdir(INTEL_CHECKS_DIRECTORY)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of module sets to use for Intel checks.
    with open(INTEL_CHECKS_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Read the template for the PBS script used for the memory tests.
    with open(MEM_CHECK_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    mem_check_pbs_template = Template(template_content)
    if debug:
        print(f"mem_check_pbs_template = {mem_check_pbs_template}")

    # Read the template for the PBS script used for the thread tests.
    with open(THREAD_CHECK_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    thread_check_pbs_template = Template(template_content)
    if debug:
        print(f"thread_check_pbs_template = {thread_check_pbs_template}")

    # Read the template for the PBS script used for the reporting.
    with open(REPORT_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    report_pbs_template = Template(template_content)
    if debug:
        print(f"report_pbs_template = {report_pbs_template}")

    # ------------------------------------------------------------------------

    # Run the Intel Inspector checks with each set of modules.

    # Create the common make command for all module sets.
    make_cmd = 'make gamera_mpi voltron_mpi'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Run Intel checks with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing Intel Inspector checks with module set '
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

        # Add the additional flags needed for Intel Inspector checks.
        cmake_options += ' -DDISABLE_DEBUG_BOUNDS_CHECKS=ON'
        cmake_options += ' -DCMAKE_BUILD_TYPE=DEBUG'
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Make a directory for this test, and go there.
        dir_name = f"{INTEL_CHECKS_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(INTEL_CHECKS_DIRECTORY, dir_name)
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
                f"Skipping remaining steps for module set {module_set_name}.",
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
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: make for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'make.out')}"
                ' for output from make.\n'
                f"Skipping remaining steps for module set {module_set_name}.",
                file=sys.stderr
            )
            continue

        # Go to the bin directory for testing.
        os.chdir(BUILD_BIN_DIR)

        # Copy in the files used by the tests.
        test_files = [
            'tinyCase.xml',
            'bcwind.h5',
            'lfmD.h5',
            'rcmconfig.h5',
            'memSuppress.sup',
            'threadSuppress.sup',
        ]
        for filename in test_files:
            from_path = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
            to_path = os.path.join('.', filename)
            shutil.copyfile(from_path, to_path)

        # Assemble common data to fill in the PBS templates.
        pbs_options = {}
        pbs_options['account'] = DERECHO_TESTING_ACCOUNT
        pbs_options['queue'] = os.environ['DERECHO_TESTING_QUEUE']
        pbs_options['job_priority'] = os.environ['DERECHO_TESTING_PRIORITY']
        pbs_options['modules'] = module_names
        pbs_options['kaijuhome'] = KAIJUHOME
        pbs_options['tmpdir'] = os.environ['TMPDIR']
        pbs_options['slack_bot_token'] = os.environ['SLACK_BOT_TOKEN']
        pbs_options['mage_test_root'] = os.environ['MAGE_TEST_ROOT']
        pbs_options['derecho_testing_account'] = os.environ['DERECHO_TESTING_ACCOUNT']

        # Set options specific to the memory check, then render the template.
        pbs_options['job_name'] = 'mage_intelCheckSubmitMem'
        pbs_options['walltime'] = '12:00:00'
        pbs_content = mem_check_pbs_template.render(pbs_options)
        with open(MEM_CHECK_PBS_FILENAME, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Set options specific to the thread check, then render the template.
        pbs_options['job_name'] = 'mage_intelCheckSubmitThread'
        pbs_options['walltime'] = '12:00:00'
        pbs_content = thread_check_pbs_template.render(pbs_options)
        with open(THREAD_CHECK_PBS_FILENAME, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Set options specific to the report generator, then render the
        # template.
        pbs_options['job_name'] = 'mage_intelCheckReportSubmit'
        pbs_options['walltime'] = '02:00:00'
        pbs_content = report_pbs_template.render(pbs_options)
        with open(REPORT_PBS_FILENAME, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Initialize the job ID list.
        job_ids = []

        # Run the memory check job.
        cmd = f"qsub {MEM_CHECK_PBS_FILENAME}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(
                'ERROR: qsub failed.\n'
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See test log for output.\n'
                f"Skipping remaining steps for module set {module_set_name}.",
                file=sys.stderr
            )
            continue
        job_id = cproc.stdout.split('.')[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids.append(job_id)

        # Run the thread check job.
        cmd = f"qsub {THREAD_CHECK_PBS_FILENAME}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(
                'ERROR: qsub failed.\n'
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See test log for output.\n'
                f"Skipping remaining steps for module set {module_set_name}.",
                file=sys.stderr
            )
            continue
        job_id = cproc.stdout.split('.')[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids.append(job_id)

        # Run the report job when the other two jobs are complete.
        cmd = f"qsub -W depend=afterok:{':'.join(job_ids)} {REPORT_PBS_FILENAME}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(
                'ERROR: qsub failed.\n'
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See test log for output.\n'
                f"Skipping remaining steps for module set {module_set_name}.",
                file=sys.stderr
            )
            continue
        job_id = cproc.stdout.split('.')[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids.append(job_id)

        # Record the job IDs.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            for job_id in job_ids:
                f.write(f"{job_id}\n")

    # -----------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Test results are in {os.getcwd()}.\n"
    )
    pbs_files = [
        MEM_CHECK_PBS_FILENAME, THREAD_CHECK_PBS_FILENAME,
        REPORT_PBS_FILENAME
    ]
    for (pbs_file, job_id) in zip(pbs_files, job_ids):
        test_report_details_string += (
            f"{pbs_file} submitted as job {job_id}.\n"
        )

    # Summarize the test results
    test_report_summary_string = (
        'Intel Inspector tests submitted (`intelChecks.py`).'
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
