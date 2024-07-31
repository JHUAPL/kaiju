#!/usr/bin/env python

"""Run MAGE python unit tests.

This script runs a series of unit tests of the MAGE python software.

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
DESCRIPTION = 'Script for MAGE python unit testing'


# Paths under this test session directory.

# Path to directory for this set of python unit tests.
PYUNIT_TEST_DIRECTORY = os.path.join(os.environ['MAGE_TEST_SET_ROOT'],
                                     'pyunitTest')

# Name of PBS file to create from the jinja2 template for the tests.
PBS_FILE = 'pyunitTest.pbs'


# Paths under the kaiju installation to test.

# Top of installation tree for kaiju installation to test.
KAIJUHOME = os.environ['KAIJUHOME']

# Path to testing script directory.
KAIJU_TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to pytests directory.
KAIJU_PYTESTS_DIRECTORY = os.path.join(KAIJUHOME, 'pytests')

# Path to jinja2 template file for PBS script for weekly dash runs.
PBS_TEMPLATE = os.path.join(
    KAIJU_TEST_SCRIPTS_DIRECTORY, 'pyunitTest-template.pbs'
)


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

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Make a directory to hold the python unit tests, and go there.
    if verbose:
        print(f"Creating {PYUNIT_TEST_DIRECTORY}.")
    os.mkdir(PYUNIT_TEST_DIRECTORY)
    os.chdir(PYUNIT_TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Create the PBS script for this test session from the template.

    # Read the template for the PBS script used for the tests.
    with open(PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    pbs_template = Template(template_content)
    if debug:
        print(f"pbs_template = {pbs_template}")

    # Set the values for the template fields.
    pbs_options = {}
    pbs_options['job_name'] = 'pyunitTest'
    pbs_options['account'] = os.environ['DERECHO_TESTING_ACCOUNT']
    pbs_options['queue'] = 'main'
    pbs_options['job_priority'] = 'economy'
    pbs_options['walltime'] = '01:00:00'
    pbs_options['select'] = '1'
    pbs_options['ncpus'] = '128'
    pbs_options['mpiprocs'] = '1'
    pbs_options['ompthreads'] = '128'
    pbs_options['mage_test_root'] = os.environ['MAGE_TEST_ROOT']
    pbs_options['cdf_setup_script'] = (
        f"{os.environ['MAGE_TEST_ROOT']}/local/cdf/3.9.0/bin/definitions.B"
    )
    pbs_options['condarc'] = os.environ['CONDARC']
    pbs_options['conda_envs_path'] = os.environ['CONDA_ENVS_PATH']
    pbs_options['conda_environment'] = 'kaiju-3.8-testing'
    pbs_options['kaijuhome'] = os.environ['KAIJUHOME']
    pbs_options['tmpdir'] = os.environ['TMPDIR']
    pbs_options['slack_bot_token'] = os.environ['SLACK_BOT_TOKEN']
    pbs_options['pytest_output_file'] = 'kaiju-pyunit.txt'

    # Render the template and write it to a file.
    pbs_content = pbs_template.render(pbs_options)
    with open(PBS_FILE, 'w', encoding='utf-8') as f:
        f.write(pbs_content)

    # ------------------------------------------------------------------------

    # Copy the python unit test files from the source tree.
    for filename in ['test_satcomp_cdasws.py']:
        from_path = os.path.join(KAIJU_PYTESTS_DIRECTORY, filename)
        to_path = os.path.join('.', filename)
        shutil.copyfile(from_path, to_path)

    # ------------------------------------------------------------------------

    # Run the python unit tests as a PBS job.

    # Submit the unit test script for python.
    cmd = f"qsub {PBS_FILE}"
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    readString = cproc.stdout.rstrip()
    if debug:
        print(f"readString = {readString}")
    job_name_1 = readString.split('.')[0]
    if debug:
        print(f"job_name_1 = {job_name_1}")
    if verbose:
        print(
            f"Python unit test PBS script {PBS_FILE} submitted as "
            f"job {job_name_1}."
        )

    # -------------------------------------------------------------------------

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
    test_report_details_string += (
        'Python unit test PBS job script `pyunit.pbs` submitted as job '
        f"{job_name_1}.\n"
    )

    # Summarize the test results
    test_summary_message = (
        'Python unit tests submitted by `pyunitTest.py`'
        f" for branch or commit or tag {os.environ['BRANCH_OR_COMMIT']}: "
    )

    # Print the test results summary and details.
    print(test_summary_message)
    print(test_report_details_string)

    # If loud mode is on, post report to Slack.
    if be_loud:
        test_summary_message += 'Details in thread for this messsage.\n'
        slack_response_summary = common.slack_send_message(
            slack_client, test_summary_message, is_test=is_test
        )
        if slack_response_summary['ok']:
            thread_ts = slack_response_summary['ts']
            slack_response_details = common.slack_send_message(
                slack_client, test_report_details_string, thread_ts=thread_ts,
                is_test=is_test
            )
        else:
            print('*ERROR* Unable to post test result summary to Slack.')

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    main()
