#!/usr/bin/env python

"""Run the MAGE weekly dash tests.

This script runs the MAGE weekly dash tests.

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
DESCRIPTION = 'Run the MAGE weekly dash tests.'

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for weekly dash results
WEEKLY_DASH_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'weeklyDash')

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for weekly dash
WEEKLY_DASH_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'weekly_dash.lst')

# Prefix for weekly dash directory name
WEEKLY_DASH_DIRECTORY_PREFIX = 'weeklyDash_'

# Subdirectory of build directory containing compiled products to use in tests
BIN_DIR = 'bin'

# List of weekly dash test files to copy
WEEKLY_DASH_TEST_FILES = [
    'weeklyDashGo.xml',
]

# Paths to jinja2 template file for PBS script.
WEEKLY_DASH_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'weeklyDash-template.pbs'
)

# Name of rendered PBS script.
WEEKLY_DASH_PBS_SCRIPT = 'weeklyDash.pbs'


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

    # Set up for communication with Slack.
    if verbose:
        print('Creating Slack client.')
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the weekly dash tests.
    if verbose:
        print(f"Creating {WEEKLY_DASH_DIRECTORY}.")
    os.mkdir(WEEKLY_DASH_DIRECTORY)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    with open(WEEKLY_DASH_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Read the template for the PBS script used for the test data generation.
    with open(WEEKLY_DASH_PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    weekly_dash_pbs_template = Template(template_content)
    if debug:
        print(f"weekly_dash_pbs_template = {weekly_dash_pbs_template}")

    # ------------------------------------------------------------------------

    # Create the make command to build the code.
    make_cmd = 'make voltron_mpi.x'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # ------------------------------------------------------------------------

    # Create the list for submit results. Only set to True if all build and
    # qsub commands for a set are OK.
    submit_ok = [False]*len(module_list_files)
    if debug:
        print(f"submit_ok = {submit_ok}")

    # Create the list of job IDs.
    job_ids = [None]*len(module_list_files)
    if debug:
        print(f"jobs_ids = {job_ids}")

    # Run the weekly dash with each set of modules.
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing weekly dash with module set '
                  f"{module_list_file}")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_set_name = {module_set_name}.")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
        if debug:
            print(f"path = {path}")
        if verbose:
            print(f"Reading module list file {path}.")
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

        # Add the cmake option for the weekly dash build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=Release'
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Make a directory for this build, and go there.
        dir_name = f"{WEEKLY_DASH_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(WEEKLY_DASH_DIRECTORY, dir_name)
        if verbose:
            print(f"Creating and moving to build directory {build_directory}.")
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

        # Move into the bin directory to run the tests.
        os.chdir(BIN_DIR)

        # Generate the LFM grid file.
        if verbose:
            print('Creating LFM grid file.')
        cmd = 'genLFM.py -gid Q'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: Unable to create LFM grid file for module set '
                  f"{module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See testing log for output from genLFM.py.\n'
                  'Skipping remaining steps for module set'
                  f"{module_set_name}\n")
            continue

        # Generate the solar wind boundary condition file.
        if verbose:
            print('Creating solar wind initial conditions file.')
        cmd = 'cda2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: Unable to create solar wind boundary conditions file'
                  f" for module set {module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See testing log for output from cda2wind.py.\n'
                  'Skipping remaining steps for module set'
                  f"{module_set_name}\n")
            continue

        # Generate the RCM configuration file.
        if verbose:
            print('Creating RCM configuration file.')
        cmd = 'genRCM.py'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: Unable to create RCM configuration file'
                  f" for module set {module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See testing log for output from genRCM.py.\n'
                  'Skipping remaining steps for module set '
                  f"{module_set_name}\n")
            continue

        # Copy the PBS and XML files needed for the weekly dash job.
        if verbose:
            print('Copying PBS and XML files needed for weekly dash.')
        for filename in WEEKLY_DASH_TEST_FILES:
            from_file = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
            to_file = os.path.join('.', filename)
            shutil.copyfile(from_file, to_file)

        # Assemble data to fill in the PBS template.
        pbs_options = {}
        pbs_options['job_name'] = dir_name
        pbs_options['account'] = os.environ['DERECHO_TESTING_ACCOUNT']
        pbs_options['queue'] = os.environ['DERECHO_TESTING_QUEUE']
        pbs_options['job_priority'] = os.environ['DERECHO_TESTING_PRIORITY']
        pbs_options['walltime'] = '08:00:00'
        pbs_options['modules'] = module_names
        pbs_options['kaijuhome'] = KAIJUHOME
        pbs_options['tmpdir'] = os.environ['TMPDIR']
        pbs_options['slack_bot_token'] = os.environ['SLACK_BOT_TOKEN']
        pbs_options['mage_test_root'] = os.environ['MAGE_TEST_ROOT']
        pbs_options['mage_test_set_root'] = os.environ['MAGE_TEST_SET_ROOT']
        pbs_options['branch_or_commit'] = os.environ['BRANCH_OR_COMMIT']
        pbs_options['report_options'] = ''
        if debug:
            pbs_options['report_options'] += ' -d'
        pbs_options['report_options'] += ' -l'  # Always report.
        if slack_on_fail:
            pbs_options['report_options'] += ' -s'
        if is_test:
            pbs_options['report_options'] += ' -t'
        if verbose:
            pbs_options['report_options'] += ' -v'

        # Render the job template.
        pbs_content = weekly_dash_pbs_template.render(pbs_options)
        with open(WEEKLY_DASH_PBS_SCRIPT, 'w', encoding='utf-8') as f:
            f.write(pbs_content)

        # Submit the weekly dash job.
        if verbose:
            print('Submitting weekly dash model run.')
        cmd = f"qsub {WEEKLY_DASH_PBS_SCRIPT}"
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

        # Record the job ID.
        job_ids[i_module_set] = job_id

        # Record successful submission.
        submit_ok[i_module_set] = True

        # Save the job number in a file.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            f.write(f"{job_id}\n")

        # End of loop over module sets.

    if debug:
        print(f"submit_ok = {submit_ok}")
        print(f"job_ids = {job_ids}")

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Test results are in {WEEKLY_DASH_DIRECTORY}.\n"
    )
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if not submit_ok[i_module_set]:
            test_report_details_string += (
                f"Module set `{module_list_file}` submission: *FAILED*\n"
            )
            continue
        test_report_details_string += (
            f"Weekly dash for module set `{module_list_file}` submitted as "
            f"job {job_ids[i_module_set]}.\n"
        )


    # Summarize the test results.
    test_report_summary_string = (
        f"Weekly dash submission for `{os.environ['BRANCH_OR_COMMIT']}`: "
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
