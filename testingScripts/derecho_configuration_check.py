#!/usr/bin/env python

"""Run the derecho configuration change check.

Run the derecho configuration change check. This script runs a PBS job which
builds and runs a MAGE model and numerically compares it to a reference run.

Authors
-------
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
DESCRIPTION = "Run the derecho configuration change check."

# Home directory of kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Default module set file
DEFAULT_MODULE_SET_FILE = os.path.join(
    KAIJUHOME, "testingScripts", "mage_build_test_modules", "intel_mpich.lst"
)

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ["MAGE_TEST_SET_ROOT"]

# Directory for derecho configuration check results
DERECHO_CONFIGURATION_CHECK = os.path.join(MAGE_TEST_SET_ROOT,
                                           "derecho_configuration_check")

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, "testingScripts")

# Path to jinja2 template file for PBS script.
PBS_TEMPLATE_FILE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "derecho_configuration_check-template.pbs"
)

# Prefix for weekly dash directory name
DERECHO_CONFIGURATION_CHECK_PREFIX = "derecho_configuration_check_"

# Name of rendered PBS script.
DERECHO_CONFIGURATION_CHECK_PBS_SCRIPT = "derecho_configuration_check.pbs"

# List of weekly dash test files to copy
WEEKLY_DASH_TEST_FILES = [
    "weeklyDashGo.xml",
]


def derecho_configuration_check(args: dict):
    """Run the derecho configuration change check.

    Run the derecho configuration change check.

    Parameters
    ----------
    args : dict
        Dictionary of command-line and other options.

    Returns
    -------
    None

    Raises
    ------
    subprocess.CalledProcessError
        If an exception occurs in subprocess.run()
    """
    # Local convenience variables.
    debug = args.get("debug", False)
    loud = args.get("loud", False)
    slack_on_fail = args.get("slack_on_fail", False)
    test = args.get("test", False)
    verbose = args.get("verbose", False)
    module_set_file = args.get("module_set_file", DEFAULT_MODULE_SET_FILE)

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Read the module list file, extracting the cmake environment and cmake
    # options, if any.
    if verbose:
        print(f"Reading module list file {module_set_file}.")
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(module_set_file)
    )

    # Extract the name of the list.
    module_set_name = os.path.split(module_set_file)[-1].rstrip(".lst")

    # Add the cmake option for the weekly dash build.
    cmake_options += " -DCMAKE_BUILD_TYPE=Release"

    # ------------------------------------------------------------------------

    # Read the template for the PBS script used for the weekly dash.
    with open(PBS_TEMPLATE_FILE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # ------------------------------------------------------------------------

    # Make a directory for this build, and go there.
    dir_name = f"{DERECHO_CONFIGURATION_CHECK_PREFIX}{module_set_name}"
    build_directory = os.path.join(DERECHO_CONFIGURATION_CHECK, dir_name)
    if verbose:
        print(f"Creating and moving to build directory {build_directory}.")
    os.makedirs(build_directory)
    os.chdir(build_directory)

    # ------------------------------------------------------------------------

    # Assemble commands needed in the PBS script.

    # Create the cmake command to build the Makefile.
    cmake_cmd = f"{cmake_environment} cmake {cmake_options} {KAIJUHOME}"

    # Create the make command to build the code.
    make_cmd = "make voltron_mpi.x"

    # Create the command to generate the LFM grid.
    genLFM_cmd = "genLFM.py -gid Q"

    # Create the command to generate the solar wind boundary condition file.
    cda2wind_cmd = (
        "cda2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00"
    )

    # Create the command to generate the RCM configuration.
    genRCM_cmd = "genRCM.py"

    # Create the command for launching an MPI program.
    mpiexec_cmd = f"mpiexec {KAIJUHOME}/scripts/preproc/pinCpuCores.sh"

    # Create the command to run the model.
    voltron_cmd = "./voltron_mpi.x weeklyDashGo.xml"

    # ------------------------------------------------------------------------

    # Copy the input files for the weekly dash job to the build directory.
    if verbose:
        print("Copying XML files needed for weekly dash.")
    for filename in WEEKLY_DASH_TEST_FILES:
        from_file = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
        to_file = os.path.join(".", filename)
        shutil.copyfile(from_file, to_file)

    # ------------------------------------------------------------------------

    # Assemble data to fill in the PBS template.
    pbs_options = {}
    pbs_options["job_name"] = dir_name
    pbs_options["account"] = os.environ["DERECHO_TESTING_ACCOUNT"]
    pbs_options["queue"] = os.environ["DERECHO_TESTING_QUEUE"]
    pbs_options["job_priority"] = os.environ["DERECHO_TESTING_PRIORITY"]
    pbs_options["walltime"] = "08:00:00"
    pbs_options["modules"] = module_names
    pbs_options["condarc"] = os.environ["CONDARC"]
    pbs_options["conda_envs_path"] = os.environ["CONDA_ENVS_PATH"]
    pbs_options["conda_environment"] = os.environ["CONDA_ENVIRONMENT"]
    pbs_options["mage_test_root"] = os.environ["MAGE_TEST_ROOT"]
    pbs_options["mage_test_set_root"] = os.environ["MAGE_TEST_SET_ROOT"]
    pbs_options["kaijuhome"] = KAIJUHOME
    pbs_options["kaipy_private_root"] = os.environ["KAIPY_PRIVATE_ROOT"]
    pbs_options["tmpdir"] = os.environ["TMPDIR"]
    pbs_options["slack_bot_token"] = os.environ["SLACK_BOT_TOKEN"]
    pbs_options["branch_or_commit"] = os.environ["BRANCH_OR_COMMIT"]
    pbs_options["cmake_cmd"] = cmake_cmd
    pbs_options["make_cmd"] = make_cmd
    pbs_options["genLFM_cmd"] = genLFM_cmd
    pbs_options["cda2wind_cmd"] = cda2wind_cmd
    pbs_options["genRCM_cmd"] = genRCM_cmd
    pbs_options["mpiexec_cmd"] = mpiexec_cmd
    pbs_options["voltron_cmd"] = voltron_cmd

    # Render the job template.
    pbs_content = pbs_template.render(pbs_options)
    with open(DERECHO_CONFIGURATION_CHECK_PBS_SCRIPT, "w",
              encoding="utf-8") as f:
        f.write(pbs_content)

    # Submit the weekly dash job.
    if verbose:
        print("Submitting weekly dash job "
              f"{DERECHO_CONFIGURATION_CHECK_PBS_SCRIPT}.")
    cmd = f"qsub {DERECHO_CONFIGURATION_CHECK_PBS_SCRIPT}"
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    jobid = cproc.stdout.split(".")[0]

    # Save the job ID.
    with open("jobid.txt", "w", encoding="utf-8") as f:
        f.write(f"{jobid}\n")

    # ------------------------------------------------------------------------

    # Summarize the test.

    # Set up for communication with Slack.
    if verbose:
        print("Creating Slack client.")
    slack_client = common.slack_create_client()

    # Detail the test results
    test_details = ""
    test_details += (
        f"Test results are in `{build_directory}`.\n"
    )

    # Summarize the test results.
    test_summary = (
        f"`derecho` configuration check for `{os.environ['BRANCH_OR_COMMIT']}`"
        f" using module set `{module_set_file}` submitted as job {jobid}."
    )

    # Print the test results summary and details.
    if verbose:
        print(test_summary)
        print(test_details)

    # If loud mode is on, post report to Slack. The initial message is the
    # test summary, and the thread contains the test details.
    if loud:
        slack_response = common.slack_send_message(
            slack_client, test_summary, is_test=test
        )
        thread_ts = slack_response["ts"]
        slack_response = common.slack_send_message(
            slack_client, test_details, thread_ts=thread_ts,
            is_test=test
        )

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


def main():
    """Driver for command-line version of code."""
    # Set up the command-line parser.
    parser = common.create_command_line_parser(DESCRIPTION)

    # Add additional arguments specific to this script.
    parser.add_argument(
        "--module_set_file", "-f", default=DEFAULT_MODULE_SET_FILE,
        help=(
            "Path to text file containing set of modules to build with "
            "(default: %(default)s)"
        )
    )

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")

    # Convert the arguments from Namespace to dict.
    args = vars(args)

    # Pass the command-line arguments to the main function as a dict.
    derecho_configuration_check(args)


if __name__ == "__main__":
    main()
