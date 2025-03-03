#!/usr/bin/env python

"""Run the MAGE reproducibility check.

Run the MAGE reproducibility check. The reproducibility check makes two
duplicate runs, and then numerically compares the results.

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
DESCRIPTION = "Run a MAGE reproducibility check."

# Home directory of kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Default module set file
DEFAULT_MODULE_SET_FILE = os.path.join(
    KAIJUHOME, "testingScripts", "mage_build_test_modules", "intel_mpich.lst"
)

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ["MAGE_TEST_SET_ROOT"]

# Directory for reproducibility check results
REPRODUCIBILITY_CHECK_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT,
                                               "reproducibility_check")

# Prefix for reproducibility check directory name
REPRODUCIBILITY_CHECK_DIRECTORY_PREFIX = "reproducibility_check_"

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, "testingScripts")

# List of weekly dash test files to copy
WEEKLY_DASH_TEST_FILES = [
    "weeklyDashGo.xml",
]

# Path to jinja2 template file for PBS script for build job.
BUILD_MAGE_PBS_TEMPLATE_FILE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "build_mage-template.pbs"
)

# Name of rendered PBS script to build MAGE.
BUILD_MAGE_PBS_SCRIPT = "build_mage.pbs"

# Path to jinja2 template file for PBS script for run jobs.
RUN_MAGE_PBS_TEMPLATE_FILE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "run_mage-template.pbs"
)

# Name of rendered PBS script for MAGE runs.
RUN1_MAGE_PBS_SCRIPT = "run1_mage.pbs"
RUN2_MAGE_PBS_SCRIPT = "run2_mage.pbs"

# Path to jinja2 template file for PBS script for comparison.
MAGE_REPRODUCIBILITY_CHECK_PBS_TEMPLATE_FILE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "mage_reproducibility_check-template.pbs"
)

# Name of rendered PBS script for MAGE run comparison.
MAGE_REPRODUCIBILITY_CHECK_PBS_SCRIPT = "mage_reproducibility_check.pbs"


def mage_reproducibility_check(args: dict):
    """Perform a MAGE reproducibility check.

    Perform a MAGE reproducibility check. A reproducibility check is composed
    of the following steps:

        1. Build the MAGE software.
        2. Make a single run using the standard weekly dash settings and
           inputs.
        3. Perform a second run which is a duplicate of the run from step 2.
        4. Perform a detailed numerical comparison of the results from the
           two runs.

    The reproducibility check is split into multiple PBS jobs - one job per
    step in the procedure shown above. Job #1 runs first, and if it finishes
    successfully, jobs #2 and #3 are run in parallel (in separate directories).
    Assuming #2 and #3 finish successfully, job #4 is run.

    Note that this script does not make an assumption about the version of
    the code to test, or the modules used to build it. Typically, for weekly
    reproducibility checks, the module set will be the set which uses the
    Intel compiler but does *not* use MKL, and the nominal version of the code
    used in the reproducibililty check is the latest commit on the development
    branch.

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
    print(f"module_names = {module_names}")

    # Extract the name of the list.
    module_set_name = os.path.split(module_set_file)[-1].rstrip(".lst")

    # Add the cmake option for the weekly dash build.
    cmake_options += " -DCMAKE_BUILD_TYPE=Release"

    # ------------------------------------------------------------------------

    # Make a directory for this test, and go there.
    dir_name = f"{REPRODUCIBILITY_CHECK_DIRECTORY_PREFIX}{module_set_name}"
    build_directory = os.path.join(REPRODUCIBILITY_CHECK_DIRECTORY, dir_name)
    if verbose:
        print(f"Creating and moving to build directory {build_directory}.")
    os.makedirs(build_directory)
    os.chdir(build_directory)

    # ------------------------------------------------------------------------

    # Assemble and run the job to build the software.
    if verbose:
        print("Creating PBS job to build MAGE software.")

    # Read the template for the PBS script.
    with open(BUILD_MAGE_PBS_TEMPLATE_FILE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Assemble commands needed in the PBS script.
    cmake_cmd = f"{cmake_environment} cmake {cmake_options} {KAIJUHOME}"
    make_cmd = "make voltron_mpi.x"

    # Assemble data to fill in the PBS template.
    pbs_options = {}
    pbs_options["job_name"] = f"build-{dir_name}"
    pbs_options["account"] = os.environ["DERECHO_TESTING_ACCOUNT"]
    pbs_options["queue"] = os.environ["DERECHO_TESTING_QUEUE"]
    pbs_options["job_priority"] = os.environ["DERECHO_TESTING_PRIORITY"]
    pbs_options["walltime"] = "00:20:00"
    pbs_options["modules"] = module_names
    pbs_options["cmake_cmd"] = cmake_cmd
    pbs_options["make_cmd"] = make_cmd

    # Render the job template.
    pbs_content = pbs_template.render(pbs_options)
    with open(BUILD_MAGE_PBS_SCRIPT, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Submit the job.
    cmd = f"qsub {BUILD_MAGE_PBS_SCRIPT}"
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    jobid_build = cproc.stdout.split(".")[0]

    # ------------------------------------------------------------------------

    # Assemble and run the job to execute the first MAGE run.
    if verbose:
        print("Creating PBS job for first MAGE run.")

    # Read the template for the PBS script.
    with open(RUN_MAGE_PBS_TEMPLATE_FILE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Assemble commands needed in the PBS script.
    genLFM_cmd = "genLFM.py -gid Q"
    cda2wind_cmd = (
        "cda2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00"
    )
    genRCM_cmd = "genRCM.py"
    mpiexec_cmd = f"mpiexec {KAIJUHOME}/scripts/preproc/pinCpuCores.sh"
    voltron_cmd = "../bin/voltron_mpi.x weeklyDashGo.xml"

    # Make a directory for this run and go there.
    run_dir = "run1"
    os.mkdir(run_dir)
    os.chdir(run_dir)

    # Copy the input files for the run.
    for filename in WEEKLY_DASH_TEST_FILES:
        from_file = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
        to_file = os.path.join(".", filename)
        shutil.copyfile(from_file, to_file)

    # Assemble data to fill in the PBS template.
    pbs_options = {}
    pbs_options["job_name"] = f"run1-{dir_name}"
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
    with open(RUN1_MAGE_PBS_SCRIPT, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Submit the job to run after the build job.
    cmd = f"qsub -W depend=afterok:{jobid_build} {RUN1_MAGE_PBS_SCRIPT}"
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    jobid_run1 = cproc.stdout.split(".")[0]

    # ------------------------------------------------------------------------

    # Assemble and run the job to execute the second MAGE run.
    if verbose:
        print("Creating PBS job for second MAGE run.")

    # Make a directory for this run and go there.
    os.chdir(build_directory)
    run_dir = "run2"
    os.mkdir(run_dir)
    os.chdir(run_dir)

    # Copy the input files for the run.
    for filename in WEEKLY_DASH_TEST_FILES:
        from_file = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
        to_file = os.path.join(".", filename)
        shutil.copyfile(from_file, to_file)

    # Assemble data to fill in the PBS template - just change the name.
    pbs_options["job_name"] = f"run2-{dir_name}"

    # Render the job template.
    pbs_content = pbs_template.render(pbs_options)
    with open(RUN2_MAGE_PBS_SCRIPT, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Submit the job to run after the build job.
    cmd = f"qsub -W depend=afterok:{jobid_build} {RUN2_MAGE_PBS_SCRIPT}"
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    jobid_run2 = cproc.stdout.split(".")[0]

    # ------------------------------------------------------------------------

    # Assemble and run the job to execute the comparison.
    if verbose:
        print("Creating PBS job for run comparison.")

    # Move to the directory containing both run directories.
    os.chdir(build_directory)

    # Read the template for the PBS script.
    with open(MAGE_REPRODUCIBILITY_CHECK_PBS_TEMPLATE_FILE, "r",
              encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Assemble data to fill in the PBS template.
    pbs_options = {}
    pbs_options["job_name"] = f"comparison-{dir_name}"
    pbs_options["account"] = os.environ["DERECHO_TESTING_ACCOUNT"]
    pbs_options["queue"] = os.environ["DERECHO_TESTING_QUEUE"]
    pbs_options["job_priority"] = os.environ["DERECHO_TESTING_PRIORITY"]
    pbs_options["walltime"] = "02:00:00"
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
    pbs_options["xml1"] = os.path.join(build_directory, "run1",
                                       "weeklyDashGo.xml")
    pbs_options["xml2"] = os.path.join(build_directory, "run2",
                                       "weeklyDashGo.xml")

    # Render the job template.
    pbs_content = pbs_template.render(pbs_options)
    with open(MAGE_REPRODUCIBILITY_CHECK_PBS_SCRIPT, "w",
              encoding="utf-8") as f:
        f.write(pbs_content)

    # Submit the job to run after both run jobs complete OK.
    cmd = (f"qsub -W depend=afterok:{':'.join([jobid_run1, jobid_run2])} "
           f"{MAGE_REPRODUCIBILITY_CHECK_PBS_SCRIPT}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    jobid_comparison = cproc.stdout.split(".")[0]

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
        f"Reproducibillity check for `{os.environ['BRANCH_OR_COMMIT']}` "
        f"using module set `{module_set_file}` submitted as jobs "
        f"`{jobid_build}`, `{jobid_run1}`, `{jobid_run2}`, "
        f"`{jobid_comparison}`."
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

#     # Add additional arguments specific to this script.
#     parser.add_argument(
#         "--module_set_file", "-f", default=DEFAULT_MODULE_SET_FILE,
#         help=(
#             "Path to text file containing set of modules to build with "
#             "(default: %(default)s)"
#         )
#     )

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")

    # Convert the arguments from Namespace to dict.
    args = vars(args)

    # Pass the command-line arguments to the main function as a dict.
    mage_reproducibility_check(args)


if __name__ == "__main__":
    main()
