#!/usr/bin/env python

"""Run MAGE Intel Inspector tests.

This script runs a series of tests of the MAGE software using Intel tools.

The Intel Inspector memory checks run in about an hour on two derecho nodes.

The Intel Inspector thread checks run in about 45 minutes on two derecho nodes.

The report script runs in about 2 minutes on one derecho node.

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


# # Program constants

# Program description.
DESCRIPTION = "Script for MAGE checks with Intel Inspector tools"

# Home directory of kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ["MAGE_TEST_SET_ROOT"]

# Directory for Intel Inspector checks
INTEL_CHECKS_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, "intelChecks")

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, "testingScripts")

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     "mage_build_test_modules")

# Name of file containing names of modules lists to use for Intel checks
INTEL_CHECKS_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY,
                                      "intelChecks.lst")

# Path to jinja2 template file for PBS script for the memory tests.
MEM_CHECK_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "intelCheckSubmitMem-template.pbs"
)

# Path to jinja2 template file for PBS script for the thread tests.
THREAD_CHECK_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "intelCheckSubmitThread-template.pbs"
)

# Path to jinja2 template file for PBS script for the reporting script.
REPORT_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "intelCheckSubmitReport-template.pbs"
)

# Prefix for naming Intel Inspector checks directories
INTEL_CHECKS_DIRECTORY_PREFIX = "intelChecks_"

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = "bin"

# Data and configuration files used by the Intel Inspector tests.
TEST_INPUT_FILES = [
    "tinyCase.xml",
    "bcwind.h5",
    "memSuppress.sup",
    "threadSuppress.sup",
]

# Name of PBS file for memory checks.
MEM_CHECK_PBS_FILENAME = "intelCheckSubmitMem.pbs"

# Name of PBS file for thread checks.
THREAD_CHECK_PBS_FILENAME = "intelCheckSubmitThread.pbs"

# Name of PBS file for report generation.
REPORT_PBS_FILENAME = "intelCheckSubmitReport.pbs"

# Branch or commit (or tag) used for testing.
BRANCH_OR_COMMIT = os.environ["BRANCH_OR_COMMIT"]


def intelChecks(args: dict):
    """Begin main program.

    This is the main program code.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Local convenience variables.
    debug = args.get("debug", False)
    be_loud = args.get("loud", False)
    slack_on_fail = args.get("slack_on_fail", False)
    is_test = args.get("test", False)
    verbose = args.get("verbose", False)

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
    with open(INTEL_CHECKS_LIST_FILE, encoding="utf-8") as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Read the template for the PBS script used for the memory tests.
    with open(MEM_CHECK_PBS_TEMPLATE, encoding="utf-8") as f:
        template_content = f.read()
    mem_check_pbs_template = Template(template_content)
    if debug:
        print(f"mem_check_pbs_template = {mem_check_pbs_template}")

    # Read the template for the PBS script used for the thread tests.
    with open(THREAD_CHECK_PBS_TEMPLATE, encoding="utf-8") as f:
        template_content = f.read()
    thread_check_pbs_template = Template(template_content)
    if debug:
        print(f"thread_check_pbs_template = {thread_check_pbs_template}")

    # Read the template for the PBS script used for the report generation.
    with open(REPORT_PBS_TEMPLATE, encoding="utf-8") as f:
        template_content = f.read()
    report_pbs_template = Template(template_content)
    if debug:
        print(f"report_pbs_template = {report_pbs_template}")

    # ------------------------------------------------------------------------

    # Run the Intel Inspector checks with each set of modules.

    # Create the common make command for all module sets.
    make_cmd = "make gamera_mpi voltron_mpi"
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Create the list for submit results. Only set to True if all qsub
    # commands for a set are OK.
    submit_ok = [False]*len(module_list_files)
    if debug:
        print(f"submit_ok = {submit_ok}")

    # Create a list of lists for job IDs. There are 3 job IDs per set - one for
    # memory check, one for thread check, and one for the test report.
    job_ids = [[None, None, None]]*len(module_list_files)
    if debug:
        print(f"job_ids = {job_ids}")

    # Run Intel checks with each set of modules.
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if verbose:
            print("Performing Intel Inspector checks with module set "
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip(".lst")
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

        # Add the additional flags needed for Intel Inspector checks.
        cmake_options += " -DDISABLE_DEBUG_BOUNDS_CHECKS=ON"
        cmake_options += " -DCMAKE_BUILD_TYPE=DEBUG"
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # --------------------------------------------------------------------

        # Prepare the directory and commands for the test.

        # Make a directory for this test, and go there.
        dir_name = f"{INTEL_CHECKS_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(INTEL_CHECKS_DIRECTORY, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Run cmake to build the Makefile.
        if verbose:
            print("Running cmake to create Makefile for module set"
                  f" {module_set_name}.")
        cmd = (f"{module_cmd}; {cmake_environment} cmake {cmake_options}"
               f" {KAIJUHOME} >& cmake.out")
        if debug:
            print(f"cmd = {cmd}")
        try:
            # NOTE: stdout and stderr goes cmake.out.
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: cmake for module set {module_set_name} failed.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  f"See {os.path.join(build_directory, 'cmake.out')}"
                  " for output from cmake.\n"
                  "Skipping remaining steps for module set "
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        if debug:
            print(f"cproc = {cproc}")

        # Run the build.
        if verbose:
            print("Running make to build kaiju for module set"
                  f" {module_set_name}.")
        cmd = f"{module_cmd}; {make_cmd} >& make.out"
        if debug:
            print(f"cmd = {cmd}")
        try:
            # NOTE: stdout and stderr go into make.out.
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: make for module set {module_set_name} failed.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  f"See {os.path.join(build_directory, 'make.out')}"
                  " for output from make.\n"
                  "Skipping remaining steps for module set "
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        if debug:
            print(f"cproc = {cproc}")

        # --------------------------------------------------------------------

        # Copy or create input files for the tests.

        # Go to the bin directory for testing.
        os.chdir(BUILD_BIN_DIR)

        # Copy in the files used by the tests.
        if verbose:
            print("Copying files needed for Intel checks.")
        for filename in TEST_INPUT_FILES:
            from_path = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
            to_path = os.path.join(".", filename)
            shutil.copyfile(from_path, to_path)

        # Generate the LFM grid file.
        if verbose:
            print("Creating LFM grid file.")
        cmd = "genLFM.py -gid D"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print("ERROR: Unable to create LFM grid file for module set "
                  f"{module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  "See testing log for output from genLFM.py.\n"
                  "Skipping remaining steps for module set"
                  f"{module_set_name}\n")
            continue
        if debug:
            print(f"cproc = {cproc}")

        # Generate the RCM configuration file.
        if verbose:
            print("Creating RCM configuration file.")
        cmd = "genRCM.py"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print("ERROR: Unable to create RCM configuration file"
                  f" for module set {module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  "See testing log for output from genRCM.py.\n"
                  "Skipping remaining steps for module set "
                  f"{module_set_name}\n")
            continue
        if debug:
            print(f"cproc = {cproc}")

        # --------------------------------------------------------------------

        # Create the PBS scripts for the tests.

        # Assemble common data to fill in the PBS templates.
        pbs_options = {}
        pbs_options["account"] = os.environ["DERECHO_TESTING_ACCOUNT"]
        pbs_options["queue"] = os.environ["DERECHO_TESTING_QUEUE"]
        pbs_options["job_priority"] = os.environ["DERECHO_TESTING_PRIORITY"]
        pbs_options["modules"] = module_names
        pbs_options["kaijuhome"] = KAIJUHOME
        pbs_options["kaipyhome"] = os.environ["KAIPYHOME"]
        pbs_options["tmpdir"] = os.environ["TMPDIR"]
        pbs_options["slack_bot_token"] = os.environ["SLACK_BOT_TOKEN"]
        pbs_options["mage_test_root"] = os.environ["MAGE_TEST_ROOT"]
        pbs_options["branch_or_commit"] = BRANCH_OR_COMMIT

        # Set options specific to the memory check, then render the template.
        pbs_options["job_name"] = "mage_intelCheckSubmitMem"
        pbs_options["walltime"] = "04:00:00"
        pbs_content = mem_check_pbs_template.render(pbs_options)
        with open(MEM_CHECK_PBS_FILENAME, "w", encoding="utf-8") as f:
            f.write(pbs_content)

        # Set options specific to the thread check, then render the template.
        pbs_options["job_name"] = "mage_intelCheckSubmitThread"
        pbs_options["walltime"] = "04:00:00"
        pbs_content = thread_check_pbs_template.render(pbs_options)
        with open(THREAD_CHECK_PBS_FILENAME, "w", encoding="utf-8") as f:
            f.write(pbs_content)

        # Set options specific to the report generator, then render the
        # template.
        pbs_options["job_name"] = "mage_intelCheckReportSubmit"
        pbs_options["walltime"] = "02:00:00"
        pbs_options["report_options"] = ""
        if debug:
            pbs_options["report_options"] += " -d"
        if slack_on_fail:
            pbs_options["report_options"] += " -s"
        if is_test:
            pbs_options["report_options"] += " -t"
        if verbose:
            pbs_options["report_options"] += " -v"
        pbs_content = report_pbs_template.render(pbs_options)
        with open(REPORT_PBS_FILENAME, "w", encoding="utf-8") as f:
            f.write(pbs_content)

        # --------------------------------------------------------------------

        # Run the PBS scripts for the tests.

        # Run the memory check job.
        cmd = f"qsub {MEM_CHECK_PBS_FILENAME}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print("ERROR: qsub failed for memory check.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  "See test log for output.\n"
                  "Skipping remaining steps for module set "
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split(".")[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids[i_module_set][0] = job_id

        # Run the thread check job.
        cmd = f"qsub {THREAD_CHECK_PBS_FILENAME}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print("ERROR: qsub failed for thread check.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  "See test log for output.\n"
                  "Skipping remaining steps for module set "
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split(".")[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids[i_module_set][1] = job_id

        # Run the report job when the other two jobs are complete.
        cmd = (
            f"qsub -W depend=afterok:{':'.join(job_ids[i_module_set][:-1])} "
            f"{REPORT_PBS_FILENAME}"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print("ERROR: qsub failed for report generation.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  "See test log for output.\n"
                  "Skipping remaining steps for module set "
                  f"{module_set_name}.",
                  file=sys.stderr)
            continue
        job_id = cproc.stdout.split(".")[0]
        if debug:
            print(f"job_id = {job_id}")
        job_ids[i_module_set][2] = job_id

        # Record the job IDs for this module set in a file.
        with open("jobs.txt", "w", encoding="utf-8") as f:
            for job_id in job_ids[i_module_set]:
                f.write(f"{job_id}\n")

        # This module set worked.
        submit_ok[i_module_set] = True

    # End of loop over module sets.

    # -----------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ""
    test_report_details_string += (
        f"Test results are on `derecho` in `{INTEL_CHECKS_DIRECTORY}`.\n"
    )
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        test_report_details_string = (
            "Submit Intel Inspector tests for module set "
            f"`{module_list_file}`: "
        )
        if submit_ok[i_module_set]:
            test_report_details_string += "*PASSED*"
        else:
            test_report_details_string += "*FAILED*"

    # Summarize the test results
    test_report_summary_string = (
        f"Intel Inspector test submission for `{BRANCH_OR_COMMIT}`: "
    )
    if "FAILED" in test_report_details_string:
        test_report_summary_string += "*FAILED*"
    else:
        test_report_summary_string += "*PASSED*"

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If a test failed, or loud mode is on, post report to Slack.
    if (slack_on_fail and "FAILED" in test_report_summary_string) or be_loud:
        slack_client = common.slack_create_client()
        if debug:
            print(f"slack_client = {slack_client}")
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=is_test
        )
        if debug:
            print(f"slack_response_summary = {slack_response_summary}")
        thread_ts = slack_response_summary["ts"]
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_details_string, thread_ts=thread_ts,
            is_test=is_test
        )
        if debug:
            print(f"slack_response_summary = {slack_response_summary}")

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


def main():
    """Driver for command-line version of code."""
    # Set up the command-line parser.
    parser = common.create_command_line_parser(DESCRIPTION)

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")

    # Convert the arguments from Namespace to dict.
    args = vars(args)

    # Pass the command-line arguments to the main function as a dict.
    intelChecks(args)


if __name__ == "__main__":
    main()
