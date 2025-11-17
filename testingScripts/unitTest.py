#!/usr/bin/env python

"""Run MAGE Fortran unit tests.

This script runs a series of unit tests of the MAGE Fortran software. These
tests are run as PBS jobs on derecho. There will be one job which builds the
code, one job which generates the data for testing, 3 jobs that use the
newly-generated data for unit testing, and a job for the test report.

There are 6 PBS job scripts used per module set. Each is generated from a
jinja2 template.

1. unitTest-build.pbs - Build the kaiju code for unit testing. Runs in about
   21 minutes on a single derecho node. Output in PBS job file
   unitTest-build.o*, cmake.out, and make.out.

2. genTestData.pbs - Data generation. Runs after the build job in about 4-5
   minutes on 5 derecho nodes. Output in PBS job file genTestData.o*, and
   cmiD_deep_8_genRes.out.

3. runCaseTests.pbs - Runs after the data generation job in about 35 minutes
   on 1 derecho node. Output in PBS log file runCaseTests.o*, caseTests.out,
   and caseMpiTests.out.

4. runNonCaseTests1.pbs - Runs after the data generation job in about 6
   minutes on 1 derecho node. Output in PBS log file runNonCaseTests1.o*,
   gamTests.out, mixTests.out, voltTests.out, baseMpiTests.out,
   gamMpiTests.out. shgrTests.out.

5. runNonCaseTests2.pbs - Runs after the data generation job in about XX
   minutes on 2 derecho nodes. Output in PBS log file runNonCaseTests2.o*, and
   voltMpiTests.out.

6. unitTestReport.pbs - Report generation. Runs after all test jobs in about
   XX minutes on 1 derecho node. Output in PBS log file unitTestReport.o*, and
   unitTestReport.out.

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
import subprocess
import sys

# Import 3rd-party modules.
from jinja2 import Template

# Import project modules.
import common


# Program constants

# Program description.
DESCRIPTION = "Script for MAGE Fortran unit testing"

# Home directory of kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ["MAGE_TEST_SET_ROOT"]

# Directory for unit tests
UNIT_TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, "unitTest")

# Top-level directory for testing on derecho.
MAGE_TEST_ROOT = os.environ["MAGE_TEST_ROOT"]

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, "testingScripts")

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     "mage_build_test_modules")

# Name of file containing names of modules lists to use for unit tests
UNIT_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, "unit_test.lst")

# Path to directory containing the unit test CODE.
UNIT_TEST_CODE_DIRECTORY = os.path.join(KAIJUHOME, "tests")

# Paths to jinja2 template files for PBS scripts.
BUILD_PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, "unitTest-build-template.pbs"
)
DATA_GENERATION_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_CODE_DIRECTORY, "genTestData-template.pbs"
)
RUN_CASE_TESTS_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_CODE_DIRECTORY, "runCaseTests-template.pbs"
)
RUN_NON_CASE_TESTS_1_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_CODE_DIRECTORY, "runNonCaseTests1-template.pbs"
)
RUN_NON_CASE_TESTS_2_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_CODE_DIRECTORY, "runNonCaseTests2-template.pbs"
)
UNIT_TEST_REPORT_PBS_TEMPLATE = os.path.join(
    UNIT_TEST_CODE_DIRECTORY, "unitTestReport-template.pbs"
)

# Prefix for naming unit test directories
UNIT_TEST_DIRECTORY_PREFIX = "unitTest_"

# Names of PBS scripts to create from templates.
BUILD_PBS_SCRIPT = "unitTest-build.pbs"
DATA_GENERATION_PBS_SCRIPT = "genTestData.pbs"
RUN_CASE_TESTS_PBS_SCRIPT = "runCaseTests.pbs"
RUN_NON_CASE_TESTS_1_PBS_SCRIPT = "runNonCaseTests1.pbs"
RUN_NON_CASE_TESTS_2_PBS_SCRIPT = "runNonCaseTests2.pbs"
UNIT_TEST_REPORT_PBS_SCRIPT = "unitTestReport.pbs"

# Name of file to hold job list.
JOB_LIST_FILE = "jobs.txt"


def create_build_pbs_script(module_list_file: str):
    """Create the PBS script to build the code.

    Create the PBS script to build the code.

    Parameters
    ----------
    module_list_file : str
        Path to module list file for the build.

    Returns
    -------
    pbs_script_name : str
        Name of PBS script file.

    Raises
    ------
    None
    """
    # Read this module list file, extracting cmake environment and
    # options, if any.
    path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(path)
    )

    # <HACK>
    # Extra argument needed for unit test build.
    cmake_options += " -DCMAKE_BUILD_TYPE=RELWITHDEBINFO"
    # </HACK>

    # Read the template for the PBS script.
    with open(BUILD_PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Create the options dictionary for the template.
    options = {
        "job_name": "unitTest-build",
        "account": os.environ["DERECHO_TESTING_ACCOUNT"],
        "queue": os.environ["DERECHO_TESTING_QUEUE"],
        "job_priority": os.environ["DERECHO_TESTING_PRIORITY"],
        "walltime": "00:30:00",
        "modules": module_names,
        "mage_test_root": MAGE_TEST_ROOT,
        "kaijuhome": KAIJUHOME,
        "cmake_cmd": f"{cmake_environment} cmake {cmake_options} {KAIJUHOME}",
        "make_cmd": "make gamera_mpi voltron_mpi allTests",
    }

    # Render the template.
    pbs_content = pbs_template.render(options)

    # Write the rendered file.
    pbs_script_name = BUILD_PBS_SCRIPT
    with open(pbs_script_name, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Return the name of the script.
    return pbs_script_name


def create_data_generation_pbs_script(module_list_file: str):
    """Create the PBS script to generate the test data.

    Create the PBS script to generate the test data.

    Parameters
    ----------
    module_list_file : str
        Path to module list file for the build.

    Returns
    -------
    pbs_script_name : str
        Name of PBS script file.

    Raises
    ------
    None
    """
    # Read this module list file, extracting cmake environment and
    # options, if any.
    path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(path)
    )

    # Read the template for the PBS script.
    with open(DATA_GENERATION_PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Create the options dictionary for the template.
    options = {
        "job_name": "unitTest-genTestData",
        "account": os.environ["DERECHO_TESTING_ACCOUNT"],
        "queue": os.environ["DERECHO_TESTING_QUEUE"],
        "job_priority": os.environ["DERECHO_TESTING_PRIORITY"],
        "walltime": "00:30:00",
        "modules": module_names,
        "kaijuhome": KAIJUHOME,
        "mage_test_root": MAGE_TEST_ROOT,
    }

    # Render the template.
    pbs_content = pbs_template.render(options)

    # Write the rendered file.
    pbs_script_name = DATA_GENERATION_PBS_SCRIPT
    with open(pbs_script_name, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Return the name of the script.
    return pbs_script_name


def create_case_tests_pbs_script(module_list_file: str):
    """Create the PBS script to run the case tests.

    Create the PBS script to run the case tests.

    Parameters
    ----------
    module_list_file : str
        Path to module list file for the build.

    Returns
    -------
    pbs_script_name : str
        Name of PBS script file.

    Raises
    ------
    None
    """
    # Read this module list file, extracting cmake environment and
    # options, if any.
    path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(path)
    )

    # Read the template for the PBS script.
    with open(RUN_CASE_TESTS_PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Create the options dictionary for the template.
    options = {
        "job_name": "unitTest-caseTests",
        "account": os.environ["DERECHO_TESTING_ACCOUNT"],
        "queue": os.environ["DERECHO_TESTING_QUEUE"],
        "job_priority": os.environ["DERECHO_TESTING_PRIORITY"],
        "walltime": "00:40:00",
        "modules": module_names,
        "kaijuhome": KAIJUHOME,
        "mage_test_root": MAGE_TEST_ROOT,
    }

    # Render the template.
    pbs_content = pbs_template.render(options)

    # Write the rendered file.
    pbs_script_name = RUN_CASE_TESTS_PBS_SCRIPT
    with open(pbs_script_name, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Return the name of the script.
    return pbs_script_name


def create_noncase_tests1_pbs_script(module_list_file: str):
    """Create the PBS script to run the first noncase tests.

    Create the PBS script to run the first noncase tests.

    Parameters
    ----------
    module_list_file : str
        Path to module list file for the build.

    Returns
    -------
    pbs_script_name : str
        Name of PBS script file.

    Raises
    ------
    None
    """
    # Read this module list file, extracting cmake environment and
    # options, if any.
    path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(path)
    )

    # Read the template for the PBS script.
    with open(RUN_NON_CASE_TESTS_1_PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Create the options dictionary for the template.
    options = {
        "job_name": "unitTest-noncaseTests1",
        "account": os.environ["DERECHO_TESTING_ACCOUNT"],
        "queue": os.environ["DERECHO_TESTING_QUEUE"],
        "job_priority": os.environ["DERECHO_TESTING_PRIORITY"],
        "walltime": "01:00:00",
        "modules": module_names,
        "kaijuhome": KAIJUHOME,
    }

    # Render the template.
    pbs_content = pbs_template.render(options)

    # Write the rendered file.
    pbs_script_name = RUN_NON_CASE_TESTS_1_PBS_SCRIPT
    with open(pbs_script_name, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Return the name of the script.
    return pbs_script_name


def create_noncase_tests2_pbs_script(module_list_file: str):
    """Create the PBS script to run the second noncase tests.

    Create the PBS script to run the second noncase tests.

    Parameters
    ----------
    module_list_file : str
        Path to module list file for the build.

    Returns
    -------
    pbs_script_name : str
        Name of PBS script file.

    Raises
    ------
    None
    """
    # Read this module list file, extracting cmake environment and
    # options, if any.
    path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(path)
    )

    # Read the template for the PBS script.
    with open(RUN_NON_CASE_TESTS_2_PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Create the options dictionary for the template.
    options = {
        "job_name": "unitTest-noncaseTests2",
        "account": os.environ["DERECHO_TESTING_ACCOUNT"],
        "queue": os.environ["DERECHO_TESTING_QUEUE"],
        "job_priority": os.environ["DERECHO_TESTING_PRIORITY"],
        "walltime": "12:00:00",
        "modules": module_names,
        "kaijuhome": KAIJUHOME,
    }

    # Render the template.
    pbs_content = pbs_template.render(options)

    # Write the rendered file.
    pbs_script_name = RUN_NON_CASE_TESTS_2_PBS_SCRIPT
    with open(pbs_script_name, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Return the name of the script.
    return pbs_script_name


def create_report_pbs_script(module_list_file: str, args: dict):
    """Create the PBS script to run the test report.

    Create the PBS script to run the test report.

    Parameters
    ----------
    module_list_file : str
        Path to module list file for the build.
    args : dict
        Command-line options and values.

    Returns
    -------
    pbs_script_name : str
        Name of PBS script file.

    Raises
    ------
    None
    """
    # Read this module list file, extracting cmake environment and
    # options, if any.
    path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(path)
    )

    # Read the template for the PBS script.
    with open(UNIT_TEST_REPORT_PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    pbs_template = Template(template_content)

    # Assemble the report options string.
    report_options = ""
    if args["debug"]:
        report_options += " -d"
    if args["loud"]:
        report_options += " -l"
    if args["slack_on_fail"]:
        report_options += " -s"
    if args["test"]:
        report_options += " -t"
    if args["verbose"]:
        report_options += " -v"

    # Create the options dictionary for the template.
    options = {
        "job_name": "unitTest-report",
        "account": os.environ["DERECHO_TESTING_ACCOUNT"],
        "queue": os.environ["DERECHO_TESTING_QUEUE"],
        "job_priority": os.environ["DERECHO_TESTING_PRIORITY"],
        "walltime": "00:10:00",
        "modules": module_names,
        "conda_environment": os.environ["CONDA_ENVIRONMENT"],
        "kaijuhome": KAIJUHOME,
        "mage_test_set_root": MAGE_TEST_SET_ROOT,
        "slack_bot_token": os.environ["SLACK_BOT_TOKEN"],
        "branch_or_commit": os.environ["BRANCH_OR_COMMIT"],
        "report_options": report_options,
    }

    # Render the template.
    pbs_content = pbs_template.render(options)

    # Write the rendered file.
    pbs_script_name = UNIT_TEST_REPORT_PBS_SCRIPT
    with open(pbs_script_name, "w", encoding="utf-8") as f:
        f.write(pbs_content)

    # Return the name of the script.
    return pbs_script_name


def qsub(pbs_script: str, qsub_options: str = ""):
    """Submit a PBS script.

    Submit a PBS script.

    Parameters
    ----------
    pbs_script : str
        Path to script to submit.
    qsub_options : str
        Options for qsub command.

    Returns
    -------
    job_id : int
        PBS job ID

    Raises
    ------
    subprocess.CalledProcessError
        If an error occures when running the qsub command.
    """
    # Assemble the command.
    cmd = f"qsub {qsub_options} {pbs_script}"
    print(f"{cmd=}")

    # Submit the job.
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"qsub failed with return code {e.returncode} for script "
              f"{pbs_script}.\n", file=sys.stderr)
        print(e.stderr)
        raise
    job_id = cproc.stdout.split(".")[0]

    # Return the job ID.
    return job_id


def unitTest(args: dict = None):
    """Run the unit tests.

    Run the unit tests.

    Parameters
    ----------
    args : dict
        Dictionary of program options.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Local convenience variables.
    debug = args["debug"]
    be_loud = args["loud"]
    slack_on_fail = args["slack_on_fail"]
    is_test = args["test"]
    verbose = args["verbose"]

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the Fortran unit tests.
    if verbose:
        print(f"Creating {UNIT_TEST_DIRECTORY}.")
    os.mkdir(UNIT_TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.
    if verbose:
        print(f"Reading module set list from {UNIT_TEST_LIST_FILE}.")

    # Read the list of  module sets to use for unit tests.
    with open(UNIT_TEST_LIST_FILE, encoding="utf-8") as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Create a list for submission status.
    submit_ok = [False]*len(module_list_files)

    # Run the unit tests with each set of modules.
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if verbose:
            print(f"Running unit tests with module set {module_list_file}.")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip(".lst")
        if debug:
            print(f"{module_set_name=}")

        # Make a directory for this test, and go there.
        dir_name = f"{UNIT_TEST_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(UNIT_TEST_DIRECTORY, dir_name)
        if debug:
            print(f"{build_directory=}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Create the PBS script for the build job.
        build_pbs_script = create_build_pbs_script(module_list_file)

        # Create the PBS script for data generation.
        data_generation_pbs_script = create_data_generation_pbs_script(
            module_list_file
        )

        # Create the PBS script for the case tests.
        run_case_tests_pbs_script = create_case_tests_pbs_script(
            module_list_file
        )

        # Create the PBS script for the first non-case tests.
        run_noncase_tests1_pbs_script = create_noncase_tests1_pbs_script(
            module_list_file
        )

        # Create the PBS script for the second non-case tests.
        run_noncase_tests2_pbs_script = create_noncase_tests2_pbs_script(
            module_list_file
        )

        # Create the PBS script for the unit test report.
        run_report_pbs_script = create_report_pbs_script(
            module_list_file, args
        )

        # Submit the build job.
        qsub_options = ""
        build_job_id = qsub(
            build_pbs_script, qsub_options
        )

        # Submit the data generation job, to run after the build job is done.
        qsub_options = f"-W depend=afterany:{build_job_id}"
        data_generation_job_id = qsub(
            data_generation_pbs_script, qsub_options
        )

        # Submit the case tests job, after the data generation job finishes.
        qsub_options = f"-W depend=afterany:{data_generation_job_id}"
        case_tests_job_id = qsub(
            run_case_tests_pbs_script, qsub_options
        )

        # Submit the first noncase tests job, after the data generation job
        # finishes.
        qsub_options = f"-W depend=afterany:{data_generation_job_id}"
        noncase_tests1_job_id = qsub(
            run_noncase_tests1_pbs_script, qsub_options
        )

        # Submit the second noncase tests job, after the data generation job
        # finishes.
        qsub_options = f"-W depend=afterany:{data_generation_job_id}"
        noncase_tests2_job_id = qsub(
            run_noncase_tests2_pbs_script, qsub_options
        )

        # Submit the test report job, which will runs after all test jobs.
        qsub_options = (
            "-W depend=afterany"
            + f":{case_tests_job_id}"
            + f":{noncase_tests1_job_id}"
            + f":{noncase_tests2_job_id}"
        )
        report_job_id = qsub(
            run_report_pbs_script, qsub_options
        )

        # Record the job IDs for this module set in a file.
        with open(JOB_LIST_FILE, "w", encoding="utf-8") as f:
            f.write(f"{build_job_id}\n")
            f.write(f"{data_generation_job_id}\n")
            f.write(f"{case_tests_job_id}\n")
            f.write(f"{noncase_tests1_job_id}\n")
            f.write(f"{noncase_tests2_job_id}\n")
            f.write(f"{report_job_id}\n")

        # Record the submit status for this module set.
        submit_ok[i_module_set] = True

    # End of loop over module sets

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ""
    test_report_details_string += (
        f"Unit test results are on `derecho` in `{UNIT_TEST_DIRECTORY}`.\n"
    )
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if not submit_ok[i_module_set]:
            test_report_details_string += (
                f"Unit test submit for module set `{module_list_file}`: "
                "*FAILED*"
            )

    # Summarize the test results
    test_report_summary_string = (
        f"Unit test submission for `{os.environ["BRANCH_OR_COMMIT"]}`: "
    )
    if False in submit_ok:
        test_report_summary_string += "*FAILED*"
    else:
        test_report_summary_string += "*PASSED*"

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If a test failed, or loud mode is on, post report to Slack.
    if (slack_on_fail and "FAILED" in test_report_details_string) or be_loud:
        slack_client = common.slack_create_client()
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=is_test
        )
        thread_ts = slack_response_summary["ts"]
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_details_string, thread_ts=thread_ts,
            is_test=is_test
        )

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


def main():
    """Main program code.

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
        print(f"{args=}")

    # ------------------------------------------------------------------------

    # Call the main program logic. Note that the Namespace object (args)
    # returned from the option parser is converted to a dict using vars().
    unitTest(vars(args))


if __name__ == "__main__":
    """Begin main program."""
    main()
