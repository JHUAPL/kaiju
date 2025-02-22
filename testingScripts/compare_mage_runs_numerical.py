#!/usr/bin/env python


"""Numerically compare a pair of MAGE model runs.

This script generates a report which numerically compares the results of two
MAGE model runs. A run is specified using the path to the XML file which
defined the run. Comparisons are made using the h5diff tool.

NOTE: This script does not check any attributes in the datasets in the HDF5
files. This approach is taken to avoid false test failures when variable
attributes like kzcsMHD and kzcsTOT are compared.

NOTE: This script only compares groups with names of the form "Step#NNN" where
NNN is the step number. This script does *not* compare other top-level groups.

Authors
-------
Eric Winter

"""


# Import standard modules.
import datetime
import glob
import os
import platform
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
import common
from kaipy import kaiH5


# Program constants

# Program description.
DESCRIPTION = "Compare MAGE model runs numerically."

# Strings to represent test pass and fail.
TEST_PASS = "PASS"
TEST_FAIL = "FAIL"


def compare_GAMERA_results(runxml1: str, runxml2: str, verbose: bool = False):
    """Numerically compare the GAMERA output files from two runs.

    Numerically compare the GAMERA output files from two runs.

    Parameters
    ----------
    runxm1 : str
        Path to XML file describing 1st run.
    runxm2 : str
        Path to XML file describing 2nd run.
    verbose : bool
        Set to True to print verbose information during comparison.

    Returns
    -------
    TEST_PASS or TEST_FAIL : str
        Description of result of comparison.

    Raises
    ------
    None
    """
    # Determine the directories containing the sets of results.
    dir1 = os.path.split(runxml1)[0]
    dir2 = os.path.split(runxml2)[0]

    # Generate a sorted list of output files for the 1st run.
    pattern1 = os.path.join(dir1, "*.gam.h5")
    files1 = glob.glob(pattern1)
    files = [os.path.split(f)[1] for f in files1]
    files.sort()

    # Compare each output file in the two directories.
    # Comparisons are done with h5diff, which must be in the PATH.
    # Attributes of the steps and other top-level groups are excluded from
    # comparison.
    for filename in files:
        file1 = os.path.join(dir1, filename)
        file2 = os.path.join(dir2, filename)
        if verbose:
            print(f"Numerically comparing {file1} to {file2}.")

        # Compare top-level groups, without attributes.
        # group_names = ["Bx0", "BxD", "By0", "ByD", "Bz0", "BzD",
        #                "X", "Y", "Z", "dV"]
        # for group_name in group_names:
        #     group_path = f"/{group_name}"
        #     if verbose:
        #         print(f"  Comparing {group_path}.")
        #     cmd = (
        #         f"h5diff --exclude-attribute {group_path} {file1} {file2} "
        #         f"{group_path}"
        #     )
        #     cproc = subprocess.run(cmd, shell=True, check=True)
        #     if cproc.returncode != 0:
        #         return TEST_FAIL

        # Compare each step, without attributes.
        _, step_ids = kaiH5.cntSteps(file1)
        for step_id in step_ids:
            step_path = f"/Step#{step_id}"
            if verbose:
                print(f"  Comparing {step_path}.")
            cmd = (
                f"h5diff --exclude-attribute {step_path} {file1} {file2} "
                f"{step_path}"
            )
            cproc = subprocess.run(cmd, shell=True, check=True)
            if cproc.returncode != 0:
                return TEST_FAIL

    # Return the result of the comparison.
    return TEST_PASS


def compare_MHDRCM_results(runxml1: str, runxml2: str, verbose: bool = False):
    """Numerically compare the MHD RCM output files from two runs.

    Numerically compare the MHD RCM output files from two runs.

    Parameters
    ----------
    runxm1 : str
        Path to XML file describing 1st run.
    runxm2 : str
        Path to XML file describing 2nd run.
    verbose : bool
        Set to True to print verbose information during comparison.

    Returns
    -------
    TEST_PASS or TEST_FAIL : str
        Description of result of comparison.

    Raises
    ------
    None
    """
    # Determine the directories containing the sets of results.
    dir1 = os.path.split(runxml1)[0]
    dir2 = os.path.split(runxml2)[0]

    # Generate a sorted list of output files for the 1st run.
    pattern1 = os.path.join(dir1, "*.mhdrcm.h5")
    files1 = glob.glob(pattern1)
    files = [os.path.split(f)[1] for f in files1]
    files.sort()

    # Compare each output file in the two directories.
    # Comparisons are done with h5diff, which must be in the PATH.
    # Attributes of the steps and other top-level groups are excluded from
    # comparison.
    for filename in files:
        file1 = os.path.join(dir1, filename)
        file2 = os.path.join(dir2, filename)
        if verbose:
            print(f"Numerically comparing {file1} to {file2}.")

        # Compare each step, without attributes.
        _, step_ids = kaiH5.cntSteps(file1)
        for step_id in step_ids:
            step_path = f"/Step#{step_id}"
            if verbose:
                print(f"  Comparing {step_path}.")
            cmd = (
                f"h5diff --exclude-attribute {step_path} {file1} {file2} "
                f"{step_path}"
            )
            cproc = subprocess.run(cmd, shell=True, check=True)
            if cproc.returncode != 0:
                return TEST_FAIL

    # Return the result of the comparison.
    return TEST_PASS


def compare_REMIX_results(runxml1: str, runxml2: str, verbose: bool = False):
    """Numerically compare the REMIX output files from two runs.

    Numerically compare the REMIX output files from two runs.

    Parameters
    ----------
    runxm1 : str
        Path to XML file describing 1st run.
    runxm2 : str
        Path to XML file describing 2nd run.
    verbose : bool
        Set to True to print verbose information during comparison.

    Returns
    -------
    TEST_PASS or TEST_FAIL : str
        Description of result of comparison.

    Raises
    ------
    None
    """
    # Determine the directories containing the sets of results.
    dir1 = os.path.split(runxml1)[0]
    dir2 = os.path.split(runxml2)[0]

    # Generate a sorted list of output files for the 1st run.
    pattern1 = os.path.join(dir1, "*.mix.h5")
    files1 = glob.glob(pattern1)
    files = [os.path.split(f)[1] for f in files1]
    files.sort()

    # Compare each result file in the two directories.
    # Comparisons are done with h5diff, which must be in the PATH.
    # Attributes of the steps and other top-level groups are excluded from
    # comparison.
    for filename in files:
        file1 = os.path.join(dir1, filename)
        file2 = os.path.join(dir2, filename)
        if verbose:
            print(f"Numerically comparing {file1} to {file2}.")

        # Compare each step, without attributes.
        _, step_ids = kaiH5.cntSteps(file1)
        for step_id in step_ids:
            step_path = f"/Step#{step_id}"
            if verbose:
                print(f"  Comparing {step_path}.")
            cmd = (
                f"h5diff --exclude-attribute {step_path} {file1} {file2} "
                f"{step_path}"
            )
            cproc = subprocess.run(cmd, shell=True, check=True)
            if cproc.returncode != 0:
                return TEST_FAIL

    # Return the result of the comparison.
    return TEST_PASS


def compare_RCM_results(runxml1: str, runxml2: str, verbose: bool = False):
    """Numerically compare the RCM output files from two runs.

    Numerically compare the RCM output files from two runs.

    Parameters
    ----------
    runxm1 : str
        Path to XML file describing 1st run.
    runxm2 : str
        Path to XML file describing 2nd run.
    verbose : bool
        Set to True to print verbose information during comparison.

    Returns
    -------
    TEST_PASS or TEST_FAIL : str
        Description of result of comparison.

    Raises
    ------
    None
    """
    # Determine the directories containing the sets of results.
    dir1 = os.path.split(runxml1)[0]
    dir2 = os.path.split(runxml2)[0]

    # Generate a sorted list of output files for the 1st run.
    pattern1 = os.path.join(dir1, "*.rcm.h5")
    files1 = glob.glob(pattern1)
    files = [os.path.split(f)[1] for f in files1]
    files.sort()

    # Compare each result file in the two directories.
    # Comparisons are done with h5diff, which must be in the PATH.
    # Attributes of the steps and other top-level groups are excluded from
    # comparison.
    for filename in files:
        file1 = os.path.join(dir1, filename)
        file2 = os.path.join(dir2, filename)
        if verbose:
            print(f"Numerically comparing {file1} to {file2}.")

        # Compare each step, without attributes.
        _, step_ids = kaiH5.cntSteps(file1)
        for step_id in step_ids:
            step_path = f"/Step#{step_id}"
            if verbose:
                print(f"  Comparing {step_path}.")
            cmd = (
                f"h5diff --exclude-attribute {step_path} {file1} {file2} "
                f"{step_path}"
            )
            cproc = subprocess.run(cmd, shell=True, check=True)
            if cproc.returncode != 0:
                return TEST_FAIL

    # Return the result of the comparison.
    return TEST_PASS


def compare_VOLTRON_results(runxml1: str, runxml2: str, verbose: bool = False):
    """Numerically compare the VOLTRON output files from two runs.

    Numerically compare the VOLTRON output files from two runs.

    Parameters
    ----------
    runxm1 : str
        Path to XML file describing 1st run.
    runxm2 : str
        Path to XML file describing 2nd run.
    verbose : bool
        Set to True to print verbose information during comparison.

    Returns
    -------
    TEST_PASS or TEST_FAIL : str
        Description of result of comparison.

    Raises
    ------
    None
    """
    # Determine the directories containing the sets of results.
    dir1 = os.path.split(runxml1)[0]
    dir2 = os.path.split(runxml2)[0]

    # Generate a sorted list of output files for the 1st run.
    pattern1 = os.path.join(dir1, "*.volt.h5")
    files1 = glob.glob(pattern1)
    files = [os.path.split(f)[1] for f in files1]
    files.sort()

    # Compare each result file in the two directories.
    # Comparisons are done with h5diff, which must be in the PATH.
    # Attributes of the steps and other top-level groups are excluded from
    # comparison.
    for filename in files:
        file1 = os.path.join(dir1, filename)
        file2 = os.path.join(dir2, filename)
        if verbose:
            print(f"Numerically comparing {file1} to {file2}.")

        # Compare each step, without attributes.
        _, step_ids = kaiH5.cntSteps(file1)
        for step_id in step_ids:
            step_path = f"/Step#{step_id}"
            if verbose:
                print(f"  Comparing {step_path}.")
            cmd = (
                f"h5diff --exclude-attribute {step_path} {file1} {file2} "
                f"{step_path}"
            )
            cproc = subprocess.run(cmd, shell=True, check=True)
            if cproc.returncode != 0:
                return TEST_FAIL

    # Return the result of the comparison.
    return TEST_PASS


def compare_mage_runs_numerical(args: dict):
    """Numerically compare a pair of MAGE model runs.

    Numerically compare a pair of MAGE model runs.

    Parameters
    ----------
    args : dict
        Dictionary of command-line and other options.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Local convenience variables.
    debug = args.get("debug", False)
    loud = args.get("loud", False)
    run_description = args.get("run_description", None)
    slack_on_fail = args.get("slack_on_fail", False)
    test = args.get("test", False)
    verbose = args.get("verbose", False)
    runxml1 = args["runxml1"]
    runxml2 = args["runxml2"]

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}.")
        print(f"Current directory is {os.getcwd()}.")

    # ------------------------------------------------------------------------

    # Initialize the list of comparison results.
    comparison_results = []

    # Compare the output files from GAMERA.
    if verbose:
        print("Comparing output files from GAMERA.")
    comparison_result = compare_GAMERA_results(runxml1, runxml2, verbose)
    if debug:
        print(f"comparison_result = {comparison_result}")
    comparison_results.append(comparison_result)

    # Compare the MHD RCM output files.
    if verbose:
        print("Comparing MHD RCM output files.")
    comparison_result = compare_MHDRCM_results(runxml1, runxml2, verbose)
    if debug:
        print(f"comparison_result = {comparison_result}")
    comparison_results.append(comparison_result)

    # Compare the REMIX output files.
    if verbose:
        print("Comparing REMIX output files.")
    comparison_result = compare_REMIX_results(runxml1, runxml2, verbose)
    if debug:
        print(f"comparison_result = {comparison_result}")
    comparison_results.append(comparison_result)

    # Compare the RCM output files.
    if verbose:
        print("Comparing RCM output files.")
    comparison_result = compare_RCM_results(runxml1, runxml2, verbose)
    if debug:
        print(f"comparison_result = {comparison_result}")
    comparison_results.append(comparison_result)

    # Compare the VOLTRON output files.
    if verbose:
        print("Comparing VOLTRON output files.")
    comparison_result = compare_VOLTRON_results(runxml1, runxml2, verbose)
    if debug:
        print(f"comparison_result = {comparison_result}")
    comparison_results.append(comparison_result)

    # ------------------------------------------------------------------------

    # Detail the test results.
    test_report_details_string = ""
    test_report_details_string += f"GAMERA: *{comparison_results[0]}*\n"
    test_report_details_string += f"MHD RCM: *{comparison_results[1]}*\n"
    test_report_details_string += f"REMIX: *{comparison_results[2]}*\n"
    test_report_details_string += f"RCM: *{comparison_results[3]}*\n"
    test_report_details_string += f"VOLTRON: *{comparison_results[4]}*\n"

    # Summarize the test results.
    if run_description is not None:
        test_report_summary_string = run_description
    else:
        test_report_summary_string = "Numerical comparison results"
    test_report_summary_string += "\n"
    if "FAIL" in test_report_details_string:
        test_report_summary_string += "*FAIL*"
    else:
        test_report_summary_string += "*PASS*"

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If a test failed, or loud mode is on, post report to Slack.
    if (slack_on_fail and "FAIL" in test_report_summary_string) or loud:
        slack_client = common.slack_create_client()
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=test
        )
        thread_ts = slack_response_summary["ts"]
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_details_string, thread_ts=thread_ts,
            is_test=test
        )

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")


def main():
    """Driver for command-line version of code."""
    # Set up the command-line parser.
    parser = common.create_command_line_parser(DESCRIPTION)

    # Add additional arguments specific to this script.
    parser.add_argument(
        "--run_description", default=None,
        help=(
            "Descriptive string to use in messages from script "
            "(default: %(default)s)"
        )
    )
    parser.add_argument(
        "runxml1", help="Path to XML file describing 1st run."
    )
    parser.add_argument(
        "runxml2", help="Path to XML file describing 2nd run."
    )

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")

    # Convert the arguments from Namespace to dict.
    args = vars(args)
    if args["debug"]:
        print(f"args = {args}")

    # Pass the command-line arguments to the main function as a dict.
    compare_mage_runs_numerical(args)


if __name__ == "__main__":
    main()
