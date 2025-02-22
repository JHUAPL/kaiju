#!/usr/bin/env python

"""Report on the MAGE Fortran unit test results.

This script generates a report on the results of the most recent run of
the MAGE Fortran unit tests.

This script should be run by the PBS script unitTestReport.pbs, to ensure
that the report is only generated when the other unit test jobs are complete.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import glob
import os
import sys

# Import 3rd-party modules.

# Import project modules.
import common


# Program constants

# Program description.
DESCRIPTION = 'Report on the MAGE Fortran unit test results.'

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for unit tests
UNIT_TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'unitTest')

# glob pattern for naming unit test directories
UNIT_TEST_DIRECTORY_GLOB_PATTERN = 'unitTest_*'

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'

# Name of file containing job IDs for each unit test directory.
JOB_ID_LIST_FILE = 'jobs.txt'


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
    slack_on_fail = args.slack_on_fail
    is_test = args.test
    verbose = args.verbose

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Move to the unit test directory.
    os.chdir(UNIT_TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Get list of unit test directories.
    unit_test_directories = glob.glob(UNIT_TEST_DIRECTORY_GLOB_PATTERN)
    if debug:
        print(f"unit_test_directories = {unit_test_directories}")

    # ------------------------------------------------------------------------

    # Initialize result flags.
    myError = False
    jobKilled = False
    okFailure = False
    okCount = 0

    # Check the results in each unit test directory.
    for unit_test_directory in unit_test_directories:
        if verbose:
            print(f"Checking unit test results in {unit_test_directory}.")

        # Move to the directory containing the unit test results.
        path = os.path.join(UNIT_TEST_DIRECTORY, unit_test_directory,
                            BUILD_BIN_DIR)
        if debug:
            print(f"path = {path}")
        os.chdir(path)

        # Check for the job ID list file. If not found, skip the rest of this
        # loop.
        if not os.path.isfile(JOB_ID_LIST_FILE):
            print(
                f"Job list file {JOB_ID_LIST_FILE} in {unit_test_directory}"
                ' not found, skipping report for this directory.'
            )
            continue

        # Read the job IDs from the job ID list file.
        with open(JOB_ID_LIST_FILE, 'r', encoding='utf-8') as f:
            lines = f. readlines()
        job_ids = [_.rstrip() for _ in lines]
        if debug:
            print(f"job_ids = {job_ids}")

        # NOTE: This needs to be reorganized.

        # Compute the names of the job log files.
        job_file_0 = f"genTestData.o{job_ids[0]}"  # 0 OKs
        job_file_1 = f"runCaseTests.o{job_ids[1]}" # 2 OKs
        job_file_2 = f"runNonCaseTests1.o{job_ids[2]}"  # 7 OKs
        job_file_3 = f"runNonCaseTests2.o{job_ids[3]}"  # 1 OK
        if debug:
            print(f"job_file_0 = {job_file_0}")
            print(f"job_file_1 = {job_file_1}")
            print(f"job_file_2 = {job_file_2}")
            print(f"job_file_3 = {job_file_3}")

        # Combine the results of each test log file.
        bigFile = []
        job_files = [job_file_0, job_file_1, job_file_2, job_file_3]
        for job_file in job_files:
            with open(job_file, 'r', encoding='utf-8') as f:
                bigFile += f.readlines()
            bigFile.append('\n\n\n')

        # Scan through for some key things like "error" and "job killed"
        for line in bigFile:
            line = line.rstrip()
            if line == ' OK':
                okCount += 1
            if 'error' in line:
                myError = True
            elif 'job killed' in line:
                jobKilled = True

        # There should be exactly 10 OKs.
        OK_COUNT_EXPECTED = 10
        if verbose:
            print(f"Found {okCount} OKs, expected {OK_COUNT_EXPECTED}.")
        if okCount != OK_COUNT_EXPECTED:
            okFailure = True
        else:
            okFailure = False

        if debug:
            print(f"myError = {myError}")
            print(f"jobKilled = {jobKilled}")
            print(f"okFailure = {okFailure}")
            print(f"okCount = {okCount}")

        # End of unit test directory loop

    # If no tests were complete, say so.
    if len(unit_test_directories) == 0:
        print(
            'Results of unit test report `unitTestReport.py`:\n'
            'No Fortran unit test results were available.'
        )

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Test results are on `derecho` in {os.getcwd()}.\n"
    )
    if myError:
        test_report_details_string += 'Errors occurred during testing.\n'
    if jobKilled:
        test_report_details_string += 'The testing job was killed.\n'
    if okFailure:
        test_report_details_string += 'There was not the correct OK count.\n'
    else:
        test_report_details_string += 'Individual tests: *ALL PASSED*\n'

    # Summarize the test results.
    test_report_summary_string = (
        f"Unit test results for `{os.environ['BRANCH_OR_COMMIT']}`: "
    )
    if myError or jobKilled or okFailure:
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
