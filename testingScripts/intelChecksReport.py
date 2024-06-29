#!/usr/bin/env python

"""Create the report for MAGE tests run using Intel Inspector.

This script creates the report for MAGE tests run using Intel Inspector.

This script assumes it is run in the directory containing the test results.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import re
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Create report for Intel Inspector tests.'


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
    is_test = args.test

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Check for for the job list file.
    if not os.path.exists('jobs.txt'):
        print('Nothing to report on.')
        sys.exit(0)

    # Read in the jobs.txt file to get the job numbers
    with open('jobs.txt', 'r', encoding='utf-8') as f:
        job1 = f.readline().rstrip()
        job2 = f.readline().rstrip()
    if debug:
        print(f"job1 = {job1}")
        print(f"job2 = {job2}")

    # Combine the job output files.
    job_file_1 = f"mage_intelCheckSubmitMem.o{job1}"
    job_file_2 = f"mage_intelCheckSubmitThread.o{job2}"
    if debug:
        print(f"job_file_1 = {job_file_1}")
        print(f"job_file_2 = {job_file_2}")

    # Check the the logs exist.
    if not os.path.exists(job_file_1):
        print(f"{job_file_1} not found, tests are incomplete.")
        sys.exit(0)
    if not os.path.exists(job_file_2):
        print(f"{job_file_2} not found, tests are incomplete.")
        sys.exit(0)

    # Set up regex for counting problems found.
    problem_pattern = r'(\d+) new problem\(s\) found'

    # ------------------------------------------------------------------------

    # Search the output logs for errors.

    # Memory
    memory_errors_file = 'combinedMemErrs.txt'
    for _, dirs, _ in os.walk('.'):
        for d in dirs:
            if 'memResults' in d:
                try:
                    memory_check_output = subprocess.check_output(
                        ['inspxe-cl', '-report summary', '-result-dir ' + d,
                         '-s-f memSuppress.sup'],
                        stderr=subprocess.STDOUT, universal_newlines=True)
                except subprocess.CalledProcessError as e:
                    # we need to handle non-zero error code
                    memory_check_output = e.output
                problem_match = re.search(problem_pattern, memory_check_output)
                if not problem_match or int(problem_match.group(1)) > 0:
                    try:
                        memory_check_output = subprocess.check_output(
                            ['inspxe-cl', '-report problems',
                             '-result-dir ' + d, '-s-f memSuppress.sup',
                             '-report-all'],
                            stderr=subprocess.STDOUT, universal_newlines=True)
                    except subprocess.CalledProcessError as e:
                        # we need to handle non-zero error code
                        memory_check_output = e.output
                    with open(
                        memory_errors_file, 'a', encoding='utf-8'
                    ) as f:
                        f.write(memory_check_output)
                        f.write('\n')

    # Thread
    thread_errors_file = 'combinedThreadErrs.txt'
    for _, dirs, _ in os.walk('.'):
        for d in dirs:
            if 'threadResults' in d:
                try:
                    thread_check_output = subprocess.check_output(
                        ['inspxe-cl', '-report summary', '-result-dir ' + d,
                         '-s-f threadSuppress.sup'],
                        stderr=subprocess.STDOUT, universal_newlines=True
                    )
                except subprocess.CalledProcessError as e:
                    # we need to handle non-zero error code
                    thread_check_output = e.output
                problem_match = re.search(problem_pattern, thread_check_output)
                if not problem_match or int(problem_match.group(1)) > 0:
                    try:
                        thread_check_output = subprocess.check_output([
                            'inspxe-cl', '-report problems',
                            '-result-dir ' + d, '-s-f threadSuppress.sup',
                            '-report-all'],
                            stderr=subprocess.STDOUT, universal_newlines=True)
                    except subprocess.CalledProcessError as e:
                        # we need to handle non-zero error code
                        thread_check_output = e.output
                    with open(
                        thread_errors_file, 'a', encoding='utf-8'
                    ) as f:
                        f.write(thread_check_output)
                        f.write('\n')

    # -----------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Intel check results are in {os.getcwd()}.\n"
    )
    test_report_details_string += 'Results of memory tests:\n'
    with open(memory_errors_file, 'r', encoding='utf-8') as f:
        test_report_details_string += f.read()
    test_report_details_string += '\n'
    test_report_details_string += 'Results of thread tests:\n'
    with open(thread_errors_file, 'r', encoding='utf-8') as f:
        test_report_details_string += f.read()

    # Summarize the test results
    test_report_summary_string = (
        'Intel Inspector tests complete (`intelChecksReport.py`).'
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
