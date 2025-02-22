#!/usr/bin/env python

"""Create the MAGE weekly dash test report.

This script creates the MAGE weekly dash test report. This script assumes the
result files are in the current directory.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import platform
import re
import subprocess
import sys

# Import 3rd-party modules.
from astropy.time import Time
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

# Import project modules.
import common
import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv


# Program constants

# Program description.
DESCRIPTION = 'Create the MAGE weekly dash test report.'

# Root of directory tree for all tests.
MAGE_TEST_ROOT = os.environ['MAGE_TEST_ROOT']

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for unit tests
WEEKLY_DASH_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'weeklyDash')

# Glob pattern for individual weekly dash directories
WEEKLY_DASH_DIRECTORY_GLOB_PATTERN = 'weeklyDash_*'

# Regular expression for git hash read from weekly dash output log.
GIT_HASH_PATTERN = 'Git hash   = (.{8})'

# Path to directory containing master-branch reference results.
REFERENCE_RESULTS_DIRECTORY_MASTER = os.path.join(
    MAGE_TEST_ROOT, 'weekly_dash_files', 'reference_results', 'master'
)

# Compute the path to the log file for the master branch reference
# results.
REFERENCE_LOG_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'voltron_mpi.out'
)

# Path to directory containing development-branch reference results.
REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT = os.path.join(
    MAGE_TEST_ROOT, 'weekly_dash_files', 'reference_results', 'development'
)

# Compute the path to the log file for the development branch reference
# results.
REFERENCE_LOG_DEVELOPMENT = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT, 'voltron_mpi.out'
)

# Name of file containg PBS job IDs.
JOB_LIST_FILE = 'jobs.txt'

# String naming branch or commit used in this test.
BRANCH_OR_COMMIT = os.environ['BRANCH_OR_COMMIT']

# Name of voltron output file.
VOLTRON_OUTPUT_FILE = 'msphere.volt.h5'

# Compute the path to the voltron output file for the master branch reference
# results.
VOLTRON_OUTPUT_FILE_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, VOLTRON_OUTPUT_FILE
)

# Compute the path to the voltron output file for the development branch
# reference results.
VOLTRON_OUTPUT_FILE_DEVELOPMENT = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT, VOLTRON_OUTPUT_FILE
)

# Name of remix output file.
REMIX_OUTPUT_FILE = 'msphere.mix.h5'

# Compute the path to the remix output file for the master branch reference
# results.
REMIX_OUTPUT_FILE_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, REMIX_OUTPUT_FILE
)

# Compute the path to the remix output file for the development branch
# reference results.
REMIX_OUTPUT_FILE_DEVELOPMENT = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT, REMIX_OUTPUT_FILE
)

# Compute the paths to the quicklook plots for the master branch.
MAGNETOSPHERE_QUICKLOOK_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'qkmsphpic.png'
)
REMIX_NORTH_QUICKLOOK_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'remix_n.png'
)
REMIX_SOUTH_QUICKLOOK_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'remix_s.png'
)
RCM_QUICKLOOK_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'qkrcmpic.png'
)

# Compute the paths to the quicklook plots for the development branch.
MAGNETOSPHERE_QUICKLOOK_DEVELOPMENT = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'qkmsphpic.png'
)
REMIX_NORTH_QUICKLOOK_DEVELOPMENT = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT, 'remix_n.png'
)
REMIX_SOUTH_QUICKLOOK_DEVELOPMENT = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT, 'remix_s.png'
)
RCM_QUICKLOOK_DEVELOPMENT = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT, 'qkrcmpic.png'
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
    # slack_on_fail = args.slack_on_fail
    is_test = args.test
    verbose = args.verbose

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Read reference results for the master branch.
    if verbose:
        print('Reading reference results for real-time performance for master '
              f" branch from {REFERENCE_LOG_MASTER}.")

    # Read the git hash from the log file.
    if verbose:
        print(f"Reading git hash from {REFERENCE_LOG_MASTER}.")
    git_hash_master = 'XXXXXXXX'
    with open(REFERENCE_LOG_MASTER, 'r', encoding='utf-8') as f:
        for line in f:
            git_hash_match = re.search(GIT_HASH_PATTERN, line)
            if git_hash_match:
                git_hash_master = git_hash_match.group(1)
                break
    if debug:
        print(f"git_hash_master = {git_hash_master}")

    # Read the output times for each Voltron output message (as UT strings)
    # from the log file.
    if verbose:
        print(f"Reading UT from {REFERENCE_LOG_MASTER}.")
    cmd = (
        'sed --quiet "s/^ \\+UT \\+= \\+\\([0-9-]\\+ [0-9:]\\+\\).*$/\\1/p" '
        f"{REFERENCE_LOG_MASTER}"
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: Unable to read UT from {REFERENCE_LOG_MASTER}.\n"
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n",
            file=sys.stderr
        )
        sys.exit(1)
    UT_log_str_master = cproc.stdout
    if debug:
        print(f"UT_log_str_master = {UT_log_str_master}")

    # Convert the UT string to a list of datetime objects.
    UT_log_dt_master = [datetime.datetime.strptime(ut, '%Y-%m-%d %H:%M:%S')
                        for ut in UT_log_str_master.splitlines()]
    if debug:
        print(f"UT_log_dt_master = {UT_log_dt_master}")

    # Read % real-time performance values from the log file.
    if verbose:
        print(f"Reading performance data from {REFERENCE_LOG_MASTER}.")
    cmd = (
        'sed --quiet "s/^ \\+Running @ *\\([0-9]\\+\\.\\?[0-9]*\\)% of '
        f'real-time.*$/\\1/p" {REFERENCE_LOG_MASTER}'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to read performance data from'
            f"{REFERENCE_LOG_MASTER}.\n"
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n",
            file=sys.stderr
        )
        sys.exit(1)
    RT_log_str_master = cproc.stdout
    if debug:
        print(f"RT_log_str_master = {RT_log_str_master}")

    # Convert % real-time values to floats.
    RT_log_f_master = [float(x) for x in RT_log_str_master.splitlines()]
    if debug:
        print(f"RT_log_f_master = {RT_log_f_master}")

    # <HACK>
    # Make sure the lists of UT and RT are of equal length now (console
    # output not always reliable). Equalize the lengths by truncating the
    # lists at the end.
    if len(UT_log_dt_master) > len(RT_log_f_master):
        if verbose:
            print(f"WARNING: UT data from {REFERENCE_LOG_MASTER} is longer"
                  f" than RT data ({len(UT_log_dt_master)} > "
                  f"{len(RT_log_f_master)}); truncating UT data.")
        UT_log_dt_master = UT_log_dt_master[:len(RT_log_f_master)]
    elif len(RT_log_f_master) > len(UT_log_dt_master):
        if verbose:
            print(f"WARNING: RT data from {REFERENCE_LOG_MASTER} is longer"
                  f" than UT data ({len(RT_log_f_master)} > "
                  f"{len(UT_log_dt_master)}); truncating UT data.")
        RT_log_f_master = RT_log_f_master[:len(UT_log_dt_master)]
    # </HACK>

    # ------------------------------------------------------------------------

    # Read reference results for the development branch.
    if verbose:
        print(
            'Reading reference results for real-time performance for '
            'development branch from '
            f"{REFERENCE_RESULTS_DIRECTORY_DEVELOPMENT}."
        )

    # Read the git hash from the log file.
    if verbose:
        print(f"Reading git hash from {REFERENCE_LOG_DEVELOPMENT}.")
    git_hash_development = 'XXXXXXXX'
    with open(REFERENCE_LOG_DEVELOPMENT, 'r', encoding='utf-8') as f:
        for line in f:
            git_hash_match = re.search(GIT_HASH_PATTERN, line)
            if git_hash_match:
                git_hash_development = git_hash_match.group(1)
                break
    if debug:
        print(f"git_hash_development = {git_hash_development}")

    # Read the output times for each Voltron output message (as UT strings)
    # from the log file.
    if verbose:
        print(f"Reading UT from {REFERENCE_LOG_DEVELOPMENT}.")
    cmd = (
        'sed --quiet "s/^ \\+UT \\+= \\+\\([0-9-]\\+ [0-9:]\\+\\).*$/\\1/p" '
        f"{REFERENCE_LOG_DEVELOPMENT}"
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: Unable to read UT from {REFERENCE_LOG_DEVELOPMENT}.\n"
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n",
            file=sys.stderr
        )
        sys.exit(1)
    UT_log_str_development = cproc.stdout
    if debug:
        print(f"UT_log_str_development = {UT_log_str_development}")

    # Convert the UT string to a list of datetime objects.
    UT_log_dt_development = [
        datetime.datetime.strptime(ut, '%Y-%m-%d %H:%M:%S')
        for ut in UT_log_str_development.splitlines()
    ]
    if debug:
        print(f"UT_log_dt_development = {UT_log_dt_development}")

    # Read % real-time performance values from the reference results log.
    if verbose:
        print(f"Reading performance data from {REFERENCE_LOG_DEVELOPMENT}.")
    cmd = (
        'sed --quiet "s/^ \\+Running @ *\\([0-9]\\+\\.\\?[0-9]*\\)% of '
        f'real-time.*$/\\1/p" {REFERENCE_LOG_DEVELOPMENT}'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to read performance data from'
            f"{REFERENCE_LOG_DEVELOPMENT}.\n"
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n",
            file=sys.stderr
        )
        sys.exit(1)
    RT_log_str_development = cproc.stdout
    if debug:
        print(f"RT_log_str_development = {RT_log_str_development}")

    # Convert % real-time values to floats.
    RT_log_f_development = [
        float(x) for x in RT_log_str_development.splitlines()
    ]
    if debug:
        print(f"RT_log_f_development = {RT_log_f_development}")

    # <HACK>
    # Make sure the lists of UT and RT are of equal length now (console
    # output not always reliable). Equalize the lengths by truncating the
    # lists at the end.
    if len(UT_log_dt_development) > len(RT_log_f_development):
        if verbose:
            print(f"WARNING: UT data from {REFERENCE_LOG_DEVELOPMENT}"
                  f"is longer than RT data ({len(UT_log_dt_development)} > "
                  f"{len(RT_log_f_development)}); truncating UT data.")
        UT_log_dt_development = (
            UT_log_dt_development[:len(RT_log_f_development)]
        )
    elif len(RT_log_f_development) > len(UT_log_dt_development):
        if verbose:
            print(f"WARNING: RT data from {REFERENCE_LOG_DEVELOPMENT}"
                  f" is longer than UT data ({len(RT_log_f_development)} > "
                  f"{len(UT_log_dt_development)}); truncating UT data.")
        RT_log_f_development = (
            RT_log_f_development[:len(UT_log_dt_development)]
        )
    # </HACK>

    # ------------------------------------------------------------------------

    # Read results from the latest run.
    if verbose:
        print(f"Reading results for latest run in {os.getcwd()}.")

    # Read in the jobs.txt file to get the job number.
    try:
        with open(JOB_LIST_FILE, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(
            f"ERROR: Unable to open {JOB_LIST_FILE} for weekly dash.\n"
            'See job log for output.\n'
            'Aborting weekly dash report.\n',
            file=sys.stderr
        )
        sys.exit(1)
    job_id = lines[0].rstrip()
    if debug:
        print(f"job_id = {job_id}")

    # <HACK>
    weekly_dash_log_latest = 'weeklyDashGo.out'
    # </HACK>

    # Read the git hash from the log file.
    if verbose:
        print(f"Reading git hash from {weekly_dash_log_latest}.")
    git_hash_latest = 'XXXXXXXX'
    with open(weekly_dash_log_latest, 'r', encoding='utf-8') as f:
        for line in f:
            git_hash_match = re.search(GIT_HASH_PATTERN, line)
            if git_hash_match:
                git_hash_latest = git_hash_match.group(1)
                break
    if debug:
        print(f"git_hash_latest = {git_hash_latest}")

    # Read the output times for each Voltron output message (as UT strings)
    # from the log file.
    if verbose:
        print(f"Reading UT from {weekly_dash_log_latest}.")
    cmd = (
        'sed --quiet "s/^ \\+UT \\+= \\+\\([0-9-]\\+ [0-9:]\\+\\).*$/\\1/p" '
        f"{weekly_dash_log_latest}"
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: Unable to read UT from {weekly_dash_log_latest}.\n"
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n",
            file=sys.stderr
        )
        sys.exit(1)
    UT_log_str_latest = cproc.stdout
    if debug:
        print(f"UT_log_str_latest = {UT_log_str_latest}")

    # Convert the UT string to a list of datetime objects.
    UT_log_dt_latest = [
        datetime.datetime.strptime(ut, '%Y-%m-%d %H:%M:%S')
        for ut in UT_log_str_latest.splitlines()
    ]
    if debug:
        print(f"UT_log_dt_latest = {UT_log_dt_latest}")

    # Read % real-time performance values from the results log.
    if verbose:
        print(f"Reading performance data from {weekly_dash_log_latest}.")
    cmd = (
        'sed --quiet "s/^ \\+Running @ *\\([0-9]\\+\\.\\?[0-9]*\\)% of '
        f'real-time.*$/\\1/p" {weekly_dash_log_latest}'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to read performance data from'
            f"{weekly_dash_log_latest}.\n"
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n",
            file=sys.stderr
        )
        sys.exit(1)
    RT_log_str_latest = cproc.stdout
    if debug:
        print(f"RT_log_str_latest = {RT_log_str_latest}")

    # Convert % real-time values to floats.
    RT_log_f_latest = [
        float(x) for x in RT_log_str_latest.splitlines()
    ]
    if debug:
        print(f"RT_log_f_latest = {RT_log_f_latest}")

    # <HACK>
    # Make sure the lists of UT and RT are of equal length now (console
    # output not always reliable). Equalize the lengths by truncating the
    # lists at the end.
    if len(UT_log_dt_latest) > len(RT_log_f_latest):
        if verbose:
            print(f"WARNING: UT data from {weekly_dash_log_latest}"
                  f"is longer than RT data ({len(UT_log_dt_latest)} > "
                  f"{len(RT_log_f_latest)}); truncating UT data.")
        UT_log_dt_latest = (
            UT_log_dt_latest[:len(RT_log_f_latest)]
        )
    elif len(RT_log_f_latest) > len(UT_log_dt_latest):
        if verbose:
            print(f"WARNING: RT data from {weekly_dash_log_latest}"
                  f" is longer than UT data ({len(RT_log_f_latest)} > "
                  f"{len(UT_log_dt_latest)}); truncating UT data.")
        RT_log_f_latest = (
            RT_log_f_latest[:len(UT_log_dt_latest)]
        )
    # </HACK>

    # ------------------------------------------------------------------------

    # Create all plots in a memory buffer.
    mpl.use('Agg')

    # ------------------------------------------------------------------------

    # Make the real-time performance plot.
    if verbose:
        print('Creating real-time performance plot.')

    # Plot parameters
    line_width = 0.75
    alpha = 0.25  # Transparency
    grid_color = 'slategrey'
    figsize = (14, 7)

    # Create the figure to hold the plot.
    fig = plt.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(1, 1, hspace=0.05, wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])

    # Create the plot.
    ax.plot(UT_log_dt_master, RT_log_f_master,
            label=f"Master ({git_hash_master})",
            linewidth=line_width)
    ax.plot(UT_log_dt_development, RT_log_f_development,
            label=f"Development ({git_hash_development})",
            linewidth=line_width)
    ax.plot(UT_log_dt_latest, RT_log_f_latest,
            label=f"{BRANCH_OR_COMMIT} ({git_hash_latest})",
            linewidth=line_width)
    ax.legend(loc='lower right', fontsize='small')

    # Decorate the x-axis.
    ax.set_xlabel('Date (UTC)')
    # Major ticks on the hour.
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M\n%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    # Minor ticks every 15 minutes.
    # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H:%M'))
    # ax.xaxis.set_minor_locator(mdates.MinuteLocator([15, 30, 45]))
    # Major and minor grid lines.
    ax.xaxis.grid(True, which='major', linewidth=line_width, alpha=alpha,
                  color=grid_color)
    # ax.xaxis.grid(True, which='minor', linewidth=line_width/4, alpha=alpha,
    #               color=grid_color)

    # Decorate the y-axis.
    ax.set_ylabel('Percent of Real-Time [%]')
    # Major grid lines only.
    ax.yaxis.grid(True, which='major', linewidth=line_width, alpha=alpha,
                  color=grid_color)

    # Decorate the plot as a whole
    ax.set_title('Real-Time Performance')

    # Save the plot to a file.
    fOut = 'perfPlots.png'
    kv.savePic(fOut)
    plt.close('all')

    # ------------------------------------------------------------------------

    # Read Dst data from the master-branch reference results.
    if verbose:
        print('Reading reference Dst for master branch from'
              f"{VOLTRON_OUTPUT_FILE_MASTER}.")

    # Read the git hash from the voltron output file.
    if verbose:
        print(f"Reading git hash from {VOLTRON_OUTPUT_FILE_MASTER}.")
    git_hash_master = kh5.GetHash(VOLTRON_OUTPUT_FILE_MASTER)
    if debug:
        print(f"git_hash_master = {git_hash_master}")

    # Read the step count and step IDs from the voltron output file.
    n_steps_master, step_IDs_master = kh5.cntSteps(VOLTRON_OUTPUT_FILE_MASTER)
    if debug:
        print(f"n_steps_master = {n_steps_master}")
        print(f"step_IDs_master = {step_IDs_master}")

    # Read the MJD for each output simulation step.
    MJD_master = kh5.getTs(VOLTRON_OUTPUT_FILE_MASTER, step_IDs_master, 'MJD')
    if debug:
        print(f"MJD_master = {MJD_master}")

    # Convert the floating-point MJD values to UT strings in ISO format.
    UT_str_master = Time(MJD_master, format='mjd').isot
    if debug:
        print(f"UT_str_master = {UT_str_master}")

    # Convert the UT strings to datetime objects.
    UT_dt_master = [datetime.datetime.strptime(ut, '%Y-%m-%dT%H:%M:%S.%f')
                    for ut in UT_str_master]
    if debug:
        print(f"UT_dt_master = {UT_dt_master}")

    # Read the Dst values from the voltron output file.
    Dst_master = kh5.getTs(VOLTRON_OUTPUT_FILE_MASTER, step_IDs_master,
                           'BSDst')
    if debug:
        print(f"Dst_master = {Dst_master}")

    # ------------------------------------------------------------------------

    # Read Dst data from the development-branch reference results.
    if verbose:
        print('Reading reference Dst for development branch from'
              f"{VOLTRON_OUTPUT_FILE_DEVELOPMENT}.")

    # Read the git hash from the voltron output file.
    if verbose:
        print(f"Reading git hash from {VOLTRON_OUTPUT_FILE_DEVELOPMENT}.")
    git_hash_development = kh5.GetHash(VOLTRON_OUTPUT_FILE_DEVELOPMENT)
    if debug:
        print(f"git_hash_development = {git_hash_development}")

    # Read the step count and step IDs from the voltron output file.
    n_steps_development, step_IDs_development = kh5.cntSteps(
        VOLTRON_OUTPUT_FILE_DEVELOPMENT
    )
    if debug:
        print(f"n_steps_development = {n_steps_development}")
        print(f"step_IDs_development = {step_IDs_development}")

    # Read the MJD for each output simulation step.
    MJD_development = kh5.getTs(VOLTRON_OUTPUT_FILE_DEVELOPMENT,
                                step_IDs_development, 'MJD')
    if debug:
        print(f"MJD_development = {MJD_development}")

    # Convert the floating-point MJD values to UT strings in ISO format.
    UT_str_development = Time(MJD_development, format='mjd').isot
    if debug:
        print(f"UT_str_development = {UT_str_development}")

    # Convert the UT strings to datetime objects.
    UT_dt_development = [datetime.datetime.strptime(ut, '%Y-%m-%dT%H:%M:%S.%f')
                         for ut in UT_str_development]
    if debug:
        print(f"UT_dt_development = {UT_dt_development}")

    # Read the Dst values from the voltron output file.
    Dst_development = kh5.getTs(VOLTRON_OUTPUT_FILE_DEVELOPMENT,
                                step_IDs_development, 'BSDst')
    if debug:
        print(f"Dst_development = {Dst_development}")

    # ------------------------------------------------------------------------

    # Read Dst from the latest run.
    if verbose:
        print(f"Reading Dst for latest run in {VOLTRON_OUTPUT_FILE}.")

    # Read the git hash from the voltron output file.
    if verbose:
        print(f"Reading git hash from {VOLTRON_OUTPUT_FILE}.")
    git_hash_latest = kh5.GetHash(VOLTRON_OUTPUT_FILE)
    if debug:
        print(f"git_hash_latest = {git_hash_latest}")

    # Read the step count and step IDs from the voltron output file.
    n_steps_latest, step_IDs_latest = kh5.cntSteps(VOLTRON_OUTPUT_FILE)
    if debug:
        print(f"n_steps_latest = {n_steps_latest}")
        print(f"step_IDs_latest = {step_IDs_latest}")

    # Read the MJD for each output simulation step.
    MJD_latest = kh5.getTs(VOLTRON_OUTPUT_FILE, step_IDs_latest, 'MJD')
    if debug:
        print(f"MJD_latest = {MJD_latest}")

    # Convert the floating-point MJD values to UT strings in ISO format.
    UT_str_latest = Time(MJD_latest, format='mjd').isot
    if debug:
        print(f"UT_str_latest = {UT_str_latest}")

    # Convert the UT strings to datetime objects.
    UT_dt_latest = [datetime.datetime.strptime(ut, '%Y-%m-%dT%H:%M:%S.%f')
                    for ut in UT_str_latest]
    if debug:
        print(f"UT_dt_latest = {UT_dt_latest}")

    # Read the Dst values from the voltron output file.
    Dst_latest = kh5.getTs(VOLTRON_OUTPUT_FILE, step_IDs_latest, 'BSDst')
    if debug:
        print(f"Dst_latest = {Dst_latest}")

    # ------------------------------------------------------------------------

    # Make the DST plot.

    # Create the figure to hold the plot.
    fig = mpl.pyplot.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(1, 1, hspace=0.05, wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])

    # Create the plot.
    ax.plot(UT_dt_master, Dst_master,
            label=f"Master ({git_hash_master})",
            linewidth=2*line_width)
    ax.plot(UT_dt_development, Dst_development,
            label=f"Development ({git_hash_development})",
            linewidth=2*line_width)
    ax.plot(UT_dt_latest, Dst_latest,
            label=f"{BRANCH_OR_COMMIT} ({git_hash_latest})",
            linewidth=2*line_width)
    ax.legend(loc='upper right', fontsize='small')

    # Decorate the x-axis.
    ax.set_xlabel('Date (UTC)')
    # Major ticks on the hour.
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M\n%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    # Minor ticks every 15 minutes.
    # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H:%M'))
    # ax.xaxis.set_minor_locator(mdates.MinuteLocator([15, 30, 45]))
    # Major and minor grid lines.
    ax.xaxis.grid(True, which='major', linewidth=line_width, alpha=alpha,
                  color=grid_color)
    # ax.xaxis.grid(True, which='minor', linewidth=line_width/4, alpha=alpha,
    #               color=grid_color)

    # Decorate the y-axis.
    ax.set_ylabel('Dst [nT]')
    # Major grid lines only.
    ax.yaxis.grid(True, which='major', linewidth=line_width, alpha=alpha,
                  color=grid_color)

    # Decorate the plot as a whole
    ax.set_title('BSDst')

    # Save the plot to a file.
    fOut = 'Dst.png'
    kv.savePic(fOut)
    plt.close('all')

    # ------------------------------------------------------------------------

    # Read CPCP (north and south) data from the master-branch reference
    # results.
    if verbose:
        print('Reading reference CPCP (north and south) for master branch '
              f"from {REMIX_OUTPUT_FILE_MASTER}.")

    # Read the CPCP values from the voltron output file.
    CPCP_north_master = kh5.getTs(REMIX_OUTPUT_FILE_MASTER, step_IDs_master,
                                  'nCPCP')
    CPCP_south_master = kh5.getTs(REMIX_OUTPUT_FILE_MASTER, step_IDs_master,
                                  'sCPCP')
    if debug:
        print(f"CPCP_north_master = {CPCP_north_master}")
        print(f"CPCP_south_master = {CPCP_south_master}")

    # ------------------------------------------------------------------------

    # Read CPCP (north and south) data from the development-branch reference
    # results.
    if verbose:
        print('Reading reference CPCP (north and south) for development '
              f"branch from {REMIX_OUTPUT_FILE_DEVELOPMENT}.")

    # Read the CPCP values from the voltron output file.
    CPCP_north_development = kh5.getTs(REMIX_OUTPUT_FILE_DEVELOPMENT,
                                       step_IDs_development, 'nCPCP')
    CPCP_south_development = kh5.getTs(REMIX_OUTPUT_FILE_DEVELOPMENT,
                                       step_IDs_development, 'sCPCP')
    if debug:
        print(f"CPCP_north_development = {CPCP_north_development}")
        print(f"CPCP_south_development = {CPCP_south_development}")

    # ------------------------------------------------------------------------

    # Read CPCP (north and south) data from the latest run.
    if verbose:
        print('Reading CPCP (north and south) for latest run'
              f" from {REMIX_OUTPUT_FILE}.")

    # Read the CPCP values from the voltron output file.
    CPCP_north_latest = kh5.getTs(REMIX_OUTPUT_FILE,
                                  step_IDs_latest, 'nCPCP')
    CPCP_south_latest = kh5.getTs(REMIX_OUTPUT_FILE,
                                  step_IDs_latest, 'sCPCP')
    if debug:
        print(f"CPCP_north_latest = {CPCP_north_latest}")
        print(f"CPCP_south_latest = {CPCP_south_latest}")

    # ------------------------------------------------------------------------

    # Make the CPCP plot.

    # Create the figure to hold the plot.
    fig = mpl.pyplot.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(1, 1, hspace=0.05, wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])

    # Fetch the defaut color cycle.
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    # Create the plot.
    ax.plot(UT_dt_master, CPCP_north_master,
            label=f"Master (north) ({git_hash_master})",
            color=colors[0], linewidth=2*line_width)
    ax.plot(UT_dt_master, CPCP_south_master,
            label=f"Master (south) ({git_hash_master})",
            color=colors[0], linewidth=2*line_width, linestyle='dashed')
    ax.plot(UT_dt_development, CPCP_north_development,
            label=f"Development (north) ({git_hash_development})",
            color=colors[1], linewidth=2*line_width)
    ax.plot(UT_dt_development, CPCP_south_development,
            label=f"Development (south) ({git_hash_development})",
            color=colors[1], linewidth=2*line_width, linestyle='dashed')
    ax.plot(UT_dt_latest, CPCP_north_latest,
            label=f"{BRANCH_OR_COMMIT} (north) ({git_hash_latest})",
            color=colors[2], linewidth=2*line_width)
    ax.plot(UT_dt_latest, CPCP_south_latest,
            label=f"{BRANCH_OR_COMMIT} (south) ({git_hash_latest})",
            color=colors[2], linewidth=2*line_width, linestyle='dashed')
    ax.legend(loc='upper right', fontsize='small')

    # Decorate the x-axis.
    ax.set_xlabel('Date (UTC)')
    # Major ticks on the hour.
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M\n%Y-%m-%d'))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    # Minor ticks every 15 minutes.
    # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H:%M'))
    # ax.xaxis.set_minor_locator(mdates.MinuteLocator([15, 30, 45]))
    # Major and minor grid lines.
    ax.xaxis.grid(True, which='major', linewidth=line_width, alpha=alpha,
                  color=grid_color)
    # ax.xaxis.grid(True, which='minor', linewidth=line_width/4, alpha=alpha,
    #               color=grid_color)

    # Decorate the y-axis.
    ax.set_ylabel('CPCP [kV]')
    # Major grid lines only.
    ax.yaxis.grid(True, which='major', linewidth=line_width, alpha=alpha,
                  color=grid_color)

    # Decorate the plot as a whole
    ax.set_title('CPCP')

    # Save the plot to a file.
    fOut = 'CPCP.png'
    kv.savePic(fOut)
    plt.close('all')

    # ------------------------------------------------------------------------

    # Make the magnetosphere quick-look plot.
    if verbose:
        print('Creating magnetosphere quicklook plot for '
              f"{os.getcwd()}.")

    # Create the plot.
    cmd = 'msphpic.py'
    if debug:
        print(f"cmd = {cmd}")
    try:
        _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to create magnetosphere quicklook plot.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )

    # ------------------------------------------------------------------------

    # Make the REMIX quick-look plots.
    if verbose:
        print(f"Creating REMIX quicklook plots for {os.getcwd()}.")

    # Create the plot.
    cmd = 'mixpic.py'
    if debug:
        print(f"cmd = {cmd}")
    try:
        _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to create REMIX quicklook plots.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )

    # ------------------------------------------------------------------------

    # Make the RCM quick-look plot.
    if verbose:
        print(f"Creating RCM quicklook plot for {os.getcwd()}.")

    # Create the plot.
    cmd = 'rcmpic.py'
    if debug:
        print(f"cmd = {cmd}")
    try:
        _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to create RCM quicklook plot.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )

    # ------------------------------------------------------------------------

    # Create merged images for the quicklook plots.

    # Merge magnetosphere quicklooks.
    cmd = (
        f"convert {MAGNETOSPHERE_QUICKLOOK_MASTER}"
        f" {MAGNETOSPHERE_QUICKLOOK_DEVELOPMENT}"
        ' qkmsphpic.png -append combined_msphpic.png'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to combine magnetosphere quicklook plots.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )

    # Merge REMIX (north) quicklooks.
    cmd = (
        f"convert {REMIX_NORTH_QUICKLOOK_MASTER}"
        f" {REMIX_NORTH_QUICKLOOK_DEVELOPMENT}"
        ' remix_n.png -append combined_remix_n.png'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to combine REMIX (north) quicklook plots.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )

    # Merge REMIX (south) quicklooks.
    cmd = (
        f"convert {REMIX_SOUTH_QUICKLOOK_MASTER}"
        f" {REMIX_SOUTH_QUICKLOOK_DEVELOPMENT}"
        ' remix_s.png -append combined_remix_s.png'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to combine REMIX (south) quicklook plots.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )

    # Merge RCM quicklooks.
    cmd = (
        f"convert {RCM_QUICKLOOK_MASTER}"
        f" {RCM_QUICKLOOK_DEVELOPMENT}"
        ' qkrcmpic.png -append combined_qkrcmpic.png'
    )
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to combine RCM quicklook plots.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )

    # ------------------------------------------------------------------------

    # List the files to post and their comments.
    images_to_post = [
        'perfPlots.png',
        'Dst.png',
        'CPCP.png',
        'qkmsphpic.png',
        'remix_n.png',
        'remix_s.png',
        'qkrcmpic.png',
        'combined_msphpic.png',
        'combined_remix_n.png',
        'combined_remix_s.png',
        'combined_qkrcmpic.png'
    ]
    comments_to_post = [
        'Real-Time Performance\n\n',
        'DST Plots\n\n',
        'CPCP Plots\n\n',
        'Magnetosphere Quicklook Plots\n\n',
        'REMIX (north) Quicklook Plots\n\n',
        'REMIX (south) Quicklook Plots\n\n',
        'RCM Quicklook Plots\n\n',
        'Magnetosphere Quicklook Comparison Plots\n\n',
        'REMIX (north) Quicklook Comparison Plots\n\n',
        'REMIX (south) Quicklook Comparison Plots\n\n',
        'RCM Quicklook Comparison Plots\n\n'
    ]

    # If loud mode is on, post results to Slack.
    if be_loud:
        slack_client = common.slack_create_client()
        if debug:
            print(f"slack_client = {slack_client}")
        message = (
            f"Weekly dash result plots complete for `{BRANCH_OR_COMMIT}`.\n"
        )
        slack_response = common.slack_send_message(
            slack_client, message, is_test=is_test)
        if slack_response['ok']:
            parent_ts = slack_response['ts']
            message = f"Test results are in {os.getcwd()}.\n"
            message += (
                'This was a 4x4x1 (IxJxK) decomposed Quad Resolution Run using'
                ' 8 nodes for Gamera, 1 for Voltron, and 2 Squish Helper nodes'
                ' (11 nodes total).'
            )
            slack_response = common.slack_send_message(
                slack_client, message, thread_ts=parent_ts, is_test=is_test)
            for (f, c) in zip(images_to_post, comments_to_post):
                slack_response = common.slack_send_image(
                    slack_client, f, initial_comment=c,
                    thread_ts=parent_ts, is_test=is_test
                )
        else:
            print('Failed to post parent message and images to Slack.')

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")


if __name__ == '__main__':
    main()
