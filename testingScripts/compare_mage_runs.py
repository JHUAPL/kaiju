#!/usr/bin/env python

"""Compare a set of MAGE model runs.

This script generates a report which compares the results of MAGE model runs.
A run is specified using the path to the directory containing its results. If
a single run is specified, only that run is plotted. If multiple runs are
specified, the results are compared in each individual plot.

This code is based on the original weeklyDashReport.py script.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import argparse
import datetime
import os
import platform
import subprocess
import sys

# Import 3rd-party modules.
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

# Import project modules.
import common
from kaipy import kaiH5
from kaipy import kaiTools


# Program constants

# Program description.
DESCRIPTION = "Compare MAGE model runs."


def read_real_time_performance(log_path: str):
    """Read real-time performance data from a log file.

    Read real-time performance data from a log file.

    Parameters
    ----------
    log_path : str
        Path to log file to read

    Returns
    -------
    rt_dates, rt_performance : list of datetime, list of float
        Dates and % real-time values for the run.

    Raises
    ------
    None
    """
    # Run the sed command to extract the output times, then convert to
    # datetime objects.
    cmd = (
        "sed --quiet 's/^ \\+UT \\+= \\+\\([0-9-]\\+ [0-9:]\\+\\).*$/\\1/p' "
        f"{log_path}"
    )
    cproc = subprocess.run(cmd, shell=True, check=True,
                           capture_output=True, text=True)
    ut_strings = cproc.stdout
    dates = [datetime.datetime.strptime(ut, "%Y-%m-%d %H:%M:%S")
             for ut in ut_strings.splitlines()]

    # Run the sed command to extract the % real-time performance numbers.
    cmd = (
        "sed --quiet 's/^ \\+Running @ *\\([0-9]\\+\\.\\?[0-9]*\\)% of "
        f"real-time.*$/\\1/p' {log_path}"
    )
    cproc = subprocess.run(cmd, shell=True, check=True,
                           capture_output=True, text=True)
    performance_strings = cproc.stdout
    performance = [float(p) for p in performance_strings.splitlines()]

    # <HACK>
    # Make sure the lists of dates and performance are of equal length
    # (console output not always reliable). Equalize the lengths by truncating
    # the lists at the end.
    if len(dates) > len(performance):
        dates = dates[:len(performance)]
    elif len(performance) > len(dates):
        performance = performance[:len(dates)]
    # </HACK>

    # Return the data.
    return dates, performance


def create_real_time_performance_plot(xml_files: list):
    """Plot real-time performance for each run.

    Plot real-time performance for each run.

    Parameters
    ----------
    xml_files : list of str
        List of XML files describing each run.

    Returns
    -------
    fig : matplotlilb.figure.Figure
        Figure object for plot

    Raises
    ------
    None
    """
    # Read the dates and performance data for each run from the log file.
    git_branches = []
    git_hashes = []
    dates = []
    performance = []
    for xml_file in xml_files:

        # Fetch the branch and hash for this run.
        first_file = common.determine_first_gamera_result_file(xml_file)
        git_branch = kaiH5.GetBranch(first_file)
        git_branches.append(git_branch)
        git_hash = kaiH5.GetHash(first_file)
        git_hashes.append(git_hash)

        # Compute the path to the log file for this run.
        run_directory = os.path.split(xml_file)[0]
        log_path = os.path.join(run_directory, "weeklyDashGo.out")

        # Extract the dates and real-time performance from the log.
        _d, _p = read_real_time_performance(log_path)

        # Save the current date and data.
        dates.append(_d)
        performance.append(_p)

    # Create the plot.
    figsize = (14, 7)
    fig, ax = plt.subplots(figsize=figsize)
    for (_b, _h, _d, _p) in zip(git_branches, git_hashes, dates, performance):
        label = f"{_b} (commit {_h})"
        ax.plot(_d, _p, label=label)
    ax.legend(loc="lower right", fontsize="small")

    # Decorate the x-axis.
    ax.set_xlabel("Date (UTC)")
    # Major ticks on the hour.
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M\n%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    # Minor ticks every 15 minutes.
    ax.xaxis.set_minor_locator(mdates.MinuteLocator([15, 30, 45]))
    # Major and minor grid lines.
    grid_linewidth = 0.75
    grid_alpha = 0.25
    grid_color = "slategrey"
    ax.xaxis.grid(True, which="major", linewidth=grid_linewidth,
                  alpha=grid_alpha, color=grid_color)
    ax.xaxis.grid(True, which="minor", linewidth=grid_linewidth/4,
                  alpha=grid_alpha, color=grid_color)

    # Decorate the y-axis.
    ax.set_ylabel("Percent of Real-Time [%]")
    # Major grid lines only.
    ax.yaxis.grid(True, which="major", linewidth=grid_linewidth,
                  alpha=grid_alpha, color=grid_color)

    # Decorate the plot as a whole
    ax.set_title("Real-Time Performance")

    # Return the figure containing the plot.
    return fig


def create_Dst_plot(xml_files: list):
    """Plot Dst for each run.

    Plot Dst for each run.

    Parameters
    ----------
    xml_files : list of str
        List of XML files describing each run.

    Returns
    -------
    fig : matplotlilb.figure.Figure
        Figure object for plot

    Raises
    ------
    None
    """
    # Read the data for each run from a voltron output file.
    git_branches = []
    git_hashes = []
    dates = []
    BSDst = []
    for xml_file in xml_files:

        # Fetch the branch and hash for this run.
        first_file = common.determine_first_gamera_result_file(xml_file)
        git_branch = kaiH5.GetBranch(first_file)
        git_branches.append(git_branch)
        git_hash = kaiH5.GetHash(first_file)
        git_hashes.append(git_hash)

        # Determine the path to the voltron results file.
        voltron_file = common.determine_voltron_result_file(xml_file)

        # Read the dates and BSDst values from the voltron results file.
        mjds = kaiH5.getTs(voltron_file, aID="MJD")
        dates.append(kaiTools.MJD2UT(mjds))
        BSDst.append(kaiH5.getTs(voltron_file, aID="BSDst"))

    # Create the plot.
    figsize = (14, 7)
    fig, ax = plt.subplots(figsize=figsize)
    for (_b, _h, _d, _dst) in zip(git_branches, git_hashes, dates, BSDst):
        label = f"{_b} (commit {_h})"
        ax.plot(_d, _dst, label=label)
    ax.legend(loc="lower right", fontsize="small")

    # Decorate the x-axis.
    ax.set_xlabel("Date (UTC)")
    # Major ticks on the hour.
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M\n%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    # Minor ticks every 15 minutes.
    ax.xaxis.set_minor_locator(mdates.MinuteLocator([15, 30, 45]))
    # Major and minor grid lines.
    grid_linewidth = 0.75
    grid_alpha = 0.25
    grid_color = "slategrey"
    ax.xaxis.grid(True, which="major", linewidth=grid_linewidth,
                  alpha=grid_alpha, color=grid_color)
    ax.xaxis.grid(True, which="minor", linewidth=grid_linewidth/4,
                  alpha=grid_alpha, color=grid_color)

    # Decorate the y-axis.
    ax.set_ylabel("$D_{st}$")
    # Major grid lines only.
    ax.yaxis.grid(True, which="major", linewidth=grid_linewidth,
                  alpha=grid_alpha, color=grid_color)

    # Decorate the plot as a whole
    ax.set_title("$D_{st}$")

    # Return the figure containing the plot.
    return fig


def create_CPCP_plot(xml_files: list):
    """Plot CPCP for each run.

    Plot CPCP for each run. CPCP is the Cross-Polar Cap Potential.

    Parameters
    ----------
    xml_files : list of str
        List of XML files describing each run.

    Returns
    -------
    fig : matplotlilb.figure.Figure
        Figure object for plot

    Raises
    ------
    None
    """
    # Read the data for each run from a remix output file.
    git_branches = []
    git_hashes = []
    dates = []
    CPCP_north = []
    CPCP_south = []
    for xml_file in xml_files:

        # Fetch the branch and hash for this run.
        first_file = common.determine_first_gamera_result_file(xml_file)
        git_branch = kaiH5.GetBranch(first_file)
        git_branches.append(git_branch)
        git_hash = kaiH5.GetHash(first_file)
        git_hashes.append(git_hash)

        # Determine the path to the remix results file.
        remix_file = common.determine_remix_result_file(xml_file)

        # Read the dataes and CPCP (north and south) values from the remix
        # results file.
        mjds = kaiH5.getTs(remix_file, aID="MJD")
        dates.append(kaiTools.MJD2UT(mjds))
        CPCP_north.append(kaiH5.getTs(remix_file, aID="nCPCP"))
        CPCP_south.append(kaiH5.getTs(remix_file, aID="sCPCP"))

    # Fetch the defaut color cycle.
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]

    # Other plot parameters.
    linewidth = 1.5

    # Create the plot.
    figsize = (14, 7)
    fig, ax = plt.subplots(figsize=figsize)
    ic = 0
    for (_b, _h, _d, _n, _s) in zip(git_branches, git_hashes, dates,
                                    CPCP_north, CPCP_south):
        label = f"North {_b} (commit {_h})"
        ax.plot(_d, _n, label=label, color=colors[ic], linewidth=linewidth)
        label = f"South {_b} (commit {_h})"
        ax.plot(_d, _s, label=label, color=colors[ic], linewidth=linewidth,
                linestyle="dashed")
        ic += 1
    ax.legend(loc="lower right", fontsize="small")

    # Decorate the x-axis.
    ax.set_xlabel("Date (UTC)")
    # Major ticks on the hour.
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M\n%Y-%m-%d"))
    ax.xaxis.set_major_locator(mdates.HourLocator())
    # Minor ticks every 15 minutes.
    ax.xaxis.set_minor_locator(mdates.MinuteLocator([15, 30, 45]))
    # Major and minor grid lines.
    grid_linewidth = 0.75
    grid_alpha = 0.25
    grid_color = "slategrey"
    ax.xaxis.grid(True, which="major", linewidth=grid_linewidth,
                  alpha=grid_alpha, color=grid_color)
    ax.xaxis.grid(True, which="minor", linewidth=grid_linewidth/4,
                  alpha=grid_alpha, color=grid_color)

    # Decorate the y-axis.
    ax.set_ylabel("CPCP (kV)")
    # Major grid lines only.
    ax.yaxis.grid(True, which="major", linewidth=grid_linewidth,
                  alpha=grid_alpha, color=grid_color)

    # Decorate the plot as a whole
    ax.set_title("CPCP")

    # Return the figure containing the plot.
    return fig


def create_magnetosphere_quicklook_plots(xml_files: list):
    """Create the magnetosphere quicklook plot for each run.

    Create the magnetosphere quicklook plot for each run.

    Parameters
    ----------
    xml_files : list of str
        List of XML files describing each run.

    Returns
    -------
    quicklook_plots : list of str
        Path to each quicklook file.

    Raises
    ------
    None
    """
    # Make the magnetosphere quicklook plot for each run.
    quicklook_plots = []
    for xml_file in xml_files:

        # Extract the run ID.
        runid = common.extract_runid(xml_file)

        # Save the current directory.
        cwd = os.getcwd()

        # Compute the path to the results directory, and go there.
        results_dir = os.path.split(xml_file)[0]
        os.chdir(results_dir)

        # Create the quicklook plot.
        cmd = f"msphpic.py -id {runid}"
        _ = subprocess.run(cmd, shell=True, check=True)

        # Move back to the starting directory.
        os.chdir(cwd)

        # Add the plot to the list.
        path = os.path.join(results_dir, "qkmsphpic.png")
        quicklook_plots.append(path)

    # Return the list of quicklook plots.
    return quicklook_plots


def merge_magnetosphere_quicklook_plots(quicklook_plots: list):
    """Merge the magnetosphere quicklook plots for all runs.

    Merge the magnetosphere quicklook plots for all runs.

    Parameters
    ----------
    quicklook_plots : list of str
        List of quicklook plots to merge.

    Returns
    -------
    merged_plot : str
        Path to merged quicklook file.

    Raises
    ------
    None
    """
    # Merge magnetosphere quicklooks.
    merged_plot = "combined_msphpic.png"
    cmd = f"convert {' '.join(quicklook_plots)} -append {merged_plot}"
    print(f"cmd = {cmd}")
    _ = subprocess.run(cmd, shell=True, check=True)
    return merged_plot


def create_REMIX_quicklook_plots(xml_files: list):
    """Create the REMIX quicklook plot for each run.

    Create the REMIX quicklook plot for each run.

    Parameters
    ----------
    xml_files : list of str
        List of XML files describing each run.

    Returns
    -------
    plots_north : list of str
        Path to each northern-hemisphere quicklook file.
    plots_south : list of str
        Path to each southern-hemisphere quicklook file.

    Raises
    ------
    None
    """
    # Save the current directory.
    cwd = os.getcwd()

    # Make the remix quicklook plots for each run.
    plots_north = []
    plots_south = []
    for xml_file in xml_files:

        # Extract the run ID.
        runid = common.extract_runid(xml_file)

        # Compute the path to the results directory.
        results_dir = os.path.split(xml_file)[0]

        # Move to the results directory.
        os.chdir(results_dir)

        # Create the quicklook plots.
        cmd = f"mixpic.py -id {runid}"
        _ = subprocess.run(cmd, shell=True, check=True)

        # Add the plots to the lists.
        path = os.path.join(results_dir, "remix_n.png")
        plots_north.append(path)
        path = os.path.join(results_dir, "remix_s.png")
        plots_south.append(path)

        # Return to the original directory.
        os.chdir(cwd)

    # Return the list of quicklook plots.
    return plots_north, plots_south


def merge_REMIX_quicklook_plots(plots_north: list, plots_south: list):
    """Create the magnetosphere quicklook plot for each run.

    Create the magnetosphere quicklook plot for each run.

    Parameters
    ----------
    quicklook_plots : list of str
        List of quicklook plots to merge.

    Returns
    -------
    merged_plot_north : str
        Path to merged quicklook file for northern hemisphere.
    merged_plot_south : str
        Path to merged quicklook file for southern hemisphere.

    Raises
    ------
    None
    """
    # Merge northern REMIX quicklooks.
    merged_plot_north = "combined_remix_n.png"
    cmd = f"convert {' '.join(plots_north)} -append {merged_plot_north}"
    print(f"cmd = {cmd}")
    _ = subprocess.run(cmd, shell=True, check=True)

    # Merge southern REMIX quicklooks.
    merged_plot_south = "combined_remix_s.png"
    cmd = f"convert {' '.join(plots_south)} -append {merged_plot_south}"
    print(f"cmd = {cmd}")
    _ = subprocess.run(cmd, shell=True, check=True)

    # Return both merged plots.
    return merged_plot_north, merged_plot_south


def create_RCM_quicklook_plots(xml_files: list):
    """Create the RCM quicklook plot for each run.

    Create the RCM quicklook plot for each run.

    Parameters
    ----------
    xml_files : list of str
        List of XML files describing each run.

    Returns
    -------
    quicklook_plots : list of str
        Path to each quicklook file.

    Raises
    ------
    None
    """
    # Save the current directory.
    cwd = os.getcwd()

    # Make the RCM quicklook plot for each run.
    quicklook_plots = []
    for xml_file in xml_files:

        # Extract the run ID.
        runid = common.extract_runid(xml_file)

        # Compute the path to the results directory.
        results_dir = os.path.split(xml_file)[0]

        # Move to the results directory.
        os.chdir(results_dir)

        # Create the quicklook plot.
        cmd = f"rcmpic.py -id {runid}"
        _ = subprocess.run(cmd, shell=True, check=True)

        # Add the plot to the list.
        path = os.path.join(results_dir, "qkrcmpic.png")
        quicklook_plots.append(path)

        # Return to the original directory.
        os.chdir(cwd)

    # Return the list of quicklook plots.
    return quicklook_plots


def merge_RCM_quicklook_plots(quicklook_plots: list):
    """Merge the RCM quicklook plots for all runs.

    Merge the RCM quicklook plots for all runs.

    Parameters
    ----------
    quicklook_plots : list of str
        List of quicklook plots to merge.

    Returns
    -------
    merged_plot : str
        Path to merged quicklook file.

    Raises
    ------
    None
    """
    # Merge RCM quicklooks.
    merged_plot = "combined_qkrcmpic.png"
    cmd = f"convert {' '.join(quicklook_plots)} -append {merged_plot}"
    print(f"cmd = {cmd}")
    _ = subprocess.run(cmd, shell=True, check=True)
    return merged_plot


def compare_mage_runs(args):
    """Compare a set of MAGE model runs.

    Compare a set of MAGE model runs. Generate a report which is optionally
    posted to Slack.

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
    test = args.get("test", False)
    verbose = args.get("verbose", False)

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}.")
        print(f"Current directory is {os.getcwd()}.")

    # ------------------------------------------------------------------------

    # Create the combined list (from file and command line) of XML files
    # which describe runs to compare.

    # If a file containing a list of run XML files (typically reference runs
    # for the comparison) to compare was provided, read it.
    run_xml_files = []
    if args["run_xml_list_file"]:
        with open(args["run_xml_list_file"], "r", encoding="utf-8") as f:
            lines = f.readlines()
        run_xml_files += [line.strip() for line in lines]
    if debug:
        print(f"run_xml_files = {run_xml_files}")

    # Append any run XML files listed on the command line.
    if len(args["run_xml_to_compare"]) > 0:
        run_xml_files += args["run_xml_to_compare"]
    if debug:
        print(f"run_xml_files = {run_xml_files}")

    # run_xml_files is now a list which contains the full path to the XML
    # files which describes each run to compare. The directory containing each
    # XML file is assumed to contain all of the result files for that run.

    # ------------------------------------------------------------------------

    # Create all plots in a memory buffer.
    mpl.use("Agg")

    # Create and save the real-time performance plot.
    if verbose:
        print("Creating real-time performance comparison plot.")
    fig = create_real_time_performance_plot(run_xml_files)
    fig.savefig("perfPlots.png")
    plt.close(fig)

    # Create and save the Dst plot.
    if verbose:
        print("Creating Dst comparison plot.")
    fig = create_Dst_plot(run_xml_files)
    fig.savefig("Dst.png")
    plt.close(fig)

    # Create and save the CPCP (Cross-Polar Cap Potential) plot.
    if verbose:
        print("Creating CPCP comparison plot.")
    fig = create_CPCP_plot(run_xml_files)
    fig.savefig("CPCP.png")
    plt.close(fig)

    # Create the magnetosphere quicklook plots.
    if verbose:
        print("Creating magnetosphere quicklook plots.")
    mag_ql_plots = create_magnetosphere_quicklook_plots(run_xml_files)
    if debug:
        print(f"mag_ql_plots = {mag_ql_plots}")

    # Create the merged magnetosphere quicklook plot.
    if verbose:
        print("Creating merged magnetosphere quicklook plot.")
    merged_mag_ql_plot = merge_magnetosphere_quicklook_plots(
        mag_ql_plots)
    if debug:
        print(f"merged_mag_ql_plot = {merged_mag_ql_plot}")

    # Create the REMIX quicklook plots.
    if verbose:
        print("Creating REMIX quicklook plots.")
    plots_north, plots_south = create_REMIX_quicklook_plots(run_xml_files)
    if debug:
        print(f"plots_north = {plots_north}")
        print(f"plots_south = {plots_south}")

    # Create the merged REMIX quicklook plot.
    if verbose:
        print("Creating merged REMIX quicklook plot.")
    merged_remix_plot_n, merged_remix_plot_s = merge_REMIX_quicklook_plots(
        plots_north, plots_south)
    if debug:
        print(f"merged_remix_plot_n = {merged_remix_plot_n}")
        print(f"merged_remix_plot_s = {merged_remix_plot_s}")

    # Create the RCM quicklook plots.
    if verbose:
        print("Creating RCM quicklook plots.")
    rcm_plots = create_RCM_quicklook_plots(run_xml_files)
    if debug:
        print(f"rcm_plots = {rcm_plots}")

    # Create the merged RCM quicklook plots.
    if verbose:
        print("Creating merged RCM quicklook plot.")
    rcm_merged_plots = merge_RCM_quicklook_plots(rcm_plots)
    if debug:
        print(f"rcm_merged_plots = {rcm_merged_plots}")

    # ------------------------------------------------------------------------

    # Post images to Slack.

    # List the files to post and their comments.
    images_to_post = [
        "perfPlots.png",
        "Dst.png",
        "CPCP.png",
        "combined_msphpic.png",
        "combined_remix_n.png",
        "combined_remix_s.png",
        "combined_qkrcmpic.png",
    ]
    comments_to_post = [
        "Real-Time Performance\n\n",
        "DST Plots\n\n",
        "CPCP Plots\n\n",
        "Magnetosphere Quicklook Comparison Plots\n\n",
        "REMIX (north) Quicklook Comparison Plots\n\n",
        "REMIX (south) Quicklook Comparison Plots\n\n",
        "RCM Quicklook Comparison Plots\n\n",
    ]

    # If loud mode is on, post results to Slack.
    if loud:
        slack_client = common.slack_create_client()
        if run_description is not None:
            message = run_description
        else:
            message = "Comparison plots"
        slack_response = common.slack_send_message(
            slack_client, message, is_test=test)
        if slack_response["ok"]:
            parent_ts = slack_response["ts"]
            message = f"Test results are in `{os.getcwd()}`.\n"
            slack_response = common.slack_send_message(
                slack_client, message, thread_ts=parent_ts, is_test=test)
            for (f, c) in zip(images_to_post, comments_to_post):
                if verbose:
                    print(f"Sending image {f} to Slack.")
                slack_response = common.slack_send_image(
                    slack_client, f, initial_comment=c,
                    thread_ts=parent_ts, is_test=test
                )
        else:
            print("Failed to post parent message and images to Slack.")

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
        "--run_xml_list_file", "-f", default=None,
        help=(
            "Path to text file containing list of runs to compare "
            "(default: %(default)s)"
        )
    )
    parser.add_argument(
        "run_xml_to_compare", nargs=argparse.REMAINDER,
        help="List of XML inputs for MAGE model runs to compare."
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
    compare_mage_runs(args)


if __name__ == "__main__":
    main()
