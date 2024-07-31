#!/usr/bin/env python


"""engage for the MAGE magnetosphere software.

This script performs all of the steps needed to prepare to run a coupled MAGE
with GAMEREA, RCM, and TIEGCM components. By default, this script is interactive - the user
is prompted for each decision  that must be made to prepare for the run, based
on the current "--mode" setting.

The modes are:

"BASIC" (the default) - the user is prompted to set only a small subset of MAGE
parameters. All "INTERMEDIATE"- and "EXPERT"-mode parameters are automatically
set to default values.

"INTERMEDIATE" - The user is prompted for "BASIC" and "INTERMEDIATE"
parameters, with "EXPERT" parameters set to defaults.

"EXPERT" - The user is prompted for *all* adjustable parameters.
"""

# Import standard modules.
import argparse
import copy
import datetime
import json
import os
import subprocess

# Import 3rd-party modules.
import h5py
from jinja2 import Template

# Program constants

# Program description.
DESCRIPTION = "Interactive script to prepare a MAGE magnetosphere model run."

# Indent level for JSON output.
JSON_INDENT = 4

# Path to current kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(KAIJUHOME, "scripts", "makeitso")

# Path to option descriptions file.
OPTION_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_engage_descriptions.json"
)

def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--clobber", action="store_true",
        help="Overwrite existing options file (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--mode", default="BASIC",
        help="User mode (BASIC|INTERMEDIATE|EXPERT) (default: %(default)s)."
    )
    parser.add_argument(
        "--options_path", "-o", default=None,
        help="Path to JSON file of options (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print verbose output (default: %(default)s)."
    )
    return parser

def get_run_option(name, description, mode="BASIC"):
    """Prompt the user for a single run option.

    Prompt the user for a single run option. If no user input is provided,
    the default value is returned for the option. If valids are provided, the
    new value is compared against the valids, and rejected if not in the
    valids list.

    Parameters
    ----------
    name : str, default None
        Name of option
    description : dict, default None
        Dictionary of metadata for the option.
    mode : str
        User experience mode: "BASIC", "INTERMEDIATE", or "ADVANCED".

    Returns
    -------
    value_str : str
        Value of option as a string.

    Raises
    ------
    None
    """
    # Extract prompt, default, and valids.
    level = description["LEVEL"]
    prompt = description.get("prompt", "")
    default = description.get("default", None)
    valids = description.get("valids", None)

    # Compare the current mode to the parameter level setting. If the variable
    # level is higher than the user mode, just use the default.
    if mode == "BASIC" and level in ["INTERMEDIATE", "EXPERT"]:
        return default
    if mode == "INTERMEDIATE" and level in ["EXPERT"]:
        return default

    # If provided, add the valid values in val1|val2 format to the prompt.
    if valids is not None:
        vs = "|".join(valids)
        prompt += f" ({vs})"

    # If provided, add the default to the prompt.
    if default is not None:
        prompt += f" [{default}]"

    # Prompt the user and fetch the input until a good value is provided.
    ok = False
    while not ok:

        # Fetch input from the user.
        option_value = input(f"{prompt}: ")

        # Use the default if no user input provided.
        if option_value == "":
            option_value = default

        # Validate the result. If bad, start over.
        if valids is not None and option_value not in valids:
            print(f"Invalid value for option {name}: {option_value}!")
            continue

        # Keep this result.
        ok = True

    # Return the option as a string.
    return str(option_value)

def prompt_user_for_run_options(args):
    """Prompt the user for run options.

    Prompt the user for run options.

    NOTE: In this function, the complete set of parameters is split up
    into logical groups. This is done partly to make the organization of the
    parameters more obvious, and partly to allow the values of options to
    depend upon previously-specified options.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options

    Returns
    -------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Raises
    ------
    None
    """
    # Save the user mode.
    mode = args.mode

    # Read the dictionary of option descriptions.
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)

    # Initialize the dictionary of program options.
    options = {}

    #-------------------------------------------------------------------------

    # General options for the simulation
    o = options["simulation"] = {}
    od = option_descriptions["simulation"]

    # Prompt for the name of the job.
    for on in ["job_name"]:
        o[on] = get_run_option(on, od[on], mode)


    # Prompt for the start and stop date of the run. This will also be
    # used as the start and stop date of the data in the boundary condition
    # file, which will be created using CDAWeb data.
    for on in ["start_date", "stop_date"]:
        o[on] = get_run_option(on, od[on], mode)

    # Compute the total simulation time in seconds, use as segment duration
    # default.
    date_format = '%Y-%m-%dT%H:%M:%S'
    start_date = o["start_date"]
    stop_date = o["stop_date"]
    t1 = datetime.datetime.strptime(start_date, date_format)
    t2 = datetime.datetime.strptime(stop_date, date_format)
    simulation_duration = (t2 - t1).total_seconds()
    od["segment_duration"]["default"] = str(simulation_duration)

    # Ask if the user wants to split the run into multiple segments.
    # If so, prompt for the segment duration. If not, use the default
    # for the segment duration (which is the simulation duration).
    for on in ["use_segments"]:
        o[on] = get_run_option(on, od[on], mode)
    if o["use_segments"] == "Y":
        for on in ["segment_duration"]:
            o[on] = get_run_option(on, od[on], mode)
    else:
        o["segment_duration"] = od["segment_duration"]["default"]

    # Compute the number of segments based on the simulation duration and
    # segment duration, with 1 additional segment just for spinup. Add 1 if
    # there is a remainder.
    if o["use_segments"] == "Y":
        num_segments = simulation_duration/float(o["segment_duration"])
        if num_segments > int(num_segments):
            num_segments += 1
        num_segments = int(num_segments) + 1
    else:
        num_segments = 1

    # Prompt for the remaining parameters.
    for on in ["gamera_grid_type", "hpc_system"]:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # Return the options dictionary.
    return options


def main():
    """Main program code for makeitso.

    This is the main program code for makeitso.

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
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    clobber = args.clobber
    debug = args.debug
    options_path = args.options_path
    verbose = args.verbose

    # Fetch the run options.
    if options_path:
        # Read the run options from a JSON file.
        with open(options_path, "r", encoding="utf-8") as f:
            options = json.load(f)
    else:
        # Prompt the user for the run options.
        options = prompt_user_for_run_options(args)
    if debug:
        print(f"options = {options}")

    # Save the options dictionary as a JSON file in the current directory.
    path = f"{options['simulation']['job_name']}.json"
    if os.path.exists(path):
        if not clobber:
            raise FileExistsError(f"Options file {path} exists!")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)

if __name__ == "__main__":
    """Begin main program."""
    main()