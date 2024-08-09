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
import sys
import subprocess

# Import 3rd-party modules.
import netCDF4
import h5py
from jinja2 import Template

#import tiegcmrun stuff
#sys.path.append('/glade/u/home/wiltbemj/src/tiegcm/tiegcmrun')
sys.path.append('/glade/u/home/nikhilr/kaiju_engage/tiegcm/')
from tiegcmrun import tiegcmrun
print(f'tiegcmrum from {tiegcmrun.__file__}')

#import makeitso
sys.path.append('/glade/u/home/wiltbemj/src/kaiju-private/scripts/makeitso')
import makeitso
print(f'makeitso from {makeitso.__file__}')
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

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template-gtr.pbs")

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

# def makeitso.get_run_option(name, description, mode="BASIC"):
#     """Prompt the user for a single run option.

#     Prompt the user for a single run option. If no user input is provided,
#     the default value is returned for the option. If valids are provided, the
#     new value is compared against the valids, and rejected if not in the
#     valids list.

#     Parameters
#     ----------
#     name : str, default None
#         Name of option
#     description : dict, default None
#         Dictionary of metadata for the option.
#     mode : str
#         User experience mode: "BASIC", "INTERMEDIATE", or "ADVANCED".

#     Returns
#     -------
#     value_str : str
#         Value of option as a string.

#     Raises
#     ------
#     None
#     """
#     # Extract prompt, default, and valids.
#     level = description["LEVEL"]
#     prompt = description.get("prompt", "")
#     default = description.get("default", None)
#     valids = description.get("valids", None)

#     # Compare the current mode to the parameter level setting. If the variable
#     # level is higher than the user mode, just use the default.
#     if mode == "BASIC" and level in ["INTERMEDIATE", "EXPERT"]:
#         return default
#     if mode == "INTERMEDIATE" and level in ["EXPERT"]:
#         return default

#     # If provided, add the valid values in val1|val2 format to the prompt.
#     if valids is not None:
#         vs = "|".join(valids)
#         prompt += f" ({vs})"

#     # If provided, add the default to the prompt.
#     if default is not None:
#         prompt += f" [{default}]"

#     # Prompt the user and fetch the input until a good value is provided.
#     ok = False
#     while not ok:

#         # Fetch input from the user.
#         option_value = input(f"{prompt}: ")

#         # Use the default if no user input provided.
#         if option_value == "":
#             option_value = default

#         # Validate the result. If bad, start over.
#         if valids is not None and option_value not in valids:
#             print(f"Invalid value for option {name}: {option_value}!")
#             continue

#         # Keep this result.
#         ok = True

#     # Return the option as a string.
#     return str(option_value)

def create_pbs_scripts(gr_options: dict, tiegcm_options: dict):
    """Create the PBS job scripts for the run.

    Create the PBS job scripts from a template.

    Parameters
    ----------
    gr_options : dict
        Dictionary of program options from makeitso, each entry maps str to str.
    gr_options : dict
        Dictionary of program options from tiegcmrun, each entry maps str to str.

    Returns
    -------
    pbs_scripts : list of str
        Paths to PBS job script.
    submit_all_jobs_script : str
        Path to script which submits all PBS jobs.

    Raises
    ------
    TypeError:
        For a non-integral of nodes requested
    """
    # Read the template.
    with open(PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)

    # Create a PBS script for each segment.
    pbs_scripts = []
    for job in range(int(options["pbs"]["num_segments"])):
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        segment_id = f"{runid}-{job:02d}"
        opt["simulation"]["segment_id"] = segment_id
        pbs_content = template.render(opt)
        pbs_script = os.path.join(
            opt["pbs"]["run_directory"],
            f"{opt['simulation']['segment_id']}.pbs"
        )
        pbs_scripts.append(pbs_script)
        with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)

    # Create a single script which will submit all of the PBS jobs in order.
    submit_all_jobs_script = f"{options['simulation']['job_name']}_pbs.sh"
    with open(submit_all_jobs_script, "w", encoding="utf-8") as f:
        s = pbs_scripts[0]
        cmd = f"job_id=`qsub {s}`\n"
        f.write(cmd)
        cmd = "echo $job_id\n"
        f.write(cmd)
        for s in pbs_scripts[1:]:
            cmd = "old_job_id=$job_id\n"
            f.write(cmd)
            cmd = f"job_id=`qsub -W depend=afterok:$old_job_id {s}`\n"
            f.write(cmd)
            cmd = "echo $job_id\n"
            f.write(cmd)

    # Return the paths to the PBS scripts.
    return pbs_scripts, submit_all_jobs_script

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
        o[on] = makeitso.get_run_option(on, od[on], mode)


    # Prompt for the start and stop date of the run. This will also be
    # used as the start and stop date of the data in the boundary condition
    # file, which will be created using CDAWeb data.
    for on in ["start_date", "stop_date"]:
        o[on] = makeitso.get_run_option(on, od[on], mode)

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
        o[on] = makeitso.get_run_option(on, od[on], mode)
    if o["use_segments"] == "Y":
        for on in ["segment_duration"]:
            o[on] = makeitso.get_run_option(on, od[on], mode)
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
        o[on] = makeitso.get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------
    # coupling options
    options["coupling"] = {}
    o = options["coupling"]
    od = option_descriptions["coupling"]

    # Prompt for the remaining parameters.
    for on in ["gamera_spin_up_time", "gcm_spin_up_time", 
               "root_directory", "dtOut"]:
        o[on] = makeitso.get_run_option(on, od[on], mode)
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

    makeitso_args = {'clobber': True, 'debug': True, 'verbose': True}
    makeitso_args.update(options)
    print(f"makeitso_args = {makeitso_args}")
    makeitso_select_line, makeitso_exec_command = makeitso.makeitso(makeitso_args)
    print(f"makeitso_select_line = {makeitso_select_line}")
    print(f"makeitso_exec_command = {makeitso_exec_command}") 

    # Run the TIEGCMrun
    print(options)
    arguments = [
    "--coupling",
    "--engage", json.dumps(options)
]
    tiegcmrun.tiegcmrun(arguments)
    

if __name__ == "__main__":
    """Begin main program."""
    main()