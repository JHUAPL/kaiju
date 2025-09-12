#!/usr/bin/env python


"""makeitso for the MAGE magnetosphere software.

This script performs all of the steps needed to prepare to run a MAGE
magnetosphere simulation run. By default, this script is interactive - the user
is prompted for each decision  that must be made to prepare for the run, based
on the current "--mode" setting.

The modes are:

"BASIC" (the default) - the user is prompted to set only a small subset of MAGE
parameters. All "INTERMEDIATE"- and "EXPERT"-mode parameters are automatically
set to default values.

"INTERMEDIATE" - The user is prompted for "BASIC" and "INTERMEDIATE"
parameters, with "EXPERT" parameters set to defaults.

"EXPERT" - The user is prompted for *all* adjustable parameters.

This script can be run on the command line in the usual fashion, or can be
run as a function after importing this module. In the latter case, the
makeitso() function accepts a dict of arguments containing the command-line
options, as well as optional additional settings from the caller.
"""


# Import standard modules.
import argparse
import copy
import datetime
import json
import os
import subprocess
import math

# Import 3rd-party modules.
import h5py
from jinja2 import Template

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Interactive script to prepare a MAGE magnetosphere model run"

# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(os.environ["KAIJUHOME"], "scripts",
                                       "makeitso")

# Path to option descriptions file.
OPTION_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
)

# Path to template .ini file.
INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.ini")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")

# Indent level for JSON output.
JSON_INDENT = 4

# Default values for command-line arguments when none are supplied (such as
# when makeitso() is called by external code).
args_default = {
    "clobber": False,
    "debug": False,
    "mode": "BASIC",
    "options_path": None,
    "verbose": False,
}


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
        "--clobber", action="store_true", default=args_default["clobber"],
        help="Overwrite existing options file (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", default=args_default["debug"],
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--mode",  default=args_default["mode"],
        help="User mode (BASIC|INTERMEDIATE|EXPERT) (default: %(default)s)."
    )
    parser.add_argument(
        "--options_path", "-o", default=args_default["options_path"],
        help="Path to JSON file of options (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        default=args_default["verbose"],
        help="Print verbose output (default: %(default)s)."
    )
    return parser

def update_loaded_options(options):
    """Update a loaded set of run options to account for format changes.

    Update a loaded set of run options to account for format changes.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    options : dict
        Updated dictionary of program options.

    Raises
    ------
    None
    """
    # <HACK>
    # If the tsOut attribute is found, rename it to dtCon.
    DEFAULT_GEO_DTCON = "300"  # Seconds
    if "tsOut" in options["voltron"]["output"]:
        print("WARNING: Replacing obsolete parameter tsOut with default dtCon"
              f" value of {DEFAULT_GEO_DTCON} seconds!")
        options["voltron"]["output"]["dtCon"] = DEFAULT_GEO_DTCON
        del options["voltron"]["output"]["tsOut"]
    # </HACK>

    # Return the updated options.
    return options

    
def load_option_descriptions(path: str = OPTION_DESCRIPTIONS_FILE):
    """Load the option descriptions and update as needed.

    Read the option descriptions for makeitso from a JSON file.

    Parameters
    ----------
    path : str, default OPTION_DESCRIPTIONS_FILE
        Path to file containing the descriptions for makeitso options.

    Returns
    -------
    option_descriptions : dict
        Dictionary containing descriptions for all makeitso options.

    Raises
    ------
    None
    """
    # Read the dictionary of option descriptions.
    with open(path, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)

    # Return the option descriptions.
    return option_descriptions


def update_option_descriptions(option_descriptions: dict, args: dict):
    """Update the option descriptions with data from the caller.

    Update the option descriptions for makeitso using data provided by the
    calling code. This is typically done to allow engage.py to call makeitso()
    and override parts of the option descriptions, especially defaults.

    Parameters
    ----------
    option_descriptions: dict
        Dictionary of option descriptions read from JSON file.
    args : dict
        Dictionary of command-line options and equivalent options passed from
        the calling function, and variables set by engage for makeitso.

    Returns
    -------
    option_descriptions : dict
        Dictionary containing updated descriptions for all makeitso options.

    Raises
    ------
    None
    """
    # Update the option_descriptions dict based on data provided by engage
    # (if any). If engage provides a BASIC option, use that value without
    # prompting the user. For an INTERMEDIATE or EXPERT option, change the
    # default value for the option to the value provided by engage.

    # Check for simulation options provided by engage. Since all of
    # these options are BASIC-level, copy the new values into the default
    # for the corresponding makeitso simulation option, and delete the
    # prompt string from the description so that the default for the
    # option will be automatically used.
    if "simulation" in args:
        simulation = args["simulation"]
        for k in simulation:
            od = option_descriptions["simulation"][k]
            od["prompt"] = None
            od["default"] = simulation[k]

    # If TIEGCM coupling is specified, the prompt for stat date is set to None.
    # Note that if "coupling" exists, "simulation" *MUST* also exist, 
    # and *MUST* contain the dict keys shown below.
    if "coupling" in args:
        coupling = args["coupling"]
        simulation = args["simulation"]
        gr_warm_up_time = float(coupling["gr_warm_up_time"])
        dt = datetime.timedelta(seconds=gr_warm_up_time)
        start_date = simulation["start_date"]
        t0 = datetime.datetime.fromisoformat(start_date)
        t0 -= dt
        start_date = datetime.datetime.isoformat(t0)
        option_descriptions["simulation"]["start_date"]["prompt"] = None
        option_descriptions["simulation"]["start_date"]["default"] = (
            start_date
        )

    # Incorporate any BASIC PBS options from engage.
    if "pbs" in args:
        pbs = args["pbs"]
        for k in ["account_name", "kaiju_install_directory",
                 "kaiju_build_directory", "run_directory"]:
            od = option_descriptions["pbs"]["_common"][k]
            od["prompt"] = None
            od["default"] = pbs[k]
            
        # Incorporate HPC platform-specific PBS options from engage at basic level.
        hpc_platform = args["simulation"]["hpc_system"]
        for k in args["pbs"]:
            if k in option_descriptions["pbs"][hpc_platform] and option_descriptions["pbs"][hpc_platform][k]["LEVEL"] == "BASIC":
                od = option_descriptions["pbs"][hpc_platform][k]
                od["prompt"] = None
                od["default"] = pbs[k]
        # Incorporate HPC platform-specific PBS options from engage in default.
        option_descriptions["pbs"][hpc_platform]["modules"]["default"] = pbs["modules"]
        if hpc_platform == "pleiades":
            option_descriptions["pbs"][hpc_platform]["moduledir"]["default"] = pbs["moduledir"]
            option_descriptions["pbs"][hpc_platform]["local_modules"]["default"] = pbs["local_modules"]
    # Return the option descriptions.
    return option_descriptions


def get_run_option(name: str, description: dict, mode: str):
    """Prompt the user for a single run option.

    Prompt the user for a single run option. If no user input is provided, or
    there is no prompt string in the description, return the default value for
    the option. If valids are provided, the new value is compared against the
    valids, and rejected if not in the valids list.

    Parameters
    ----------
    name : str
        Name of option
    description : dict
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
    # Extract option level, prompt, default, and valids from the option
    # description.
    level = description["LEVEL"]
    prompt = description.get("prompt", None)
    default = description.get("default", None)
    valids = description.get("valids", None)

    # If there is no prompt, use the default.
    if prompt is None:
        return default

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


def fetch_bcwind_time_range(bcwind_path: str):
    """Fetch the start and stop times for a bcwind file.

    Fetch the start and stop times for a bcwind file.

    Parameters
    ----------
    bcwind_path : str
        Path to bcwind file

    Returns
    -------
    start_date, stop_date : str
        First and last entries in UT group, as strings, in
        'YYYY-MM-DDTHH:MM:SS' format.

    Raises
    ------
    None
    """
    with h5py.File(bcwind_path, "r") as f:
        start_date = f["UT"][0].decode("utf-8")
        stop_date = f["UT"][-1].decode("utf-8")
        # <HACK> Convert from "YYYY-MM-DD HH:MM:SS" format to
        # "YYYY-MM-DDTHH:MM:SS" (ISO 8601) format.
        start_date = start_date.replace(" ", "T")
        stop_date = stop_date.replace(" ", "T")
        # </HACK>
    return start_date, stop_date


def prompt_user_for_run_options(option_descriptions: dict, args: dict):
    """Prompt the user for run options.

    Prompt the user for run options.

    NOTE: In this function, the complete set of parameters is split up
    into logical groups. This is done partly to make the organization of the
    parameters more obvious, and partly to allow the values of options to
    depend upon previously-specified options. This is messy, but it is the
    only way to account for arbitrary dependencies among the various program
    options, expecially when they are modified by data passed from engage.py.

    Parameters
    ----------
    option_descriptions : dict
        Dictionary of option descriptions
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
    # Local convenience variables.
    debug = args["debug"]
    mode = args["mode"]
    verbose = args["verbose"]

    # Initialize the dictionary of program options.
    options = {}

    # ------------------------------------------------------------------------

    # General options for the simulation
    o = options["simulation"] = {}
    od = option_descriptions["simulation"]

    # Prompt for the name of the job.
    for on in ["job_name"]:
        o[on] = get_run_option(on, od[on], mode)

    # Ask the user if a solar wind boundary condition file is available. If
    # so, read the start and stop date from the boundary condition file. If
    # not, prompt the user for the start and stop date so that a boundary
    # condition file can be generated.
    for on in ["bcwind_available"]:
        o[on] = get_run_option(on, od[on], mode)
    if o["bcwind_available"].upper() == "Y":
        for on in ["bcwind_file"]:
            o[on] = get_run_option(on, od[on], mode)
        # Fetch the start and stop date from the bcwind file.
        start_date, stop_date = fetch_bcwind_time_range(o[on])
        o["start_date"] = start_date
        o["stop_date"] = stop_date
    else:
        # Prompt for the start and stop date of the run. This will also be
        # used as the start and stop date of the data in the boundary condition
        # file, which will be created using CDAWeb data.
        # Note that the start date should already have been updated to include
        # the warmup period needed when makeitso is called from engage.
        for on in ["start_date", "stop_date"]:
            o[on] = get_run_option(on, od[on], mode)

    # Compute the total simulation time in seconds, use as segment duration
    # default, if not already specified. This should include any warmup
    # period requested by engage, but should NOT include any spinup (t < 0)
    # period.
    start_date = o["start_date"]
    stop_date = o["stop_date"]
    t1 = datetime.datetime.fromisoformat(start_date)
    t2 = datetime.datetime.fromisoformat(stop_date)
    simulation_duration = (t2 - t1).total_seconds()
    if "default" not in od["segment_duration"]:
        od["segment_duration"]["default"] = str(simulation_duration)

    # Ask if the user wants to split the run into multiple segments.
    # If so, prompt for the segment duration. If not, use the default
    # for the segment duration.
    for on in ["use_segments"]:
        o[on] = get_run_option(on, od[on], mode)
    if o["use_segments"].upper() == "Y":
        for on in ["segment_duration"]:
            o[on] = get_run_option(on, od[on], mode)
    else:
        o["segment_duration"] = od["segment_duration"]["default"]

    # Compute the number of segments based on the simulation duration and
    # segment duration, with 1 separate segment just for spinup. Add 1 if
    # there is a remainder.
    num_segments = 1
    if o["use_segments"].upper() == "Y":
        # Compute the number of simulation segments (t > 0).
        num_segments = simulation_duration/float(o["segment_duration"])
        if num_segments > int(num_segments):
            num_segments = int(num_segments) + 1
        else:
            num_segments = int(num_segments)

    # Prompt for the remaining parameters.
    for on in ["gamera_grid_type", "gamera_grid_inner_radius",
               "gamera_grid_outer_radius", "hpc_system"]:
        o[on] = get_run_option(on, od[on], mode)

    # ------------------------------------------------------------------------

    # PBS options
    o = options["pbs"] = {}

    # Common (HPC platform-independent) options
    od = option_descriptions["pbs"]["_common"]
    
    if "default" not in od.get("account_name", {}):
        od["account_name"]["default"] = os.getlogin()
    if "default" not in od.get("kaiju_install_directory", {}):    
        od["kaiju_install_directory"]["default"] = os.environ["KAIJUHOME"]
    if "default" not in od.get("kaiju_build_directory", {}):
        od["kaiju_build_directory"]["default"] = os.path.join(
        os.environ["KAIJUHOME"], "build_mpi")
    od["num_segments"]["default"] = str(num_segments)
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # HPC platform-specific options
    hpc_platform = options["simulation"]["hpc_system"]
    gamera_grid_type = options["simulation"]["gamera_grid_type"]
    od = option_descriptions["pbs"][hpc_platform]
    od["select"]["default"] = od["select"]["default"][gamera_grid_type]
    od["num_helpers"]["default"] = od["num_helpers"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # Compute the number of nodes in the second chunk.
    # Should be 1 (for voltron) + num_helpers.
    num_helpers = int(o["num_helpers"])
    select2 = 1 + num_helpers
    o["select2"] = str(select2)

    # ------------------------------------------------------------------------

    # GAMERA options
    options["gamera"] = {}

    # <sim> options
    o = options["gamera"]["sim"] = {}
    od = option_descriptions["gamera"]["sim"]
    od["H5Grid"]["default"] = f"lfm{options['simulation']['gamera_grid_type']}.h5"
    od["runid"]["default"] = options["simulation"]["job_name"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <floors> options
    o = options["gamera"]["floors"] = {}
    od = option_descriptions["gamera"]["floors"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <timestep> options
    o = options["gamera"]["timestep"] = {}
    od = option_descriptions["gamera"]["timestep"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <restart> options
    # NOTE: Update this later so that restart parameters are only
    # prompted for when doRes is "T".
    o = options["gamera"]["restart"] = {}
    od = option_descriptions["gamera"]["restart"]
    od["resID"]["default"] = options["simulation"]["job_name"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <physics> options
    o = options["gamera"]["physics"] = {}
    od = option_descriptions["gamera"]["physics"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <ring> options
    o = options["gamera"]["ring"] = {}
    od = option_descriptions["gamera"]["ring"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <wind> options
    o = options["gamera"]["wind"] = {}
    od = option_descriptions["gamera"]["wind"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <source> options
    o = options["gamera"]["source"] = {}
    od = option_descriptions["gamera"]["source"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <iPdir> options
    o = options["gamera"]["iPdir"] = {}
    od = option_descriptions["gamera"]["iPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <jPdir> options
    o = options["gamera"]["jPdir"] = {}
    od = option_descriptions["gamera"]["jPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <kPdir> options
    o = options["gamera"]["kPdir"] = {}
    od = option_descriptions["gamera"]["kPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <coupling> options
    o = options["gamera"]["coupling"] = {}
    od = option_descriptions["gamera"]["coupling"]
    od["blockHalo"]["default"] = od["blockHalo"]["default"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # ------------------------------------------------------------------------

    # VOLTRON options
    options["voltron"] = {}

    # <time> options
    o = options["voltron"]["time"] = {}
    od = option_descriptions["voltron"]["time"]
    od["tFin"]["default"] = simulation_duration
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <spinup> options
    o = options["voltron"]["spinup"] = {}
    od = option_descriptions["voltron"]["spinup"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <output> options
    o = options["voltron"]["output"] = {}
    od = option_descriptions["voltron"]["output"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <coupling> options
    o = options["voltron"]["coupling"] = {}
    od = option_descriptions["voltron"]["coupling"]
    od["doAsyncCoupling"]["default"] = od["doAsyncCoupling"]["default"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <restart> options
    o = options["voltron"]["restart"] = {}
    od = option_descriptions["voltron"]["restart"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <imag> options
    o = options["voltron"]["imag"] = {}
    od = option_descriptions["voltron"]["imag"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <helpers> options
    o = options["voltron"]["helpers"] = {}
    od = option_descriptions["voltron"]["helpers"]
    od["numHelpers"]["default"] = num_helpers
    if num_helpers == 0:
        od["useHelpers"]["default"] = "F"
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # ------------------------------------------------------------------------

    # CHIMP options
    options["chimp"] = {}

    # <units> options
    o = options["chimp"]["units"] = {}
    od = option_descriptions["chimp"]["units"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <fields> options
    o = options["chimp"]["fields"] = {}
    od = option_descriptions["chimp"]["fields"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <domain> options
    o = options["chimp"]["domain"] = {}
    od = option_descriptions["chimp"]["domain"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <tracer> options
    o = options["chimp"]["tracer"] = {}
    od = option_descriptions["chimp"]["tracer"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # ------------------------------------------------------------------------

    # REMIX options
    options["remix"] = {}

    # <conductance> options
    o = options["remix"]["conductance"] = {}
    od = option_descriptions["remix"]["conductance"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <precipitation> options
    o = options["remix"]["precipitation"] = {}
    od = option_descriptions["remix"]["precipitation"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # ------------------------------------------------------------------------

    # RAIJU options
    options["raiju"] = {}

    # <output> options
    options["raiju"]["output"] = {}
    o = options["raiju"]["output"]
    od = option_descriptions["raiju"]["output"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <grid> options
    options["raiju"]["grid"] = {}
    o = options["raiju"]["grid"]
    od = option_descriptions["raiju"]["grid"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <domain> options
    options["raiju"]["domain"] = {}
    o = options["raiju"]["domain"]
    od = option_descriptions["raiju"]["domain"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <config> options
    options["raiju"]["config"] = {}
    o = options["raiju"]["config"]
    od = option_descriptions["raiju"]["config"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <plasmasphere> options
    options["raiju"]["plasmasphere"] = {}
    o = options["raiju"]["plasmasphere"]
    od = option_descriptions["raiju"]["plasmasphere"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <losses> options
    options["raiju"]["losses"] = {}
    o = options["raiju"]["losses"]
    od = option_descriptions["raiju"]["losses"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <cpl> options
    options["raiju"]["cpl"] = {}
    o = options["raiju"]["cpl"]
    od = option_descriptions["raiju"]["cpl"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <fluidIn1> options
    options["raiju"]["fluidIn1"] = {}
    o = options["raiju"]["fluidIn1"]
    od = option_descriptions["raiju"]["fluidIn1"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # Return the options dictionary.
    return options


def validate_grid_file(path : str, grid_type : str):
    """Validate a LFM grid file.

    Validate a LFM grid file. Check that the grid sizes in each dimension
    are correct for the selected resolution.

    Parameters
    ----------
    path : str
        Path to LFM grid file to validate.
    grid_type : str
        Single character grid type specifier (D)ouble, (Q)uad, (O)ct, (H)ex

    Returns
    -------
    is_valid : bool
        True if grid file is valid, else False.

    Raises
    ------
    None
    """
    # Specify the expected grid cell counts in each dimension for each grid
    # type.
    grid_shapes = {
        'D': (73, 57, 57),
        'Q': (137, 105, 105),
        'O': (265, 201, 201),
        'H': (521, 393, 393),
    }

    # Open the file and read the grid sizes.
    is_valid = False
    with h5py.File(path, 'r') as hf:
        if (hf['X'].shape == grid_shapes[grid_type] and
            hf['Y'].shape == grid_shapes[grid_type] and
            hf['Z'].shape == grid_shapes[grid_type]):
            is_valid = True

    # Return the result of the check.
    return is_valid


def run_preprocessing_steps(options: dict, args: dict):
    """Execute any preprocessing steps required for the run.

    Execute any preprocessing steps required for the run.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.
    args : dict
        Dictionary of command-line options and options from engage.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Local convenience variables
    debug = args["debug"]

    # Check if the grid file exists and is valid. Create if needed.
    grid_file = f'lfm{options["simulation"]["gamera_grid_type"]}.h5'
    if os.path.exists(grid_file):
        if validate_grid_file(grid_file,
                              options['simulation']['gamera_grid_type']):
            print(f'Using existing grid file {grid_file}.')
        else:
            raise TypeError(f'Invalid grid file {grid_file} found, aborting.')
    else:
        # Create the LFM grid file. Assumes kaipy is installed.
        cmd = 'genLFM'
        args = [cmd, '-gid', options['simulation']['gamera_grid_type'],
                '-Rin', options['simulation']['gamera_grid_inner_radius'],
                '-Rout', options['simulation']['gamera_grid_outer_radius']]
        if debug:
            print(f"cmd = {args}")
        subprocess.run(args, check=True)

    # If needed, create the solar wind file by fetching data from CDAWeb.
    # NOTE: Assumes kaipy is installed.
    if options["simulation"]["bcwind_available"].upper() == "N":
        cmd = "cda2wind"
        args = [cmd, "-t0", options["simulation"]["start_date"], "-t1",
                options["simulation"]["stop_date"], "-interp", "-bx",
                "-f107", "100", "-kp", "3"]
        subprocess.run(args, check=True)

    # Create the RAIJU configuration file.
    # NOTE: Assumes kaipy is installed.
    cmd = "genRAIJU"
    args = [cmd]
    subprocess.run(args, check=True)


def create_ini_files(options: dict, args: dict):
    """Create the MAGE .ini files from a template.

    Create the MAGE .ini files from a template.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.
    args : dict
        Dictionary of command-line options and options from engage

    Returns
    -------
    ini_files : list of str
        Paths to the .ini files for the gamera run.

    Raises
    ------
    None
    """
    # Read and create the template.
    with open(INI_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)

    # Set default value for padding to tFin for coupling.
    tfin_padding = 0.0
    # Engage modifications to parameters.
    # If TIEGCM coupling is specified, calculate warmup segments based 
    # on gr_warm_up_time and segment_duration.
    # If gr_warm_up_time is an exact multiple of segment_duration, 
    # use segment_duration for each warmup segment.
    # If gr_warm_up_time is less than segment_duration, use 
    # gr_warm_up_time as the duration for a single warmup segment.
    # If gr_warm_up_time is greater than segment_duration but not an 
    # exact multiple, use segment_duration and round up the number of segments.
    if "coupling" in args:     
        coupling = args["coupling"]
        gr_warm_up_time = float(coupling["gr_warm_up_time"])
        segment_duration = float(options["simulation"]["segment_duration"])
        simulation_duration = float(options["voltron"]["time"]["tFin"])
        if gr_warm_up_time % segment_duration == 0:
            warmup_segment_duration = segment_duration
            i_last_warmup_ini = int(gr_warm_up_time/warmup_segment_duration)
        elif gr_warm_up_time < segment_duration:
            warmup_segment_duration = gr_warm_up_time
            i_last_warmup_ini = int(gr_warm_up_time/warmup_segment_duration)
        elif gr_warm_up_time > segment_duration:
            warmup_segment_duration = segment_duration
            i_last_warmup_ini = int(math.ceil(gr_warm_up_time/segment_duration))
        # Add padding to tFin for coupling.
        if coupling["tfin_delta"] == "T":
            tfin_coupling_padding = float(options["voltron"]["coupling"]["dtCouple"]) - 1
        else:
            tfin_coupling_padding = 0.0
        # Set doGCM to value from engage.
        if "doGCM" in coupling:
            doGCM = coupling["doGCM"]
        else:   
            doGCM = coupling["doGCM"]


    # Initialize the list of file paths.
    ini_files = []

    # Create the .ini files for each PBS job.
    if options["simulation"]["use_segments"].upper() == "Y":

        # Extract the end time (in seconds) for the simulation.
        tFin = float(options["voltron"]["time"]["tFin"])

        # Create an .ini file for the spinup segment, if requested.
        if options["voltron"]["spinup"]["doSpin"] == "T":
            opt = copy.deepcopy(options)  # Need a copy of options
            runid = opt["simulation"]["job_name"]
            segment_id = f"{runid}-SPINUP"
            opt["simulation"]["segment_id"] = segment_id
            # Just perform spinup in first segment. Use an end time of 1.0
            # to ensure the restart file from the end of spinup is properly
            # generated.
            tFin_segment = 1.0
            opt["voltron"]["time"]["tFin"] = str(tFin_segment)
            opt["gamera"]["restart"]["doRes"] = "F"
            ini_content = template.render(opt)
            ini_file = os.path.join(opt["pbs"]["run_directory"],
                                    f"{segment_id}.ini")
            ini_files.append(ini_file)
            with open(ini_file, "w", encoding="utf-8") as f:
                f.write(ini_content)

        # Create an .ini file for the warmup segment, if requested.
        # NOTE: This is a special case for the GTR runs. The
        # gr_warm_up_time is used to determine the end time of the
        # warmup segment
        num_warmup_segments = 0
        if "coupling" in args:
            tFin_warmup = float(coupling["gr_warm_up_time"]) + 1.0
            for job in range(1, i_last_warmup_ini + 1):
                opt = copy.deepcopy(options)
                runid = opt["simulation"]["job_name"]
                # NOTE: This naming scheme supports a maximum of 99 segments.
                segment_id = f"{runid}-WARMUP-{job:02d}"
                opt["simulation"]["segment_id"] = segment_id
                opt["gamera"]["restart"]["doRes"] = "T"
                dT = warmup_segment_duration #float(options["simulation"]["segment_duration"])
                # Add 1 to ensure last restart file created
                tFin_segment = job*dT + 1.0
                if tFin_segment > tFin_warmup:
                    tFin_segment = tFin_warmup
                opt["voltron"]["time"]["tFin"] = str(tFin_segment)
                dtRes = float(options["voltron"]["restart"]["dtRes"])
                if job > 1:
                    nRes = int(((tFin_segment - 1) - dT )/dtRes) 
                else:
                    nRes = 0
                opt["gamera"]["restart"]["nRes"] = str(nRes)
                ini_content = template.render(opt)
                ini_file = os.path.join(opt["pbs"]["run_directory"],
                                        f"{segment_id}.ini")
                ini_files.append(ini_file)
                with open(ini_file, "w", encoding="utf-8") as f:
                    f.write(ini_content)
            num_warmup_segments = i_last_warmup_ini 
        # Create an .ini file for each simulation segment. Files for each
        # segment will be numbered starting with 1.
        if "coupling" in args:
            num_segments = math.ceil((simulation_duration - num_warmup_segments*warmup_segment_duration)/segment_duration)
        else:
            num_segments = int(options["pbs"]["num_segments"])
        for job in range(1, num_segments + 1):
            opt = copy.deepcopy(options)  # Need a copy of options
            runid = opt["simulation"]["job_name"]
            # NOTE: This naming scheme supports a maximum of 99 segments.
            segment_id = f"{runid}-{job:02d}"
            opt["simulation"]["segment_id"] = segment_id
            opt["gamera"]["restart"]["doRes"] = "T"
            dT = float(options["simulation"]["segment_duration"])
            dtRes = float(options["voltron"]["restart"]["dtRes"])
            # tfin_delta is the warmup time to add to tFin for coupling.
            if "coupling" in args:
                # tFin for coupling is different.
                tfin_delta = float(coupling["gr_warm_up_time"])
            else:
                tfin_delta = 0.0
            # Add 1 to ensure last restart file created    
            tFin_segment = job*dT + tfin_delta + 1.0
            nRes = int(((tFin_segment - 1) - dT )/dtRes) 
            opt["gamera"]["restart"]["nRes"] = str(nRes)   
            # Engage modifications to parameters in coupled segment.
            if "coupling" in args:
                opt["voltron"]["coupling"]["doGCM"] = doGCM
                # tFin padding different for last segment.
                if job == num_segments:
                    tfin_padding = -1.0
                else:
                    # Subtract 1 from tFin padding for coupling beacuse to offset the +1.0 for restart file done above.
                    tfin_padding = tfin_coupling_padding - 1.0
            opt["voltron"]["time"]["tFin"] = str(tFin_segment + tfin_padding)
            #print(f'Creating job {job} with tFin_seg = {opt["voltron"]["time"]["tFin"]}')
            ini_content = template.render(opt)
            ini_file = os.path.join(opt["pbs"]["run_directory"],
                                    f"{segment_id}.ini")
            ini_files.append(ini_file)
            with open(ini_file, "w", encoding="utf-8") as f:
                f.write(ini_content)

    else:
        # Use a single job segment.
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        segment_id = runid
        opt["simulation"]["segment_id"] = segment_id
        ini_content = template.render(opt)
        ini_file = os.path.join(opt["pbs"]["run_directory"],
                                f"{segment_id}.ini")
        ini_files.append(ini_file)
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)

    # If a warmup period was used, rename the .ini files for the segments
    # which cover the warmup period. Then rename the remaining files to
    # account for this change.
    '''
    if "coupling" in args:

        # Compute the number of warmup segments.
        coupling = args["coupling"]
        gr_warm_up_time = float(coupling["gr_warm_up_time"])
        segment_duration = float(options["simulation"]["segment_duration"])
        i_last_warmup_ini = int(gr_warm_up_time/segment_duration)

        # Rename the warmup segments.
        for i in range(1, i_last_warmup_ini + 1):
            old_name = ini_files[i]
            new_name = f"{runid}-WARMUP-{i:02d}.ini"
            os.rename(old_name, new_name)
            ini_files[i] = new_name

        # Rename the remaining segments.
        for i in range(i_last_warmup_ini + 1, len(ini_files)):
            old_name = ini_files[i]
            new_name = f"{runid}-{i - i_last_warmup_ini:02d}.ini"
            os.rename(old_name, new_name)
            ini_files[i] = new_name
    '''
    # Return the paths to the .ini files.
    return ini_files


def convert_ini_to_xml(ini_files: list, args: dict):
    """Convert the .ini files to XML.

    Convert the .ini files describing the run to XML files. The intermediate
    .ini files are then deleted.

    Parameters
    ----------
    ini_files : list of str
        Paths to the .ini files to convert.
    args : dict
        Dictionary of command-line options and options from engage

    Returns
    -------
    xml_files : str
        Paths to the XML files.

    Raises
    ------
    None
    """
    # Local convenience variables
    debug = args["debug"]

    # Convert each .ini file to an .xml file.
    xml_files = []
    for ini_file in ini_files:

        # Put the XML file in the same directory as the .ini file.
        xml_file = ini_file.replace(".ini", ".xml")

        # Convert the .ini file to .xml.
        # NOTE: assumes kaipy is installed.
        cmd = "XMLGenerator"
        args = [cmd, ini_file, xml_file]
        subprocess.run(args, check=True)

        # Add this file to the list of XML files.
        xml_files.append(xml_file)

        # Remove the .ini file.
        os.remove(ini_file)

    # Return the paths to the XML files.
    return xml_files


def create_pbs_scripts(xml_files: list, options: dict, args: dict):
    """Create the PBS job scripts for the run.

    Create the PBS job scripts from a template.

    Parameters
    ----------
    xml_files : str
        Paths to the XML files.
    options : dict
        Dictionary of program options, each entry maps str to str.
    args : dict
        Dictionary of command-line options and options from engage

    Returns
    -------
    pbs_scripts : list of str
        Paths to PBS job scripts.
    submit_all_jobs_script : str
        Path to script which submits all PBS jobs.
    warmup_pbs_scripts : list of str
        List of PBS job scripts which encompass the MAGE warmup period used
        when coupling with TIEGCM.

    Raises
    ------
    TypeError:
        For a non-integral of nodes requested
    """
    # Local convenience variables.
    debug = args["debug"]

    # Compute the number of nodes to request based on the MPI decomposition
    # and the MPI ranks per node.
    ni = int(options["gamera"]["iPdir"]["N"])
    nj = int(options["gamera"]["jPdir"]["N"])
    nk = int(options["gamera"]["kPdir"]["N"])
    ranks_per_node = int(options["pbs"]["mpiprocs"])
    select_nodes = ni*nj*nk/ranks_per_node
    if int(select_nodes) != select_nodes:
        raise TypeError(f"Requested non-integral node count ({select_nodes})!")
    options["pbs"]["select"] = str(int(select_nodes))

    # Read the template.
    with open(PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)

    # Create a PBS script for each segment.
    pbs_scripts = []
    # for job in range(int(options["pbs"]["num_segments"])):
    for xml_file in xml_files:
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        filename = os.path.split(xml_file)[-1]
        if debug:
            print(f"filename = {filename}")
        segment_id = filename.replace(".xml", "")
        if debug:
            print(f"segment_id = {segment_id}")
        opt["simulation"]["segment_id"] = segment_id
        pbs_content = template.render(opt)
        pbs_script = os.path.join(
            opt["pbs"]["run_directory"], f"{segment_id}.pbs"
        )
        pbs_scripts.append(pbs_script)
        with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)

    # Create a single script which will submit all of the PBS jobs in order.
    # Each job only runs if all previouos jobs run successfully.
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

    # Make a list of the scripts which cover the warmup period.
    # NOTE: Assumes the warmup period is an integral multiple of the segment
    # length.
    # NOTE: Assumes job segments are in use, with only 1 spinup segment.
    warmup_pbs_scripts = []
    spinup_pbs_scripts = []
    if "coupling" in args:
        coupling = args["coupling"]
        gr_warm_up_time = float(coupling["gr_warm_up_time"])
        segment_duration = float(options["simulation"]["segment_duration"])
        simulation_duration = float(options["voltron"]["time"]["tFin"])
        if gr_warm_up_time % segment_duration == 0:
            warmup_segment_duration = segment_duration
            i_last_warmup_ini = int(gr_warm_up_time/warmup_segment_duration)
        elif gr_warm_up_time < segment_duration:
            warmup_segment_duration = gr_warm_up_time
            i_last_warmup_ini = int(gr_warm_up_time/warmup_segment_duration)
        elif gr_warm_up_time > segment_duration:
            warmup_segment_duration = segment_duration
            i_last_warmup_ini = int(math.ceil(gr_warm_up_time/segment_duration))
        num_warmup_segments = i_last_warmup_ini
        #i_last_warmup_pbs_script = int(gr_warm_up_time/segment_duration)
        spinup_pbs_scripts.append(pbs_scripts[0]) # Spinup script is first
        warmup_pbs_scripts = pbs_scripts[1:num_warmup_segments + 1] # Warmup scripts
    # Return the paths to the PBS scripts.
    return pbs_scripts, submit_all_jobs_script,spinup_pbs_scripts, warmup_pbs_scripts


def fetch_select_line(path: str):
    """Extract the select line from a PBS script.

    Extract the select line from a PBS script.

    Parameters
    ----------
    path : str
        Path to PBS script.

    Returns
    -------
    select_line : str
        The first '#PBS select' line from the PBS script.

    Raises
    ------
    None
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    select_line = None
    for line in lines:
        if line.startswith("#PBS -l select="):
            select_line = line.rstrip()
            break
    return select_line


# ----------------------------------------------------------------------------

# makeitso() is the primary entry point to this module. It will be called
# by main() when this module is run on the command line, or explicitly from a
# calling function after importing the makeitso module.


def makeitso(args: dict = None):
    """Main program code for makeitso.

    This is the main program code for makeitso. This function can be called
    from other python code.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options and equivalent options passed from
        the calling function, and variables set by engage for makeitso.

    Returns
    -------
    select_line : str
        The '#PBS -l select' line from the first PBS file created for this run.
    mpiexec_command : str
        The mpiexec command used to run MAGE software in the PBS files.

    Raises
    ------
    None
    """
    # Use defaults for unspecified arguments. Merge in additional settings
    # which are passed in by the caller.
    local_args = copy.deepcopy(args_default)
    if args is not None:
        local_args.update(args)
    args = local_args

    # Local convenience variables
    if args["debug"]:
        print(f"args = {args}")
    clobber = args["clobber"]
    debug = args["debug"]
    # mode = args["mode"]
    options_path = args["options_path"]
    verbose = args["verbose"]

    # ------------------------------------------------------------------------

    # Read the default option descriptions file.
    option_descriptions = load_option_descriptions()
    if debug:
        print(f"option_descriptions = {option_descriptions}")

    # Update the option descriptions with information passed in from the
    # calling code.
    option_descriptions = update_option_descriptions(option_descriptions, args)
    if debug:
        print(f"option_descriptions = {option_descriptions}")

    # ------------------------------------------------------------------------

    # Fetch the run options.
    if options_path:
        # Read the run options from a JSON file.
        if verbose:
            print(f"Reading run options from {options_path}.")
        with open(options_path, "r", encoding="utf-8") as f:
            options = json.load(f)
        options = update_loaded_options(options)
    else:
        # Prompt the user for the run options.
        options = prompt_user_for_run_options(option_descriptions, args)
    if debug:
        print(f"options = {options}")

    # Move to the run directory.
    os.chdir(options["pbs"]["run_directory"])

    # Save the options dictionary as a JSON file in the current directory.
    path = f"{options['simulation']['job_name']}.json"
    if os.path.exists(path):
        if not clobber:
            raise FileExistsError(f"Options file {path} exists!")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(options, args)

    # Create the .ini file(s).
    if verbose:
        print("Creating .ini file(s) for run.")
    ini_files = create_ini_files(options, args)
    if debug:
        print(f"ini_files = {ini_files}")

    # Convert the .ini file(s) to .xml files(s).
    if verbose:
        print("Converting .ini file(s) to .xml file(s).")
    xml_files = convert_ini_to_xml(ini_files, args)
    if debug:
        print(f"xml_files = {xml_files}")

    # Create the PBS job script(s).
    if verbose:
        print("Creating PBS job script(s) for run.")
    pbs_scripts, all_jobs_script, spinup_pbs_scripts, warmup_pbs_scripts = create_pbs_scripts(
        xml_files, options, args)
    if verbose:
        print(f"The PBS job scripts {pbs_scripts} are ready.")
    print(f"The PBS scripts {pbs_scripts} have been created, each with a "
          "corresponding XML file. To submit the jobs with the proper "
          "dependency (to ensure each segment runs in order), please run the "
          f"script {all_jobs_script} like this:\n"
          f"bash {all_jobs_script}")

    # Return the options dict used to create the PBS scripts, and the list
    # of PBS scripts which constitute the warmup period.
    return options, spinup_pbs_scripts, warmup_pbs_scripts


def main():
    """Main program code for the command-line version of makeitso.

    This is the main program code for the command-line version of makeitso.
    It processes command-line options, then calls the makeitso() function.

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

    # ------------------------------------------------------------------------

    # Call the main program logic. Note that the Namespace object (args)
    # returned from the option parser is converted to a dict using vars().
    makeitso(vars(args))


if __name__ == "__main__":
    main()
