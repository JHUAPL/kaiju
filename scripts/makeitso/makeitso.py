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

# Import project modules.


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
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
)

# Path to template .ini file.
INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.ini")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")


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


def fetch_bcwind_time_range(bcwind_path):
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
        # "YYYY-MM-DDTHH:MM:SS" format.
        start_date = start_date.replace(" ", "T")
        stop_date = stop_date.replace(" ", "T")
        # </HACK>
    return start_date, stop_date


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

    # Ask the user if a boundary condition file is available. If not, offer to
    # generate one from the start and end date.
    for on in ["bcwind_available"]:
        o[on] = get_run_option(on, od[on], mode)
    if o["bcwind_available"] == "Y":
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
    for on in ["gamera_grid_type", "gamera_grid_inner_radius", 
               "gamera_grid_outer_radius","hpc_system"]:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # PBS options
    options["pbs"] = {}
    o = options["pbs"]

    # Common (HPC platform-independent) options
    od = option_descriptions["pbs"]["_common"]
    od["account_name"]["default"] = os.getlogin()
    od["kaiju_install_directory"]["default"] = KAIJUHOME
    od["kaiju_build_directory"]["default"] = os.path.join(KAIJUHOME, "build_mpi")
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

    #-------------------------------------------------------------------------

    # GAMERA options
    options["gamera"] = {}

    # <sim> options
    options["gamera"]["sim"] = {}
    o = options["gamera"]["sim"]
    od = option_descriptions["gamera"]["sim"]
    od["H5Grid"]["default"] = f"lfm{options['simulation']['gamera_grid_type']}.h5"
    od["runid"]["default"] = options["simulation"]["job_name"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <floors> options
    options["gamera"]["floors"] = {}
    o = options["gamera"]["floors"]
    od = option_descriptions["gamera"]["floors"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <timestep> options
    options["gamera"]["timestep"] = {}
    o = options["gamera"]["timestep"]
    od = option_descriptions["gamera"]["timestep"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <restart> options
    # NOTE: Update this later so that restart parameters are only
    # prompted for when doRes is "T".
    options["gamera"]["restart"] = {}
    o = options["gamera"]["restart"]
    od = option_descriptions["gamera"]["restart"]
    od["resID"]["default"] = options["simulation"]["job_name"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <physics> options
    options["gamera"]["physics"] = {}
    o = options["gamera"]["physics"]
    od = option_descriptions["gamera"]["physics"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <ring> options
    options["gamera"]["ring"] = {}
    o = options["gamera"]["ring"]
    od = option_descriptions["gamera"]["ring"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <wind> options
    options["gamera"]["wind"] = {}
    o = options["gamera"]["wind"]
    od = option_descriptions["gamera"]["wind"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <source> options
    options["gamera"]["source"] = {}
    o = options["gamera"]["source"]
    od = option_descriptions["gamera"]["source"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <iPdir> options
    options["gamera"]["iPdir"] = {}
    o = options["gamera"]["iPdir"]
    od = option_descriptions["gamera"]["iPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <jPdir> options
    options["gamera"]["jPdir"] = {}
    o = options["gamera"]["jPdir"]
    od = option_descriptions["gamera"]["jPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <kPdir> options
    options["gamera"]["kPdir"] = {}
    o = options["gamera"]["kPdir"]
    od = option_descriptions["gamera"]["kPdir"]
    od["N"]["default"] = od["N"]["default"][gamera_grid_type]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <coupling> options
    options["gamera"]["coupling"] = {}
    o = options["gamera"]["coupling"]
    od = option_descriptions["gamera"]["coupling"]
    od["blockHalo"]["default"] = od["blockHalo"]["default"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # VOLTRON options
    options["voltron"] = {}

    # <time> options
    options["voltron"]["time"] = {}
    o = options["voltron"]["time"]
    od = option_descriptions["voltron"]["time"]
    od["tFin"]["default"] = simulation_duration
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <spinup> options
    options["voltron"]["spinup"] = {}
    o = options["voltron"]["spinup"]
    od = option_descriptions["voltron"]["spinup"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <output> options
    options["voltron"]["output"] = {}
    o = options["voltron"]["output"]
    od = option_descriptions["voltron"]["output"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <coupling> options
    options["voltron"]["coupling"] = {}
    o = options["voltron"]["coupling"]
    od = option_descriptions["voltron"]["coupling"]
    od["doAsyncCoupling"]["default"] = od["doAsyncCoupling"]["default"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <restart> options
    options["voltron"]["restart"] = {}
    o = options["voltron"]["restart"]
    od = option_descriptions["voltron"]["restart"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <imag> options
    options["voltron"]["imag"] = {}
    o = options["voltron"]["imag"]
    od = option_descriptions["voltron"]["imag"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <helpers> options
    options["voltron"]["helpers"] = {}
    o = options["voltron"]["helpers"]
    od = option_descriptions["voltron"]["helpers"]
    od["numHelpers"]["default"] = num_helpers
    if num_helpers == 0:
        od["useHelpers"]["default"] = "F"
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # CHIMP options
    options["chimp"] = {}

    # <units> options
    options["chimp"]["units"] = {}
    o = options["chimp"]["units"]
    od = option_descriptions["chimp"]["units"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <fields> options
    options["chimp"]["fields"] = {}
    o = options["chimp"]["fields"]
    od = option_descriptions["chimp"]["fields"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <domain> options
    options["chimp"]["domain"] = {}
    o = options["chimp"]["domain"]
    od = option_descriptions["chimp"]["domain"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <tracer> options
    options["chimp"]["tracer"] = {}
    o = options["chimp"]["tracer"]
    od = option_descriptions["chimp"]["tracer"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # REMIX options
    options["remix"] = {}

    # <conductance> options
    options["remix"]["conductance"] = {}
    o = options["remix"]["conductance"]
    od = option_descriptions["remix"]["conductance"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <precipitation> options
    options["remix"]["precipitation"] = {}
    o = options["remix"]["precipitation"]
    od = option_descriptions["remix"]["precipitation"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # RCM options
    options["rcm"] = {}

    # <conductance> options
    options["rcm"]["rcmdomain"] = {}
    o = options["rcm"]["rcmdomain"]
    od = option_descriptions["rcm"]["rcmdomain"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <ellipse> options
    # Only use if domType == "ELLIPSE"?
    options["rcm"]["ellipse"] = {}
    o = options["rcm"]["ellipse"]
    od = option_descriptions["rcm"]["ellipse"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <grid> options
    options["rcm"]["grid"] = {}
    o = options["rcm"]["grid"]
    od = option_descriptions["rcm"]["grid"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <plasmasphere> options
    options["rcm"]["plasmasphere"] = {}
    o = options["rcm"]["plasmasphere"]
    od = option_descriptions["rcm"]["plasmasphere"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------

    # Return the options dictionary.
    return options


def run_preprocessing_steps(options):
    """Execute any preprocessing steps required for the run.

    Execute any preprocessing steps required for the run.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Create the LFM grid file.
    # NOTE: Assumes genLFM.py is in PATH.
    cmd = "genLFM.py"
    args = [cmd, "-gid", options["simulation"]["gamera_grid_type"],
            '-Rin', options["simulation"]["gamera_grid_inner_radius"],
            '-Rout', options["simulation"]["gamera_grid_outer_radius"]]
    subprocess.run(args, check=True)

    # If needed, create the solar wind file by fetching data from CDAWeb.
    # NOTE: Assumes cda2wind.py is in PATH.
    if options["simulation"]["bcwind_available"] == "N":
        cmd = "cda2wind.py"
        args = [cmd, "-t0", options["simulation"]["start_date"], "-t1",
                options["simulation"]["stop_date"], "-interp", "-bx",
                "-f107", "100", "-kp", "3"]
        subprocess.run(args, check=True)

    # Create the RCM configuration file.
    # NOTE: Assumes genRCM.py is in PATH.
    cmd = "genRCM.py"
    args = [cmd]
    subprocess.run(args, check=True)


def create_ini_files(options):
    """Create the MAGE .ini files from a template.

    Create the MAGE .ini files from a template.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    ini_files : list of str
        Paths to the .ini files for the gamera run.

    Raises
    ------
    None
    """
    # Read and create the template.
    template_file = INI_TEMPLATE
    with open(template_file, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)

    # Initialize the list of file paths.
    ini_files = []

    # Create the job scripts.
    if int(options["pbs"]["num_segments"]) > 1:

        # Create an .ini file for the spinup segment.
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        job = 0
        segment_id = f"{runid}-{job:02d}"
        opt["simulation"]["segment_id"] = segment_id
        tFin = float(opt["voltron"]["time"]["tFin"])
        dT = float(options["simulation"]["segment_duration"])
        tFin_segment = 1.0  # Just perform spinup in first segment
        opt["voltron"]["time"]["tFin"] = str(tFin_segment)
        ini_content = template.render(opt)
        ini_file = os.path.join(
            opt["pbs"]["run_directory"], f"{opt['simulation']['segment_id']}.ini"
        )
        ini_files.append(ini_file)
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)

        # Create an .ini file for each simulation segment.
        for job in range(1, int(options["pbs"]["num_segments"])):
            opt = copy.deepcopy(options)  # Need a copy of options
            runid = opt["simulation"]["job_name"]
            segment_id = f"{runid}-{job:02d}"
            opt["simulation"]["segment_id"] = segment_id
            opt["gamera"]["restart"]["doRes"] = "T"
            tFin = float(opt["voltron"]["time"]["tFin"])
            dT = float(options["simulation"]["segment_duration"])
            tFin_segment = job*dT + 1  # Add 1 to ensure last restart file is created
            if tFin_segment > tFin:    # Last segment may be shorter than the others.
                tFin_segment = tFin + 1
            opt["voltron"]["time"]["tFin"] = str(tFin_segment)
            ini_content = template.render(opt)
            ini_file = os.path.join(
                opt["pbs"]["run_directory"], f"{opt['simulation']['segment_id']}.ini"
            )
            ini_files.append(ini_file)
            with open(ini_file, "w", encoding="utf-8") as f:
                f.write(ini_content)

    else:
        # Use a single job segment.
        job = 0
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        segment_id = f"{runid}-{job:02d}"
        opt["simulation"]["segment_id"] = segment_id
        ini_content = template.render(opt)
        ini_file = os.path.join(
            opt["pbs"]["run_directory"], f"{opt['simulation']['segment_id']}.ini"
        )
        ini_files.append(ini_file)
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)

    # Return the paths to the .ini files.
    return ini_files


def convert_ini_to_xml(ini_files):
    """Convert the .ini files to XML.

    Convert the .ini files describing the run to XML files. The intermediate
    .ini files are then deleted.

    Parameters
    ----------
    ini_files : list of str
        Paths to the .ini files to convert.

    Returns
    -------
    xml_files : str
        Paths to the XML files.

    Raises
    ------
    None
    """
    # Convert each .ini file to an .xml file.
    xml_files = []
    for ini_file in ini_files:

        # Put the XML file in the same directory as the .ini file.
        xml_file = ini_file.replace(".ini", ".xml")

        # Convert the .ini file to .xml.
        # NOTE: assumes XMLGenerator.py is in PATH.
        cmd = "XMLGenerator.py"
        args = [cmd, ini_file, xml_file]
        subprocess.run(args, check=True)

        # Add this file to the list of XML files.
        xml_files.append(xml_file)

        # Remove the .ini file.
        os.remove(ini_file)

    # Return the paths to the XML files.
    return xml_files


def create_pbs_scripts(options):
    """Create the PBS job scripts for the run.

    Create the PBS job scripts from a template.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

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
        cmd = f"echo $job_id\n"
        f.write(cmd)
        for s in pbs_scripts[1:]:
            cmd = "old_job_id=$job_id\n"
            f.write(cmd)
            cmd = f"job_id=`qsub -W depend=afterok:$old_job_id {s}`\n"
            f.write(cmd)
            cmd = f"echo $job_id\n"
            f.write(cmd)

    # Return the paths to the PBS scripts.
    return pbs_scripts, submit_all_jobs_script


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
    run_preprocessing_steps(options)

    # Create the .ini file(s).
    if verbose:
        print("Creating .ini file(s) for run.")
    ini_files = create_ini_files(options)
    if debug:
        print(f"ini_files = {ini_files}")

    # Convert the .ini file(s) to .xml files(s).
    if verbose:
        print("Converting .ini file(s) to .xml file(s).")
    xml_files = convert_ini_to_xml(ini_files)
    if debug:
        print(f"xml_files = {xml_files}")

    # Create the PBS job script(s).
    if verbose:
        print("Creating PBS job script(s) for run.")
    pbs_scripts, all_jobs_script = create_pbs_scripts(options)
    if verbose:
        print(f"The PBS job scripts {pbs_scripts} are ready.")
    print(f"The PBS scripts {pbs_scripts} have been created, each with a "
          "corresponding XML file. To submit the jobs with the proper "
          "dependency (to ensure each segment runs in order), please run the "
          f"script {all_jobs_script} like this:\n"
          f"bash {all_jobs_script}")

if __name__ == "__main__":
    """Begin main program."""
    main()
