#!/usr/bin/env python


"""makeitso for the GAMERA heliosphere software.

This script performs all of the steps needed to prepare to run a GAMERA MHD
heliosphere simulation run. By default, this script is interactive - the user
is prompted for each decision that must be made to prepare for the run, based
on the current "--mode" setting.

The modes are:

"BASIC" (the default) - the user is prompted to set only a small subset of GAMERA
parameters. All "INTERMEDIATE"- and "EXPERT"-mode parameters are automatically
set to default values.

"INTERMEDIATE" - The user is prompted for "BASIC" and "INTERMEDIATE"
parameters, with "EXPERT" parameters set to defaults.

"EXPERT" - The user is prompted for *all* adjustable parameters.

"COMPUTED" - The value of the option is computed internally. It is not
available to the user as an adjustable parameter.
"""


# Import standard modules.
import argparse
import copy
import datetime
import json
import os
import subprocess
import sys

# Import 3rd-party modules.
from astropy.io import fits
from jinja2 import Template

# Import project modules.
from kaipy.kdefs import JD2MJD, Tsolar_synodic
from kaipy.kaiTools import MJD2UT

# Program constants

# Program description
DESCRIPTION = "Interactive script to prepare a GAMERA heliosphere run"

# Indent level for JSON output
JSON_INDENT = 4

# Path to current kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Path to directory containing support files for makeitso
SUPPORT_FILES_DIRECTORY = os.path.join(KAIJUHOME, "scripts",
                                       "makeitso-gamhelio")

# Path to option descriptions file
OPTION_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
)

# Location of template wsa2gamera.py .ini file
WSA2GAMERA_INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY,
                                       "wsa2gamera_template.ini")

# Path to template .ini file
INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.ini")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")

# Number of seconds in an hour.
SECONDS_PER_HOUR = 3600.0

# Hours to add to tFin for a segment to ensure last restart file is created.
TFIN_NUDGE_HOURS = 0.1


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


def get_run_option(name, description, mode="BASIC", override=None):
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


def read_mjdc_from_fits(wsa_file):
    """Read the MJDc value from a WSA FITS file.

    Read the MJDc value from a WSA FITS file.

    Parameters
    ----------
    wsa_file : str
        Path to WSA FITS file

    Returns
    -------
    mjdc : str
        MJDc value from FITS file, in 'YYYY-MM-DDTHH:MM:SS' format

    Raises
    ------
    None
    """
    # Open the WSA FITS file.
    hdu = fits.open(wsa_file)

    # Read the JD for the center of the map.
    jdc = hdu[0].header['JULDATE']

    # Convert JD to MJD.
    mjdc = jdc - JD2MJD

    # Return the MJDc value.
    return mjdc


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

    # -------------------------------------------------------------------------

    # General options for the simulation
    o = options["simulation"] = {}
    od = option_descriptions["simulation"]

    # Prompt for the name of the job.
    for on in ["job_name"]:
        o[on] = get_run_option(on, od[on], mode)

    # Path to the WSA FITS file to use for initial conditions
    for on in ["wsa_file"]:
        o[on] = get_run_option(on, od[on], mode)

    # Fetch the MJDc value from the WSA FITS file.
    mjdc = read_mjdc_from_fits(o[on])
    o["mjdc"] = mjdc

    date_format = '%Y-%m-%dT%H:%M:%S'
    # Set start and end dates based on fits file MJDc
    od["start_date"]["default"] = MJD2UT(mjdc - Tsolar_synodic/2).strftime(date_format)
    od["stop_date"]["default"] = MJD2UT(mjdc + Tsolar_synodic/2).strftime(date_format)
    # Prompt for the start and stop date of the run.
    for on in ["start_date", "stop_date"]:
        o[on] = get_run_option(on, od[on], mode)

    # Compute the total simulation time in seconds, use as segment duration
    # default in hours.
    start_date = o["start_date"]
    stop_date = o["stop_date"]
    t1 = datetime.datetime.strptime(start_date, date_format)
    t2 = datetime.datetime.strptime(stop_date, date_format)
    simulation_duration = (t2 - t1).total_seconds()/SECONDS_PER_HOUR
    od["segment_duration"]["default"] = str(simulation_duration)

    # Ask if the user wants to split the run into multiple segments.
    # If so, prompt for the segment duration. If not, use the default
    # for the segment duration (which is the simulation duration).
    for on in ["use_segments"]:
        o[on] = get_run_option(on, od[on], mode)
    if o["use_segments"].upper() == "Y":
        for on in ["segment_duration"]:
            o[on] = get_run_option(on, od[on], mode)
    else:
        o["segment_duration"] = od["segment_duration"]["default"]

    # Compute the number of segments based on the simulation duration and
    # segment duration, with 1 additional segment just for spinup. Add 1 if
    # there is a remainder.
    num_segments = simulation_duration/float(o["segment_duration"])
    if num_segments > int(num_segments):
        num_segments += 1
    num_segments = int(num_segments) + 1   # For spinup segment.

    # Prompt for the remaining parameters.
    for on in ["hpc_system"]:
        o[on] = get_run_option(on, od[on], mode)

    # -------------------------------------------------------------------------

    # PBS options
    o = options["pbs"] = {}

    # Common (HPC platform-independent) options
    od = option_descriptions["pbs"]["_common"]
    od["account_name"]["default"] = os.getlogin()
    od["kaiju_install_directory"]["default"] = KAIJUHOME
    od["kaiju_build_directory"]["default"] = os.path.join(KAIJUHOME,
                                                          "build_mpi")
    # Number of segments is computed.
    o["num_segments"] = str(num_segments)
    for on in od:
        if od[on]["LEVEL"] == "COMPUTED":
            continue
        o[on] = get_run_option(on, od[on], mode)

    # HPC platform-specific options
    hpc_platform = options["simulation"]["hpc_system"]
    od = option_descriptions["pbs"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # -------------------------------------------------------------------------

    # wsa2gamera.py options
    options["wsa2gamera"] = {}

    # [Gamera] options
    o = options["wsa2gamera"]["Gamera"] = {}
    od = option_descriptions["wsa2gamera"]["Gamera"]
    od["GridDir"]["default"] = options["pbs"]["run_directory"]
    od["IbcDir"]["default"] = options["pbs"]["run_directory"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [Grid] options
    o = options["wsa2gamera"]["Grid"] = {}
    od = option_descriptions["wsa2gamera"]["Grid"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [WSA] options
    o = options["wsa2gamera"]["WSA"] = {}
    od = option_descriptions["wsa2gamera"]["WSA"]
    # Ensure the simulation setting wsa_file gets copied over 
    # This may be overridden below (duplicate prompt necessary??) 
    # in expert mode
    for on in od:
        o[on] = get_run_option(on, od[on], mode)
        if(on == "wsafile"): o[on] = options["simulation"]["wsa_file"]

    # [Constants] options
    o = options["wsa2gamera"]["Constants"] = {}
    od = option_descriptions["wsa2gamera"]["Constants"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [Normalization] options
    o = options["wsa2gamera"]["Normalization"] = {}
    od = option_descriptions["wsa2gamera"]["Normalization"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # -------------------------------------------------------------------------

    # NOTE: The options for gamhelio are listed under the gamera section of the
    # .ini file.

    # gamhelio options
    options["gamera"] = {}

    # [sim] options
    o = options["gamera"]["sim"] = {}
    od = option_descriptions["gamera"]["sim"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [time] options
    o = options["gamera"]["time"] = {}
    od = option_descriptions["gamera"]["time"]
    od["tFin"]["default"] = str(simulation_duration)
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [spinup] options
    o = options["gamera"]["spinup"] = {}
    od = option_descriptions["gamera"]["spinup"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [output] options
    o = options["gamera"]["output"] = {}
    od = option_descriptions["gamera"]["output"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [physics] options
    o = options["gamera"]["physics"] = {}
    od = option_descriptions["gamera"]["physics"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [prob] options
    o = options["gamera"]["prob"] = {}
    od = option_descriptions["gamera"]["prob"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [iPdir] options
    o = options["gamera"]["iPdir"] = {}
    od = option_descriptions["gamera"]["iPdir"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [jPdir] options
    o = options["gamera"]["jPdir"] = {}
    od = option_descriptions["gamera"]["jPdir"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [kPdir] options
    o = options["gamera"]["kPdir"] = {}
    od = option_descriptions["gamera"]["kPdir"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # <coupling> options
    options["gamera"]["coupling"] = {}
    o = options["gamera"]["coupling"]
    od = option_descriptions["gamera"]["coupling"]
    od["blockHalo"]["default"] = od["blockHalo"]["default"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [restart] options
    o = options["gamera"]["restart"] = {}
    od = option_descriptions["gamera"]["restart"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # -------------------------------------------------------------------------

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
    """
    # Read and create the template for the ini file for wsa2gamera.py, then
    # render and write it.
    with open(WSA2GAMERA_INI_TEMPLATE, encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)
    ini_content = template.render(options)
    ini_file = os.path.join(options["pbs"]["run_directory"], "wsa2gamera.ini")
    with open(ini_file, "w", encoding="utf-8") as f:
        f.write(ini_content)

    # Create the grid and inner boundary conditions files.
    cmd = os.path.join(os.environ["KAIPYHOME"], "kaipy", "scripts", "preproc",
                       "wsa2gamera.py")
    args = [cmd, "wsa2gamera.ini"]
    with open("wsa2gamera.log", "w", encoding="utf-8") as fl:
        status = subprocess.run(args, stdout=fl, stderr=fl)
        if status.returncode == 1:
            print("Error in wsa2gamera! See wsa2gamera.log for details.")
            sys.exit(1)


def create_ini_files(options):
    """Create the gamhelio .ini files from a template.

    Create the gamhelio .ini files from a template.

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

    # Create .ini files.
    if options["simulation"]["use_segments"].upper() == "Y":

        # Create an .ini file for the spinup segment.
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        job = 0
        segment_id = f"{runid}-{job:02d}"
        opt["simulation"]["segment_id"] = segment_id
        # NOTE: Convert hours to seconds
        tFin = float(opt["gamera"]["time"]["tFin"])*SECONDS_PER_HOUR
        dT = float(options["simulation"]["segment_duration"])*SECONDS_PER_HOUR
        tFin_segment = TFIN_NUDGE_HOURS  # Just perform spinup in first segment
        opt["gamera"]["time"]["tFin"] = str(tFin_segment)
        ini_content = template.render(opt)
        ini_file = os.path.join(
            opt["pbs"]["run_directory"],
            f"{opt['simulation']['segment_id']}.ini"
        )
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)
        ini_files.append(ini_file)

        # Create an .ini file for each simulation segment.
        for job in range(1, int(options["pbs"]["num_segments"])):
            opt = copy.deepcopy(options)  # Need a copy of options
            runid = opt["simulation"]["job_name"]
            segment_id = f"{runid}-{job:02d}"
            opt["simulation"]["segment_id"] = segment_id
            opt["gamera"]["restart"]["doRes"] = "T"
            tFin = float(opt["gamera"]["time"]["tFin"])
            dT = float(options["simulation"]["segment_duration"])
            # Nudge the end time so the last restart file is created.
            tFin_segment = job*dT + TFIN_NUDGE_HOURS
            if tFin_segment > tFin:    # Last segment may be shorter.
                tFin_segment = tFin + TFIN_NUDGE_HOURS
            opt["gamera"]["time"]["tFin"] = str(tFin_segment)
            ini_content = template.render(opt)
            ini_file = os.path.join(
                opt["pbs"]["run_directory"],
                f"{opt['simulation']['segment_id']}.ini"
            )
            with open(ini_file, "w", encoding="utf-8") as f:
                f.write(ini_content)
            ini_files.append(ini_file)
    else:
        # Using a single segment for spinup and simulation.
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        job = 0
        segment_id = f"{runid}-{job:02d}"
        opt["simulation"]["segment_id"] = segment_id
        ini_content = template.render(opt)
        ini_file = os.path.join(
            opt["pbs"]["run_directory"],
            f"{opt['simulation']['segment_id']}.ini"
        )
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)
        ini_files.append(ini_file)

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
        cmd = os.path.join(os.environ["KAIPYHOME"], "kaipy", "scripts",
                           "preproc", "XMLGenerator.py")
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
    if options["simulation"]["use_segments"].upper() == "Y":
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
            with open(pbs_script, "w", encoding="utf-8") as f:
                f.write(pbs_content)
            pbs_scripts.append(pbs_script)
    else:
        # Use a single job for spinup and simulation.
        job = 0
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


def main():
    """Main program code for makeitso-gamhelio.

    This is the main program code for makeitso-gamhelio.

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
            raise FileExistsError(f"Options file {path} exists!\nPlease enable the --clobber option if you intended to overwrite.")
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
