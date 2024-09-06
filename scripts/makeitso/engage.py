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
TIEGCMHOME = os.environ["TIEGCMHOME"]
#sys.path.append('/glade/u/home/wiltbemj/src/tiegcm/tiegcmrun')
sys.path.append(f'{TIEGCMHOME}/tiegcmrun')
import tiegcmrun
print(f'tiegcmrum from {tiegcmrun.__file__}')

#import makeitso
KAIJUHOME = os.environ["KAIJUHOME"]
sys.path.append(f'{KAIJUHOME}/scripts/makeitso')
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
OPTION_ENGAGE_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_engage_descriptions.json"
)

OPTION_MAKEITSO_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
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


def create_pbs_scripts(gr_options: dict, makeitso_options:dict,makeitso_pbs_scripts: list, tiegcm_options: dict,tiegcm_inp_scripts:list, tiegcm_pbs_scripts:list):
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
    
    options = copy.deepcopy(gr_options) 
    
    # GRT PBS parameters
    options["pbs"]["mpiexec_command"] = "mpiexec"

    # TIEGCM PBS parameters
    options["pbs"]["tie_nodes"] = tiegcm_options["job"]["resource"]["select"]
    options["pbs"]["tie_ncpus"] = tiegcm_options["job"]["resource"]["ncpus"]
    options["pbs"]["tie_mpiprocs"] = tiegcm_options["job"]["resource"]["mpiprocs"]
    options["pbs"]["tie_mpiranks"] = tiegcm_options["job"]["nprocs"]
    options["pbs"]["tie_exe"] = tiegcm_options["model"]["data"]["coupled_modelexe"]

    # GR PBS parameters
    options["pbs"]["gamera_nodes"] = makeitso_options["pbs"]["select"]
    options["pbs"]["gamera_ncpus"] = makeitso_options["pbs"]["ncpus"]
    options["pbs"]["gamera_mpiprocs"] = makeitso_options["pbs"]["mpiprocs"]
    options["pbs"]["gamera_ompthreads"] = makeitso_options["pbs"]["ompthreads"]
    options["pbs"]["voltron_nodes"] = makeitso_options["pbs"]["select2"]
    options["pbs"]["voltron_ncpus"] = makeitso_options["pbs"]["ncpus"]
    options["pbs"]["voltron_mpiprocs"] = makeitso_options["pbs"]["helper_mpiprocs"]
    options["pbs"]["voltron_ompthreads"] = makeitso_options["pbs"]["helper_ompthreads"]
    options["pbs"]["voltron_mpiranks"] = int(options["pbs"]["gamera_nodes"])*int(options["pbs"]["gamera_mpiprocs"])+int( options["pbs"]["voltron_nodes"])*int(options["pbs"]["voltron_mpiprocs"])
    options["pbs"]["voltron_scripts"] = makeitso_options["pbs"]["mpiexec_command"].replace("mpiexec ", "")

    # Create a PBS script for each segment.
    pbs_scripts = []
    for job in range(1,int(options["pbs"]["num_segments"])):
        opt = copy.deepcopy(options)  # Need a copy of options
        runid = opt["simulation"]["job_name"]
        segment_id = f"{runid}-{job:02d}"
        opt["simulation"]["segment_id"] = segment_id
        
        # TIEGCM input script
        opt["pbs"]["tie_inp"] = tiegcm_inp_scripts[job-1]
        
        pbs_content = template.render(opt)
        pbs_script = os.path.join(
            opt["pbs"]["run_directory"],
            f"{opt['simulation']['segment_id']}.pbs"
        )
        pbs_scripts.append(pbs_script)
        with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)

    # Create a single script which will submit all of the PBS jobs in order.
    submit_all_jobs_script = f"{gr_options['simulation']['job_name']}_pbs.sh"
    with open(submit_all_jobs_script, "w", encoding="utf-8") as f:
        cmd = f"# Standalone GR\n"
        f.write(cmd)
        makeitso_pbs = makeitso_pbs_scripts[0]
        cmd = f"makeitso_job_id=`qsub {makeitso_pbs}`\n"
        f.write(cmd)
        cmd = "echo $makeitso_job_id\n"
        f.write(cmd)
        for makeitso_pbs in makeitso_pbs_scripts[1:]:
            cmd = "old_makeitso_job_id=$makeitso_job_id\n"
            f.write(cmd)
            cmd = f"makeitso_job_id=`qsub -W depend=afterok:$old_makeitso_job_id {makeitso_pbs}`\n"
            f.write(cmd)
            cmd = "echo $makeitso_job_id\n"
            f.write(cmd)
        cmd = f"# Standalone TIEGCM\n"
        f.write(cmd)
        tiegcm_pbs = tiegcm_pbs_scripts[0]
        cmd = f"tiegcm_job_id=`qsub {tiegcm_pbs}`\n"
        f.write(cmd)
        cmd = "echo $tiegcm_job_id\n"
        f.write(cmd)
        for tiegcm_pbs in tiegcm_pbs_scripts[1:]:
            cmd = "old_tiegcm_job_id=$tiegcm_job_id\n"
            f.write(cmd)
            cmd = f"tiegcm_job_id=`qsub -W depend=afterok:$old_tiegcm_job_id {tiegcm_pbs}`\n"
            f.write(cmd)
            cmd = "echo $tiegcm_job_id\n"
            f.write(cmd)
        cmd = f"# Coupled GTR\n"
        f.write(cmd)
        s = pbs_scripts[0]
        cmd = f"job_id=`qsub -W depend=afterok:$makeitso_job_id:$tiegcm_job_id {s}`\n"
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
    os.chmod(submit_all_jobs_script, 0o755)
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
    with open(OPTION_MAKEITSO_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_makeitso_descriptions = json.load(f)
    with open(OPTION_ENGAGE_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_engage_descriptions = json.load(f)

    # Initialize the dictionary of program options.
    options = {}

    #-------------------------------------------------------------------------

    # General options for the simulation
    o = options["simulation"] = {}
    od = option_makeitso_descriptions["simulation"]

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

    # ------------------------------------------------------------------------

    # PBS options
    o = options["pbs"] = {}

    # Common (HPC platform-independent) options
    od = option_makeitso_descriptions["pbs"]["_common"]
    od["account_name"]["default"] = os.getlogin()
    od["kaiju_install_directory"]["default"] = os.environ["KAIJUHOME"]
    od["kaiju_build_directory"]["default"] = os.path.join(
        os.environ["KAIJUHOME"], "build_mpi")
    od["num_segments"]["default"] = str(num_segments)
    for on in od:
        o[on] = makeitso.get_run_option(on, od[on], mode)

    # HPC platform-specific options
    hpc_platform = options["simulation"]["hpc_system"]
    gamera_grid_type = options["simulation"]["gamera_grid_type"]
    od = option_makeitso_descriptions["pbs"][hpc_platform]
    od["select"]["default"] = od["select"]["default"][gamera_grid_type]
    od["num_helpers"]["default"] = (
        od["num_helpers"]["default"][gamera_grid_type]
    )
    if hpc_platform == "derecho":
        od["modules"]["default"] = ["ncarenv/23.09","cmake/3.26.3","craype/2.7.31","intel-classic/2023.2.1","cray-mpich/8.1.27","ncarcompilers/1.0.0","mkl/2023.2.0","hdf5-mpi/1.12.2","netcdf-mpi/4.9.2","esmf/8.6.0","conda"]
    for on in od:
        o[on] = makeitso.get_run_option(on, od[on], mode)

    #-------------------------------------------------------------------------
    # coupling options
    options["coupling"] = {}
    o = options["coupling"]
    od = option_engage_descriptions["coupling"]

    od["conda_env"]["default"] = os.environ.get('CONDA_DEFAULT_ENV')

    # Prompt for the remaining parameters.
    for on in ["gamera_spin_up_time", "gcm_spin_up_time", 
               "root_directory","conda_env","tfin_delta","doGCM"]:
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
    path = f"{options['simulation']['job_name']}-engage.json"
    if os.path.exists(path):
        if not clobber:
            raise FileExistsError(f"Options file {path} exists!")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)

    makeitso_args = {'clobber': True, 'debug': True, 'verbose': True}
    makeitso_args.update(options)
    #print(f"makeitso_args = {makeitso_args}")

    makeitso_options, makeitso_spinup_pbs_scripts, makeitso_warmup_pbs_scripts = makeitso.makeitso(makeitso_args)
    makeitso_pbs_scripts = makeitso_spinup_pbs_scripts + makeitso_warmup_pbs_scripts
    # Save the makeitso options dictionary as a JSON file in the current directory.
 
    with open('makeitso_parameters.json', 'w') as f:
        json.dump(makeitso_options, f, indent=JSON_INDENT)
        json.dump(makeitso_spinup_pbs_scripts, f, indent=JSON_INDENT)
        json.dump(makeitso_warmup_pbs_scripts, f, indent=JSON_INDENT)
    

    # Run the TIEGCMrun
    coupled_options = copy.deepcopy(options)
    coupled_options["voltron"] = makeitso_options.get("voltron", {})
    print(f"coupled_options={coupled_options}")
    
    tiegcm_args = [
    "--coupling",
    "--engage", json.dumps(coupled_options)
    ]   
    
    tiegcm_options,tiegcm_pbs_scripts,tiegcm_inp_scripts = tiegcmrun.tiegcmrun(tiegcm_args)

   
    # Save the tiegcm options dictionary as a JSON file in the current directory.
    with open('tiegcmrun_parameters.json', 'w') as f:
        json.dump(tiegcm_options, f, indent=JSON_INDENT)
        json.dump(tiegcm_pbs_scripts, f, indent=JSON_INDENT)
        json.dump(tiegcm_inp_scripts, f, indent=JSON_INDENT)

    # Create the PBS job scripts.
    pbs_scripts, submit_all_jobs_script = create_pbs_scripts(options,makeitso_options, makeitso_pbs_scripts, tiegcm_options, tiegcm_inp_scripts, tiegcm_pbs_scripts)
    print(f"pbs_scripts = {pbs_scripts}")
    print(f"submit_all_jobs_script = {submit_all_jobs_script}")
    

if __name__ == "__main__":
    """Begin main program."""
    main()