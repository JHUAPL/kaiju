#!/usr/bin/env python

"""Run the MAGE comparative test cases

This script runs the MAGE comparative test cases

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import shutil
import subprocess
import sys

# Import 3rd-party modules.
from jinja2 import Template

# Import project modules.
import common

# Program constants

# Program description.
DESCRIPTION = 'Run the MAGE comparative test cases.'

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for results
TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'compTest')

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Home directory of kaipy installation
KAIPYHOME = os.environ['KAIPYHOME']

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for tests
COMP_TESTS_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'comp_tests.lst')

# Prefix for comp tests directory name
COMP_TESTS_DIRECTORY_PREFIX = 'compTest_'

# Subdirectory of build directory containing compiled products to use in tests
BIN_DIR = 'bin'

# Path to jinja2 template file for PBS script.
PBS_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'relativeCase-template.pbs'
)

# Name of rendered PBS script.
PBS_SCRIPT = 'relativeCase.pbs'


# Path to jinja2 template file for PBS script.
POSTPROC_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'relCasePost-template.pbs'
)

# Path to jinja2 template file for combined Slack report PBS script.
SLACK_REPORT_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'relCaseReport-template.pbs'
)


# Path to jinja2 template file for XML file
XML_TEMPLATE = os.path.join(
    TEST_SCRIPTS_DIRECTORY, 'relativeCaseGo-template.xml'
)

# Name of rendered XML file
XML_FILE = 'relativeCaseGo.xml'

# debug and verbose globals
debug = False
verbose = False

def processComparativeResults(case1_jobid, case2_jobid, caseName, pbsTemplate, pbs_options):
    # Render the job template.
    pbs_content = pbsTemplate.render(pbs_options)
    pbsFilename = f"{caseName}.pbs"
    with open(pbsFilename, 'w', encoding='utf-8') as f:
        f.write(pbs_content)
    # Submit the job
    if verbose:
        print('Submitting comparative tests model post processing.')
    cmd = f"qsub -W depend=afterok:{case1_jobid}:{case2_jobid} {pbsFilename}"
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print('ERROR: qsub failed.\n'
              f"e.cmd = {e.cmd}\n"
              f"e.returncode = {e.returncode}\n"
              'See test log for output.\n',
              file=sys.stderr)
    job_id = cproc.stdout.rstrip()
    if debug:
        print(f"job_id = {job_id}")
    return job_id


def create_combined_report(postproc_job_id_1, postproc_job_id_2, caseName, pbsTemplate, pbs_options):
    # Render the job template.
    pbs_content = pbsTemplate.render(pbs_options)
    pbsFilename = f"{caseName}.pbs"
    with open(pbsFilename, 'w', encoding='utf-8') as f:
        f.write(pbs_content)
    # Submit the job
    if verbose:
        print('Submitting combined Slack report job.')
    cmd = f"qsub -W depend=afterok:{postproc_job_id_1}:{postproc_job_id_2} {pbsFilename}"
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print('ERROR: qsub failed.\n'
              f"e.cmd = {e.cmd}\n"
              f"e.returncode = {e.returncode}\n"
              'See test log for output.\n',
              file=sys.stderr)
    job_id = cproc.stdout.split('.')[0]
    if debug:
        print(f"job_id = {job_id}")

    # Return the job ID.
    return job_id


def generateAndRunCase(caseName,pbsTemplate,pbs_options,xmlTemplate,xml_options,wait_job_id=None):
    os.mkdir(caseName)
    os.chdir(caseName)
    # Render the job template.
    pbs_content = pbsTemplate.render(pbs_options)
    with open(PBS_SCRIPT, 'w', encoding='utf-8') as f:
        f.write(pbs_content)
    # Render the xml file.
    xml_content = xmlTemplate.render(xml_options)
    with open(XML_FILE, 'w', encoding='utf-8') as f:
        f.write(xml_content)
    shutil.copy2('../voltron.x', './voltron.x')
    shutil.copy2('../voltron_mpi.x', './voltron_mpi.x')
    shutil.copy2('../lfmD.h5', './lfmD.h5')
    shutil.copy2('../rcmconfig.h5', './rcmconfig.h5')
    shutil.copy2('../bcwind.h5', './bcwind.h5')
    # Submit the job
    if verbose:
        print(f"Submitting comparative tests model run for case {caseName}.")
    cmd = f"qsub {PBS_SCRIPT}"
    if wait_job_id is not None:
        cmd = f"qsub -W depend=afterok:{wait_job_id} {PBS_SCRIPT}"
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print('ERROR: qsub failed.\n'
              f"e.cmd = {e.cmd}\n"
              f"e.returncode = {e.returncode}\n"
              'See test log for output.\n',
              file=sys.stderr)
    job_id = cproc.stdout.split('.')[0]
    if debug:
        print(f"job_id = {job_id}")
    os.chdir('..')
    return job_id

def generateMpi32RestartRelease(pbsTemplate,xmlTemplate,base_pbs_options,wait_job_id):
    caseName = "relMpi32ResRelease"

    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '02:00:00'
    pbs_options['is_mpi'] = 'True' #Jinja considers this 'Truthy'
    pbs_options['gamera_ranks'] = '3'
    pbs_options['gamera_per_node'] = '2'
    pbs_options['copy_case_folder'] = 'relMpi32Release' # Jinja considers this 'Truthy'
    pbs_options['copy_case_runid'] = 'msphere_M32_R' #Jinja considers this 'Truthy'

    xml_options = {}
    xml_options['serial_coupling'] = 'F'
    xml_options['runid'] = 'msphere_M32_R'
    xml_options['do_restart'] = 'T'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '1'
    xml_options['num_i_ranks'] = '3'
    xml_options['num_j_ranks'] = '2'
    xml_options['num_g_th'] = '64'

    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options, wait_job_id)

def generateMpi32Release(pbsTemplate,xmlTemplate,base_pbs_options):
    caseName = "relMpi32Release"

    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '02:00:00'
    pbs_options['is_mpi'] = 'True' #Jinja considers this 'Truthy'
    pbs_options['gamera_ranks'] = '3'
    pbs_options['gamera_per_node'] = '2'
    pbs_options['copy_case_folder'] = '' # Jinja considers this 'Falsy'
    pbs_options['copy_case_runid'] = '' #Jinja considers this 'Falsy'

    xml_options = {}
    xml_options['serial_coupling'] = 'F'
    xml_options['runid'] = 'msphere_M32_R'
    xml_options['do_restart'] = 'F'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '-1'
    xml_options['num_i_ranks'] = '3'
    xml_options['num_j_ranks'] = '2'
    xml_options['num_g_th'] = '64'

    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options)

def generateMpi31Release(pbsTemplate,xmlTemplate,base_pbs_options):
    caseName = "relMpi31Release"

    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '02:00:00'
    pbs_options['is_mpi'] = 'True' #Jinja considers this 'Truthy'
    pbs_options['gamera_ranks'] = '3'
    pbs_options['gamera_per_node'] = '1'
    pbs_options['copy_case_folder'] = '' # Jinja considers this 'Falsy'
    pbs_options['copy_case_runid'] = '' #Jinja considers this 'Falsy'

    xml_options = {}
    xml_options['serial_coupling'] = 'F'
    xml_options['runid'] = 'msphere_M31_R'
    xml_options['do_restart'] = 'F'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '-1'
    xml_options['num_i_ranks'] = '3'
    xml_options['num_j_ranks'] = '1'
    xml_options['num_g_th'] = '128'

    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options)

def generateMpi12Release(pbsTemplate,xmlTemplate,base_pbs_options):
    caseName = "relMpi12Release"

    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '02:00:00'
    pbs_options['is_mpi'] = 'True' #Jinja considers this 'Truthy'
    pbs_options['gamera_ranks'] = '1'
    pbs_options['gamera_per_node'] = '2'
    pbs_options['copy_case_folder'] = '' # Jinja considers this 'Falsy'
    pbs_options['copy_case_runid'] = '' #Jinja considers this 'Falsy'

    xml_options = {}
    xml_options['serial_coupling'] = 'F'
    xml_options['runid'] = 'msphere_M12_R'
    xml_options['do_restart'] = 'F'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '-1'
    xml_options['num_i_ranks'] = '1'
    xml_options['num_j_ranks'] = '2'
    xml_options['num_g_th'] = '64'

    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options)

def generateMpi11SynchRelease(pbsTemplate,xmlTemplate,base_pbs_options):
    caseName = "relMpi11SynchRelease"

    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '04:00:00'
    pbs_options['is_mpi'] = 'True' #Jinja considers this 'Truthy'
    pbs_options['gamera_ranks'] = '1'
    pbs_options['gamera_per_node'] = '1'
    pbs_options['copy_case_folder'] = '' # Jinja considers this 'Falsy'
    pbs_options['copy_case_runid'] = '' #Jinja considers this 'Falsy'

    xml_options = {}
    xml_options['serial_coupling'] = 'T'
    xml_options['runid'] = 'msphere_M11S_R'
    xml_options['do_restart'] = 'F'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '-1'
    xml_options['num_i_ranks'] = '1'
    xml_options['num_j_ranks'] = '1'
    xml_options['num_g_th'] = '128'

    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options)

def generateMpi11Release(pbsTemplate,xmlTemplate,base_pbs_options):
    caseName = "relMpi11Release"

    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '04:00:00'
    pbs_options['is_mpi'] = 'True' #Jinja considers this 'Truthy'
    pbs_options['gamera_ranks'] = '1'
    pbs_options['gamera_per_node'] = '1'
    pbs_options['copy_case_folder'] = '' # Jinja considers this 'Falsy'
    pbs_options['copy_case_runid'] = '' #Jinja considers this 'Falsy'

    xml_options = {}
    xml_options['serial_coupling'] = 'F'
    xml_options['runid'] = 'msphere_M11_R'
    xml_options['do_restart'] = 'F'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '-1'
    xml_options['num_i_ranks'] = '1'
    xml_options['num_j_ranks'] = '1'
    xml_options['num_g_th'] = '128'

    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options)

def generateSerialRestartRelease(pbsTemplate,xmlTemplate,base_pbs_options,wait_job_id):
    caseName = "relSerialResRelease"

    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '04:00:00'
    pbs_options['is_mpi'] = '' #Jinja considers this 'Falsy'
    pbs_options['gamera_ranks'] = '1'
    pbs_options['gamera_per_node'] = '1'
    pbs_options['copy_case_folder'] = 'relSerialRelease' # Jinja considers this 'Truthy'
    pbs_options['copy_case_runid'] = 'msphere_S_R' #Jinja considers this 'Truthy'

    xml_options = {}
    xml_options['serial_coupling'] = 'F'
    xml_options['runid'] = 'msphere_S_R'
    xml_options['do_restart'] = 'T'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '1'
    xml_options['num_i_ranks'] = '1'
    xml_options['num_j_ranks'] = '1'
    xml_options['num_g_th'] = '128'

    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options, wait_job_id)

def generateSerialRelease(pbsTemplate,xmlTemplate,base_pbs_options):
    caseName = "relSerialRelease"
    
    pbs_options = base_pbs_options
    pbs_options['job_name'] = caseName
    pbs_options['walltime'] = '04:00:00'
    pbs_options['is_mpi'] = '' #Jinja considers this 'Falsy'
    pbs_options['gamera_ranks'] = '1'
    pbs_options['gamera_per_node'] = '1'
    pbs_options['copy_case_folder'] = '' # Jinja considers this 'Falsy'
    pbs_options['copy_case_runid'] = '' #Jinja considers this 'Falsy'
    
    xml_options = {}
    xml_options['serial_coupling'] = 'F'
    xml_options['runid'] = 'msphere_S_R'
    xml_options['do_restart'] = 'F'
    xml_options['restart_runid'] = xml_options['runid']
    xml_options['restart_number'] = '-1'
    xml_options['num_i_ranks'] = '1'
    xml_options['num_j_ranks'] = '1'
    xml_options['num_g_th'] = '128'
    
    return generateAndRunCase(caseName, pbsTemplate, pbs_options, xmlTemplate, xml_options)

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
    # set up global verbose and debug variables
    global debug,verbose

    # Set up the command-line parser.
    parser = common.create_command_line_parser(DESCRIPTION)

    # Add additional arguments
    parser.add_argument(
        '--allTests', '-a', action='store_true',
        help='Perform all comparative tests (default: %(default)s).'
    )

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    be_loud = args.loud
    slack_on_fail = args.slack_on_fail
    is_test = args.test
    verbose = args.verbose
    allTests = args.allTests

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Set up for communication with Slack.
    if verbose:
        print('Creating Slack client.')
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the tests.
    if verbose:
        print(f"Creating {TEST_DIRECTORY}.")
    os.mkdir(TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    with open(COMP_TESTS_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Read the template for the PBS script used for the test data generation.
    with open(PBS_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    pbs_template = Template(template_content)
    if debug:
        print(f"comp_test_pbs_template = {pbs_template}")

    # ------------------------------------------------------------------------

    # Read the template for the PBS script used for the test data post processing
    with open(POSTPROC_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    postproc_template = Template(template_content)
    if debug:
        print(f"comp_test_postproc_template = {postproc_template}")

    # ------------------------------------------------------------------------

    # Read the template for the PBS script used for posting reports to Slack.
    with open(SLACK_REPORT_TEMPLATE, 'r', encoding='utf-8') as f:
        template_content = f.read()
    slack_report_template = Template(template_content)
    if debug:
        print(f"slack_report_template = {slack_report_template}")

    # ------------------------------------------------------------------------


    # Read the template for the XML file used for the test data generation.
    with open(XML_TEMPLATE, 'r', encoding='utf-8') as f:
        xml_content = f.read()
    xml_template = Template(xml_content)
    if debug:
        print(f"comp_test_xml_file = {xml_template}")

    # ------------------------------------------------------------------------

    # Create the make command to build the code.
    make_cmd = 'make voltron_mpi.x; make voltron.x'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # ------------------------------------------------------------------------

    # Create the list for submit results. Only set to True if all build and
    # qsub commands for a set are OK.
    submit_ok = []

    # Create the list of job IDs.
    job_ids = []

    # Run the comparative tests with each set of modules.
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing comparative tests with module set '
                  f"{module_list_file}")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_set_name = {module_set_name}.")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
        if debug:
            print(f"path = {path}")
        if verbose:
            print(f"Reading module list file {path}.")
        module_names, cmake_environment, cmake_options = (
            common.read_build_module_list_file(path)
        )
        if debug:
            print(f"module_names = {module_names}")
            print(f"cmake_environment = {cmake_environment}")
            print(f"cmake_options = {cmake_options}")

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Add the cmake option for the weekly dash build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=Release'
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Make a directory for this build, and go there.
        dir_name = f"{COMP_TESTS_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(TEST_DIRECTORY, dir_name)
        if verbose:
            print(f"Creating and moving to build directory {build_directory}.")
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Run cmake to build the Makefile.
        if verbose:
            print(
                'Running cmake to create Makefile for module set'
                f" {module_set_name}."
            )
        cmd = (
            f"{module_cmd}; {cmake_environment} cmake {cmake_options}"
            f" {KAIJUHOME} >& cmake.out"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            # NOTE: stdout and stderr goes to stdout (into log file)
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: cmake for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'cmake.out')}"
                ' for output from cmake.\n'
                f"Skipping remaining steps for module set {module_set_name}",
                file=sys.stderr
            )
            continue

        # Run the build.
        if verbose:
            print(
                'Running make to build kaiju for module set'
                f" {module_set_name}."
            )
        cmd = f"{module_cmd}; {make_cmd} >& make.out"
        if debug:
            print(f"cmd = {cmd}")
        try:
            # NOTE: stdout and stderr go into make.out.
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: make for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'make.out')}"
                ' for output from make.\n'
                f"Skipping remaining steps for module set {module_set_name}",
                file=sys.stderr
            )
            continue

        # Move into the bin directory to run the tests.
        os.chdir(BIN_DIR)
        os.mkdir('vidData')

        # Generate the LFM grid file.
        if verbose:
            print('Creating LFM grid file.')
        cmd = 'genLFM.py -gid D'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: Unable to create LFM grid file for module set '
                  f"{module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See testing log for output from genLFM.py.\n'
                  'Skipping remaining steps for module set'
                  f"{module_set_name}\n")
            continue

        # Generate the solar wind boundary condition file.
        if verbose:
            print('Creating solar wind initial conditions file.')
        cmd = 'cda2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: Unable to create solar wind boundary conditions file'
                  f" for module set {module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See testing log for output from cda2wind.py.\n'
                  'Skipping remaining steps for module set'
                  f"{module_set_name}\n")
            continue

        # Generate the RCM configuration file.
        if verbose:
            print('Creating RCM configuration file.')
        cmd = 'genRCM.py'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('ERROR: Unable to create RCM configuration file'
                  f" for module set {module_set_name}.\n"
                  f"e.cmd = {e.cmd}\n"
                  f"e.returncode = {e.returncode}\n"
                  'See testing log for output from genRCM.py.\n'
                  'Skipping remaining steps for module set '
                  f"{module_set_name}\n")
            continue

        # Assemble data to fill in the PBS template.
        base_pbs_options = {}
        base_pbs_options['account'] = os.environ['DERECHO_TESTING_ACCOUNT']
        base_pbs_options['queue'] = os.environ['DERECHO_TESTING_QUEUE']
        base_pbs_options['job_priority'] = os.environ['DERECHO_TESTING_PRIORITY']
        base_pbs_options['modules'] = module_names
        base_pbs_options['kaijuhome'] = KAIJUHOME
        base_pbs_options['kaipyhome'] = KAIPYHOME
        base_pbs_options['tmpdir'] = os.environ['TMPDIR']
        base_pbs_options['slack_bot_token'] = os.environ['SLACK_BOT_TOKEN']
        base_pbs_options['mage_test_root'] = os.environ['MAGE_TEST_ROOT']
        base_pbs_options['mage_test_set_root'] = os.environ['MAGE_TEST_SET_ROOT']
        base_pbs_options['branch_or_commit'] = os.environ['BRANCH_OR_COMMIT']
        base_pbs_options['report_options'] = ''
        if debug:
            base_pbs_options['report_options'] += ' -d'
        base_pbs_options['report_options'] += ' -l'  # Always report.
        # if slack_on_fail:
        #     base_pbs_options['report_options'] += ' -s'
        if is_test:
            base_pbs_options['report_options'] += ' -t'
        if verbose:
            base_pbs_options['report_options'] += ' -v'
        # if allTests:
        #     base_pbs_options['report_options'] += ' -a'

        # Create space to store job ids and submit OKs
        job_ids.append([])
        submit_ok.append([])

        if not allTests:
            if debug:
                print("Performing basic tests")
            # Serial Run
            job_id_s_r = generateSerialRelease(pbs_template,xml_template,base_pbs_options)
            
            # Record the job ID.
            job_ids[i_module_set].append(job_id_s_r)
        
            # Record successful submission.
            submit_ok[i_module_set].append(True)
            
            # MPI 3x2 Run
            job_id_m32_r = generateMpi32Release(pbs_template,xml_template,base_pbs_options)
            job_ids[i_module_set].append(job_id_m32_r)
            submit_ok[i_module_set].append(True)
            
            # Restarted MPI 3x2 Run
            job_id_m32r_r = generateMpi32RestartRelease(pbs_template,xml_template,base_pbs_options,job_id_m32_r)
            job_ids[i_module_set].append(job_id_m32r_r)
            submit_ok[i_module_set].append(True)
            
            # Submit post-processing for these runs
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'SerialMpiComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relSerialRelease'
            postProcOpts['case1id'] = 'msphere_S_R'
            postProcOpts['case2F'] = 'relMpi32Release'
            postProcOpts['case2id'] = 'msphere_M32_R'
            postProcOpts['ts'] = '0'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            postproc_job_id_1 = processComparativeResults(job_id_s_r, job_id_m32_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            job_ids[i_module_set].append(postproc_job_id_1)
            submit_ok[i_module_set].append(True)

            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'MpiRestartComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relMpi32Release'
            postProcOpts['case1id'] = 'msphere_M32_R'
            postProcOpts['case2F'] = 'relMpi32ResRelease'
            postProcOpts['case2id'] = 'msphere_M32_R'
            postProcOpts['ts'] = '50'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            postproc_job_id_2 = processComparativeResults(job_id_m32_r, job_id_m32r_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            job_ids[i_module_set].append(postproc_job_id_2)
            submit_ok[i_module_set].append(True)

            # Generate the combined report.
            postProcOpts = base_pbs_options
            report_job_id = create_combined_report(postproc_job_id_1, postproc_job_id_2, "combined_report", slack_report_template, base_pbs_options)
            job_ids[i_module_set].append(report_job_id)
            submit_ok[i_module_set].append(True)

        else:
            if debug:
                print("Performing full test suite")
             # Serial Run
            job_id_s_r = generateSerialRelease(pbs_template,xml_template,base_pbs_options)
            
            # Record the job ID.
            job_ids[i_module_set].append(job_id_s_r)
            
            # Record successful submission.
            submit_ok[i_module_set].append(True)
            
            # Restarted Serial Run
            job_id_sr_r = generateSerialRestartRelease(pbs_template,xml_template,base_pbs_options,job_id_s_r)
            job_ids[i_module_set].append(job_id_sr_r)
            submit_ok[i_module_set].append(True)
            
            # Synchronous MPI 1x1 Run
            job_id_m11s_r = generateMpi11SynchRelease(pbs_template,xml_template,base_pbs_options)
            job_ids[i_module_set].append(job_id_m11s_r)
            submit_ok[i_module_set].append(True)
            
            # MPI 1x1 Run
            job_id_m11_r = generateMpi11Release(pbs_template,xml_template,base_pbs_options)
            job_ids[i_module_set].append(job_id_m11_r)
            submit_ok[i_module_set].append(True)
            
            # MPI 1x2 Run
            job_id_m12_r = generateMpi12Release(pbs_template,xml_template,base_pbs_options)
            job_ids[i_module_set].append(job_id_m12_r)
            submit_ok[i_module_set].append(True)
            
            # MPI 3x1 Run
            job_id_m31_r = generateMpi31Release(pbs_template,xml_template,base_pbs_options)
            job_ids[i_module_set].append(job_id_m31_r)
            submit_ok[i_module_set].append(True)
            
            # MPI 3x2 Run
            job_id_m32_r = generateMpi32Release(pbs_template,xml_template,base_pbs_options)
            job_ids[i_module_set].append(job_id_m32_r)
            submit_ok[i_module_set].append(True)
            
            # Restarted MPI 3x2 Run
            job_id_m32r_r = generateMpi32RestartRelease(pbs_template,xml_template,base_pbs_options,job_id_m32_r)
            job_ids[i_module_set].append(job_id_m32r_r)
            submit_ok[i_module_set].append(True)
            
            # Submit post-processing for these runs
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'SerialRestartComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relSerialRelease'
            postProcOpts['case1id'] = 'msphere_S_R'
            postProcOpts['case2F'] = 'relSerialResRelease'
            postProcOpts['case2id'] = 'msphere_S_R'
            postProcOpts['ts'] = '50'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            processComparativeResults(job_id_s_r, job_id_sr_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'SerialToMpiComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relSerialRelease'
            postProcOpts['case1id'] = 'msphere_S_R'
            postProcOpts['case2F'] = 'relMpi11SynchRelease'
            postProcOpts['case2id'] = 'msphere_M11S_R'
            postProcOpts['ts'] = '0'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            processComparativeResults(job_id_s_r, job_id_m11s_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'AsynchComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relMpi11SynchRelease'
            postProcOpts['case1id'] = 'msphere_M11S_R'
            postProcOpts['case2F'] = 'relMpi11Release'
            postProcOpts['case2id'] = 'msphere_M11_R'
            postProcOpts['ts'] = '0'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            processComparativeResults(job_id_m11s_r, job_id_m11_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'JDecompComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relMpi11Release'
            postProcOpts['case1id'] = 'msphere_M11_R'
            postProcOpts['case2F'] = 'relMpi12Release'
            postProcOpts['case2id'] = 'msphere_M12_R'
            postProcOpts['ts'] = '0'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            processComparativeResults(job_id_m11_r, job_id_m12_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'IDecompComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relMpi11Release'
            postProcOpts['case1id'] = 'msphere_M11_R'
            postProcOpts['case2F'] = 'relMpi31Release'
            postProcOpts['case2id'] = 'msphere_M31_R'
            postProcOpts['ts'] = '0'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            processComparativeResults(job_id_m11_r, job_id_m31_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'MPIDecompComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relMpi11Release'
            postProcOpts['case1id'] = 'msphere_M11_R'
            postProcOpts['case2F'] = 'relMpi32Release'
            postProcOpts['case2id'] = 'msphere_M32_R'
            postProcOpts['ts'] = '0'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            processComparativeResults(job_id_m11_r, job_id_m32_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            
            postProcOpts = base_pbs_options
            postProcOpts['caseName'] = 'MpiRestartComp'
            postProcOpts['frameFolder'] = 'vidData'
            postProcOpts['case1F'] = 'relMpi32Release'
            postProcOpts['case1id'] = 'msphere_M32_R'
            postProcOpts['case2F'] = 'relMpi32ResRelease'
            postProcOpts['case2id'] = 'msphere_M32_R'
            postProcOpts['ts'] = '50'
            postProcOpts['te'] = '120'
            postProcOpts['dt'] = '60'
            processComparativeResults(job_id_m32_r, job_id_m32r_r, postProcOpts['caseName'], postproc_template, postProcOpts)
            
        # Save the job number in a file.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            f.write(f"{job_ids[i_module_set]}\n")

        # End of loop over module sets.

    if debug:
        print(f"submit_ok = {submit_ok}")
        print(f"job_ids = {job_ids}")

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += (
        f"Test results are in {TEST_DIRECTORY}.\n"
    )
    for (i_module_set, module_list_file) in enumerate(module_list_files):
        if not all(submit_ok[i_module_set]):
            test_report_details_string += (
                f"Module set `{module_list_file}` submission: *FAILED*\n"
            )
            continue
        test_report_details_string += (
            f"Comparative tests for module set `{module_list_file}` submitted as "
            f"job {job_ids[i_module_set]}.\n"
        )

    # Summarize the test results.
    if 'FAILED' in test_report_details_string:
        test_report_summary_string = 'Comparative Tests submission: *FAILED*\n'
    else:
        test_report_summary_string = 'Comparative Tests submission: *PASSED*\n'

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If a test failed, or loud mode is on, post report to Slack.
    if (slack_on_fail and 'FAILED' in test_report_details_string) or be_loud:
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
