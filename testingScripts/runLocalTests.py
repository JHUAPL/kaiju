#!/usr/bin/env python

import argparse
import os
import subprocess
import pathlib
import datetime

def create_command_line_parser():
    """Create the command-line argument parser.
    
    Create the parser for command-line arguments.
    
    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.
    
    """
    parser = argparse.ArgumentParser(description="Script to help setup automated tests within a kaiju repo")
    
    parser.add_argument(
        "-A", default="",
        help="Charge code to use when running tests."
    )
    parser.add_argument(
        "-ce", default="",
        help="Conda environment name to load with conda module"
    )
    parser.add_argument(
        "-p", default="economy",
        help="Job priority to use when running tests (default: %(default)s)."
    )
    
    parser.add_argument(
        "--unitTests",  action='store_true',default=False,
        help="Run unit tests (default: %(default)s)."
    )
    parser.add_argument(
        "--weeklyDash",  action='store_true',default=False,
        help="Run weekly dash (default: %(default)s)."
    )
    parser.add_argument(
        "--compTests",  action='store_true',default=False,
        help="Run default subset of comparative tests (default: %(default)s)."
    )
    parser.add_argument(
        "--compTestsFull",  action='store_true',default=False,
        help="Run full suite of comparative tests (over-rides --compTests) (default: %(default)s)."
    )
    parser.add_argument(
        "--buildTests",  action='store_true',default=False,
        help="Run build tests (default: %(default)s)."
    )
    parser.add_argument(
        "--icTests",  action='store_true',default=False,
        help="Run tests to build Initial Condition files (default: %(default)s)."
    )
    parser.add_argument(
        "--intelChecks",  action='store_true',default=False,
        help="Run Intel Inspector memory and thread check tests (default: %(default)s)."
    )
    parser.add_argument(
        "--reproTests",  action='store_true',default=False,
        help="Run reproducibility tests (default: %(default)s)."
    )
    
    parser.add_argument(
        "--all",  action='store_true',default=False,
        help="Run all tests (default: %(default)s)."
    )
    
    return parser

def main():
    """Helper script to run automated tests locally inside a kaiju repository
    """
    # Set up the command-line parser.
    parser = create_command_line_parser()
    
    # Parse the command-line arguments.
    args = parser.parse_args()

    # Adjust test options
    if args.all:
        args.unitTests = True
        args.weeklyDash = True
        args.compTests = True
        args.compTestsFull = True
        args.buildTests = True
        args.icTests = True
        args.intelChecks = True
        args.reproTests = True

    if args.compTestsFull:
        args.compTests = False
    
    if not (args.unitTests or args.weeklyDash or args.compTests or args.compTestsFull or
            args.buildTests or args.icTests or args.intelChecks or args.reproTests):
        parser.print_help()
        exit()

    # find repo home directory
    called_from = os.path.dirname(os.path.abspath(__file__))
    os.chdir(called_from)
    os.chdir('..')
    homeDir = os.getcwd()

    # Check for necessary environment variables
    if len(args.ce) == 0 and 'CONDA_DEFAULT_ENV' not in os.environ:
        print("A conda environment name was not supplied, and a currently loaded conda environment could not be determined.")
        print("Please either supply the name of a conda environment with the '-ce <name>' option,")
        print("    or load an entironment before running this script, and it should be automatically found.")
        exit()
    elif len(args.ce) == 0:
        args.ce = os.environ['CONDA_DEFAULT_ENV']
        print(f"Automatically setting conda environment to {args.ce}")
    if len(args.A) == 0 and (args.unitTests or args.weeklyDash or
            args.compTests or args.compTestsFull or args.intelChecks or args.reproTests):
        print("A charge code with not supplied, but the requested tests require one.")
        print("Please supply a charge code with the -A # option.")
        exit()
    if 'KAIJUHOME' not in os.environ:
        os.environ['KAIJUHOME'] = homeDir
        print(f"Running tests out of local git repository: {homeDir}")
    if pathlib.Path(homeDir).resolve() != pathlib.Path(os.environ['KAIJUHOME']).resolve():
        print("The setupEnvironment.sh script must be sourced for the repo this script resides in before calling it.")
        exit()
    if 'KAIPYHOME' not in os.environ and (args.weeklyDash or args.compTests or args.compTestsFull):
        print("The 'KAIPYHOME' environment variable was not set, but the requested tests require it.")
        print("The setupEnvironment.sh script for ANY kaipy repo must be sourced before running these tests.")
        exit()
    elif 'KAIPYHOME' not in os.environ:
        os.environ['KAIPYHOME'] = ""
    
    # Set environment variables
    os.environ['MAGE_TEST_ROOT'] = homeDir
    os.environ['MAGE_TEST_RUNS_ROOT']=os.path.join(os.environ['MAGE_TEST_ROOT'],"test_runs")
    os.environ['DERECHO_TESTING_ACCOUNT'] = args.A
    os.environ['DERECHO_TESTING_QUEUE'] = 'main'
    os.environ['DERECHO_TESTING_PRIORITY'] = args.p
    os.environ['SLACK_BOT_TOKEN'] = '' # help ensure no accidental reporting to slack
    os.environ['PYTHONUNBUFFERED']='TRUE'
    os.environ['CONDA_ENVIRONMENT']=args.ce
    gitBranch  = subprocess.run(['git','branch','--show-current'], stdout=subprocess.PIPE).stdout.decode('utf-8')
    if 'not a git repository' in gitBranch:
        print("This script must be executed inside a kaiju git repository")
        exit()
    gitBranch = gitBranch.strip()
    os.environ['BRANCH_OR_COMMIT'] = gitBranch
    currenttime = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    test_set_dir = f"{currenttime}-{gitBranch}"
    os.environ['MAGE_TEST_SET_ROOT'] = os.path.join(os.environ['MAGE_TEST_RUNS_ROOT'],test_set_dir)
    os.makedirs(os.environ['MAGE_TEST_SET_ROOT'], exist_ok=True)
    os.chdir(os.environ['MAGE_TEST_SET_ROOT'])
    os.environ['CONDARC'] = '' # these must be specified to avoid errors
    os.environ['CONDA_ENVS_PATH'] = ''
    os.environ['KAIPY_PRIVATE_ROOT'] = os.environ['KAIPYHOME'] # some scripts use this alternate

    
    print(f"Running tests on branch {gitBranch}")
    if len(args.A) > 0:
        print(f"Using charge code {args.A} with priority {args.p}")
    print(f"Running in folder test_runs/{test_set_dir}")
    print("")
    
    # Run Tests
    if args.unitTests:
        print("Running unit tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','unitTest.py'),'-tv'])
    if args.weeklyDash:
        print("Running weekly dash")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','weeklyDash.py'),'-tv'])
    if args.compTests:
        print("Running default comparative tests subset")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','relativeTests.py'),'-tv'])
    if args.compTestsFull:
        print("Running full comparative tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','relativeTests.py'),'-tva'])
    if args.buildTests:
        print("Running build tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','buildTest.py'),'-tv'])
    if args.icTests:
        print("Running Initial Condition tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','ICtest.py'),'-tv'])
    if args.intelChecks:
        print("Running memory and thread tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','intelChecks.py'),'-tv'])
    if args.reproTests:
        print("Running reproducibility tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','mage_reproducibility_check.py'),'-tv'])

if __name__ == "__main__":
    main()

