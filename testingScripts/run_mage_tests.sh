#!/usr/bin/bash

# ############################################################################

# IMPORTANT NOTES:

# This bash script was designed to run from a cron job on the 'cron' host at
# NCAR, using ssh to execute this script on a derecho login node. See
# kaiju/testingScripts/crontab for an example of how to invoke this script.

# SSH must be configured so that ssh from cron to derecho does not require a
# password.

# ############################################################################

echo '***********************************************************************'
echo "Starting $0 at `date` on `hostname`."

# ############################################################################

# Force this script to exit on any failure.
set -e

# ----------------------------------------------------------------------------

# These *should* be the only variables you have to change when moving the
# testing environment around or changing python environments.

# Root of kaiju testing environment - *everything* goes under here, except
# for the outputs from the individual test runs.
export MAGE_TEST_ROOT='/glade/campaign/hao/msphere/automated_kaiju_tests'

# Root of kaiju test results tree.
export MAGE_TEST_RESULTS_ROOT='/glade/derecho/scratch/ewinter/mage_testing'

# Location of the miniconda installation used for testing.
export MAGE_MINICONDA='/glade/u/home/ewinter/miniconda3'

# Setup command for conda.
__conda_setup="$('/glade/u/home/ewinter/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"

# conda environment for testing
export CONDA_ENVIRONMENT='kaiju-3.12-testing'

# ----------------------------------------------------------------------------

# Other exported environment variables

# Path to directory containing the dateime-stamped directories for individual
# sets of test runs.
export MAGE_TEST_RUNS_ROOT="${MAGE_TEST_RESULTS_ROOT}/test_runs"

# PBS account to use for running tests on derecho
export DERECHO_TESTING_ACCOUNT='P28100045'

# PBS queue to use for running tests on derecho
export DERECHO_TESTING_QUEUE='main'

# PBS priority to use for running tests on derecho
export DERECHO_TESTING_PRIORITY='economy'

# Set the token for sending messages to Slack.
export SLACK_BOT_TOKEN=`cat $HOME/.ssh/slack.txt`

# IMPORTANT: Set this environment variable to force the python print()
# function to automatically flush its output in all of the testing scripts.
# This will ensure that output from the testing scripts is logged in the order
# that it is created.
export PYTHONUNBUFFERED='TRUE'

# ----------------------------------------------------------------------------

# Non-exported variables used by this script

# Path to the SSH key to use for running the tests. The key must have no
# passphrase. This is needed to allow passwordless access to BitBucket.
ssh_key_for_testing='/glade/u/home/ewinter/.ssh/id_rsa_kaiju_testing'

# Setup script for CDF code
cdf_setup_script="${MAGE_TEST_ROOT}/local/cdf/3.9.0/bin/definitions.B"

# Address of kaiju-private repository on BitBucket.
kaiju_repository='git@bitbucket.org:aplkaiju/kaiju-private.git'

# Name of local directory containing clone of repository.
local_kaiju_name='kaiju'

# Default kaiju code branch to test if no branch, commit, or tag is specified.
default_branch_to_test='development'

# Output string to separate results from different tests.
test_output_separator='------------------------------------------------------'

# ############################################################################

# Define the command-line help function.
Help()
{
    # Display Help
    echo 'Control script for running MAGE tests via cron and ssh.'
    echo
    echo "Syntax: run_mage_tests.sh [-b branch] [-c commit] [-d] [-h] [-v] 'test1[,test2,...]'"
    echo 'options:'
    echo 'b branch   Run tests using this branch'
    echo 'c commit   Run tests using this commit (or tag)'
    echo 'd          Use debug mode.'
    echo 'h          Print this help message.'
    echo 'v          Use verbose mode.'
    echo
    echo 'If no branch or commit is specified, the latest commit on the development branch will be used for the tests.'
    echo "Each test can have its own options, e.g. 'buildTest.py -d,unitTest.py -lv'"
}

# Process command-line options.
branch=''
commit=''
debug=false
verbose=false
while getopts ':b:c:dhv' option; do
    case $option in
        b) # Test a specific branch
            branch=$OPTARG;;
        c) # Test a specific commit or tag
            commit=$OPTARG;;
        d) # debug mode
            debug=true;;
        h) # display Help
            Help
            exit;;
        v) # verbose mode
            verbose=true;;
        \?) # Invalid option
            echo 'Error: Invalid option'
            exit 1;;
    esac
done

if $debug; then
    echo "branch=${branch}"
    echo "commit=${commit}"
    echo "debug=${debug}"
    echo "help=${help}"
    echo "verbose=${verbose}"
fi

# Fetch the branch or commit to test. If neither specified, use the default.
if [[ -n $branch && -n $commit ]]; then
    echo 'Cannot specify branch and commit together!'
    exit 1
fi
if [[ -z $branch && -z $commit ]]; then
    branch=$default_branch_to_test
fi

# At this point, either branch is specified, or commit/tag is specified.
# There should be no case where both are unspecified, or both are specified.
if [[ -n $branch ]]; then
    export BRANCH_OR_COMMIT=$branch
else
    export BRANCH_OR_COMMIT=$commit
fi
if $debug; then
    echo "BRANCH_OR_COMMIT=${BRANCH_OR_COMMIT}"
fi

# Fetch the list of tests to run.
tests_to_run_str=$BASH_ARGV
if $debug; then
    echo "tests_to_run_str=${tests_to_run_str}"
fi

# Split the test list string into an array.
# Tests are separated by comma, no spaces.
IFS="," tests_to_run=($tests_to_run_str)
if $debug; then
    echo 'The test commands to run are:'
    for test_cmd in "${tests_to_run[@]}"
    do
        echo $test_cmd
    done
fi

# ############################################################################

# Load the SSH private key for testing. This is needed for BitBucket access.
# Note that this key does not have a passphrase, so it can easily be used
# from a cron job.
if $verbose; then
    echo 'Loading SSH key into key agent.'
fi
eval `ssh-agent -s`
ssh-add $ssh_key_for_testing

# ############################################################################

# Set up the module system (needed when running from cron, harmless
# otherwise).
if $verbose; then
    echo 'Setting up module system.'
fi
source /etc/profile.d/z00_modules.sh

# List at-start modules.
# The 2>&1 (no spaces!) is needed since the 'module' command sends
# output to stderr by default, and so it is lost when stdout is sent
# back to the user.
if $verbose; then
    echo 'At start, the loaded modules are:'
    module list 2>&1
fi

# ############################################################################

# Activate the conda installation for MAGE testing.
if $verbose; then
    echo "Setting up conda environment for MAGE testing ${CONDA_ENVIRONMENT}."
fi

# This code is based on the code that the miniconda installer puts into
# ~/.bashrc or ~/.bash_profile when a user installs miniconda.
# This code is needed since $HOME/.bashrc is not run.
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "${MAGE_MINICONDA}/etc/profile.d/conda.sh" ]; then
        . "${MAGE_MINICONDA}/etc/profile.d/conda.sh"
    else
        export PATH="${MAGE_MINICONDA}/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate $CONDA_ENVIRONMENT
if $verbose; then
    echo "conda environment is `echo $CONDA_PREFIX`"
fi

if $verbose; then
    echo 'The installed version of kaipy in this environment is:'
    conda list kaipy
fi

# ############################################################################

# Make the CDF library available (needed by the satellite comparison
# python scripts).
if $verbose; then
    echo 'Sourcing CDF setup script.'
fi
source $cdf_setup_script

# ############################################################################

# Create the ISO 8601 datetime stamp for this set of tests.
testing_datetime=`date --iso-8601=seconds`
if $debug; then
    echo "testing_datetime=${testing_datetime}"
fi

# <HACK>
# The make tool cannot handle file or directory names that contain a colon
# (':'). Process the ISO 8601 date time string from the ISO 8601 form:
# YYYY-MM-DDTHH:mm:ss[+-]hh:mm
# into the more compact form:
# YYYYMMDD_HHMMSS
# The sed commands to do this are, in order:
# s/.\{6\}$// - Remove the time zone offset (last 6 characters)
# s/-//g - Delete all '-'.
# s/\://g - Delete all ':'.
# s/T/_/ - Convert 'T' separating date and time to underscore.
testing_datetime=`echo $testing_datetime | sed -e 's/.\{6\}$//' -e 's/-//g' -e 's/\://g' -e 's/T/_/'`
if $debug; then
    echo "testing_datetime=${testing_datetime}"
fi
# </HACK>

# Create a directory to hold all of the tests for this set of tests.
test_set_dir="${testing_datetime}-${BRANCH_OR_COMMIT}"
if $debug; then
    echo "test_set_dir=${test_set_dir}"
fi
export MAGE_TEST_SET_ROOT="${MAGE_TEST_RUNS_ROOT}/${test_set_dir}"
if $verbose; then
    echo "Creating directory for this set of tests at ${MAGE_TEST_SET_ROOT}."
fi
mkdir $MAGE_TEST_SET_ROOT

# ############################################################################

# Move to the directory for this set of tests.
cd $MAGE_TEST_SET_ROOT

# Clone the kaiju-private repository, call it kaiju locally.
if $verbose; then
    echo "Cloning repository ${kaiju_repository}."
fi
git clone $kaiju_repository $local_kaiju_name

# Move into the repository clone.
cd $local_kaiju_name

# If a branch was requested, switch to that branch. Otherwise, if a commit
# or tag was requested, check out that commit or tag.
if $verbose; then
    echo "Checking out branch/commit/tag ${BRANCH_OR_COMMIT}."
fi
git checkout $BRANCH_OR_COMMIT

# ############################################################################

# <HACK>
# Back up existing test code, since we need to use test code that is under
# development.
mkdir TEST_CODE_BACKUP
mv testingScripts TEST_CODE_BACKUP/
mv tests TEST_CODE_BACKUP/

# Copy latest test code.
new_test_code_root="${MAGE_TEST_ROOT}/kaiju-private"
if $verbose; then
    echo "Copying updated test files from ${new_test_code_root}."
fi
cp -rp $new_test_code_root/testingScripts ./testingScripts
cp -rp $new_test_code_root/tests ./tests
# </HACK>

# ############################################################################

# Set up the kaiju environment using the just-checked-out code.
kaiju_setup_script="${MAGE_TEST_SET_ROOT}/${local_kaiju_name}/scripts/setupEnvironment.sh"
if $verbose; then
    echo "Sourcing kaiju setup script ${kaiju_setup_script}."
fi
source $kaiju_setup_script
if $verbose; then
    echo "KAIJUHOME is ${KAIJUHOME}."
fi

# ############################################################################

# List at-start environment variables.
if $debug; then
    echo 'At start, the environment variables are:'
    printenv
fi

# ############################################################################

# Compute the location of the test scripts.
kaiju_test_scripts_dir="${KAIJUHOME}/testingScripts"
if $debug; then
    echo "kaiju_test_scripts_dir=${kaiju_test_scripts_dir}"
fi

# Run each test.
for mage_test in ${tests_to_run[@]}
do
    echo $test_output_separator
    if $verbose; then
        echo "Moving to ${MAGE_TEST_SET_ROOT}."
    fi
    cd $MAGE_TEST_SET_ROOT
    cmd="python ${kaiju_test_scripts_dir}/${mage_test}"
    if $verbose; then
        echo "Running test '${cmd}' at `date` on `hostname`."
    fi
    eval "${cmd}"
done
echo $test_output_separator

# ############################################################################

# Shut down the SSH agent.
if $verbose; then
    echo "Shutting down SSH key agent."
fi
ssh-agent -k

# ############################################################################

echo "Ending $0 at `date` on `hostname`."
echo '***********************************************************************'
