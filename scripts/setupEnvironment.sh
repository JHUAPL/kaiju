#!/usr/bin/env bash
# usage: setupEnvironment.sh
#
# This script will automatically set the KAIJUHOME environment variable to the root of this kaiju installation,
#   and also add the appropriate scripts folders to your PATH and PYTHONPATH environment variables

# borrowed this one-liner to get the directory containing the script:
# https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
# modified to also work with zsh
# https://stackoverflow.com/questions/9901210/bash-source0-equivalent-in-zsh
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-${(%):-%x}}" )" &> /dev/null && pwd )"

# strip off the "/scripts" folder to get to the root of the repository
ROOT_DIR="$(echo "$SCRIPT_DIR" | sed 's:/scripts$::')"

# using only export= and not export+= to work with bash and zsh
export KAIJUHOME="$ROOT_DIR"

export PATH="$PATH:$KAIJUHOME/scripts:$KAIJUHOME/scripts/preproc:$KAIJUHOME/scripts/makeitso:$KAIJUHOME/scripts/makeitso-gamhelio"
