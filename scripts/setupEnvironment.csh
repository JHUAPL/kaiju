#!/usr/bin/env csh
# usage: source setupEnvironment.csh
#
# This script will automatically set the KAIJUHOME environment variable to the root of this kaiju installation,
#   and also add the appropriate scripts folders to your PATH and PYTHONPATH environment variables
 
# borrowed this one-liner to get the directory containing the script:
# https://stackoverflow.com/questions/59895/how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itsel
 
set rootdir = `dirname $0`
set SCRIPT_DIR = `cd $rootdir && pwd`
 
# strip off the "/scripts" folder to get to the root of the repository
set ROOT_DIR = `cd $rootdir/.. && pwd`
  
setenv KAIJUHOME $ROOT_DIR
  
setenv PATH ${PATH}:${KAIJUHOME}/scripts:${KAIJUHOME}/scripts/preproc:${KAIJUHOME}/scripts/makeitso:${KAIJUHOME}/scripts/makeitso-gamhelio
 