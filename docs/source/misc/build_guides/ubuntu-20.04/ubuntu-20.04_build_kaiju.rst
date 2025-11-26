Compiling the ``kaiju`` software on Ubuntu 20.04
====================================================

Introduction
------------

These instructions will walk you through the process of building and installing the ``kaiju`` software on a Ubuntu 20.04 system.

These instructions assume that the user is using the ``bash`` shell, and that no modifications have been made to the user "dotfiles" (``$HOME/.bashrc``, ``$HOME/.bash_profile``). If you have customized either of these files for your account, please carefully save and inspect the output from each command in the build process to ensure that no unexpected problems have crept in. To facilitate this practice, all of the commands shown below will illustrate how to save command output, and how to measure how long each step takes. The latter is a useful bit of information which can help identify build problems early in the process, avoiding much wasted time and effort later.

In the instructions below, all code will be built and installed in a Ubuntu 20.04-specific subdirectory of the user home directory, i.e. ``$HOME/ubuntu-20.04``. This particular organization is not required - it is intended as an example of one possible way to segregate software that has been built for multiple systems.

Building the serial version of the ``kaiju`` software on Ubuntu 20.04
-------------------------------------------------------------------------

Step 1: Configure build tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin by configuring tools and libraries needed for the build:

.. code-block:: shell

   # Add the Intel compilers to PATH.
   source /opt/intel/oneapi/setvars.sh

   # Set variables for HDF5.
   export HOST_SYSTEM=ubuntu-20.04
   export HDF5_DIR=$HOME/$HOST_SYSTEM/local/hdf5/1.14.1-2
   export PATH=$HDF5_DIR/bin:$PATH
   export HDF5_INCLUDE_DIRS=$HDF5_DIR/include
   export HDF5_LIBRARIES=$HDF5_DIR/lib
   export CPATH=$HDF5_INCLUDE_DIRS
   export INCLUDE="-I$HDF5_INCLUDE_DIRS"

   # Configure CDF.
   source $HOME/$HOST_SYSTEM/local/cdf/3.9.0/bin/definitions.B


Step 2: Create the build directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a system-specific build directory.

.. code-block:: shell

   cd $HOME
   mkdir HOST_SYSTEM
   cd $HOST_SYSTEM


Then make an additional subdirectory level for the branch of the code you are building (the ``development`` branch is used as an example). This arrangement is useful when you need to maintain simultaneous builds of different branches.

.. code-block:: shell

   export KAIJU_BRANCH_NAME=development
   mkdir $KAIJU_BRANCH_NAME
   cd $KAIJU_BRANCH_NAME


Step 3: Clone the ``kaiju`` repository from BitBucket
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE: This step assumes you have been granted access to the ``kaiju`` repository on BitBucket, and that you have configured an SSH key pair for use with BitBucket. If you need help with these tasks, please contact a CGS team member for assistance.

Clone the ``kaiju`` repository (or "repo") from BitBucket:

.. code-block:: shell

   git clone git@bitbucket.org:aplkaiju/kaiju.git


This process should take a minute or so. When complete, verify that the ``kaiju`` code exists in your directory (the actual directory contents may differ slightly from what is shown below):

.. code-block:: shell

   ls kaiju
   analysis  cmake  CMakeLists.txt  examples  external  gitHookScripts  kaiju.sublime-project  kaipy  places  pytests  quickstart  README.md  scripts  setup.py  src  testingScripts  tests  xml


Now move down into the cloned repo, and switch to the branch of the code you wish to use. By default, the cloned repository provides the ``master`` branch, but we want the ``development`` branch:

.. code-block:: shell

   cd kaiju
   git switch $KAIJU_BRANCH_NAME


Step 4: Run ``cmake`` to create the ``Makefile`` needed to build the software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since the ``kaiju`` code can be built in serial and MPI forms, we first make a directory in which to build the serial version of the code (use whatever name you prefer, but ``build_serial`` is simple and unambiguous):

.. code-block:: shell

   export KAIJU_BUILD_NAME=build_serial
   KAIJU_BUILD_PATH=$KAIJU_HOME/$KAIJU_BUILD_NAME
   mkdir -p $KAIJU_BUILD_PATH
   cd $KAIJU_BUILD_PATH


Now run the ``cmake`` command. Save the ``cmake`` output, and use timestamps for each step. The options shown below direct the build process to use a recent version of the Intel Fortran compiler:

.. code-block:: shell

   date; time FC=`which ifort` cmake -DALLOW_INVALID_COMPILERS=ON .. >& cmake.out; date


This command usually takes 2-3 seconds, depending on system activity. Examine the output file ``cmake.out`` for problems. It *should* look something like this:

.. code-block::

   -- The Fortran compiler identification is Intel 20.2.9.20230302
   -- Check for working Fortran compiler: /opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/ifort
   -- Check for working Fortran compiler: /opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/ifort  -- works
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Checking whether /opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/ifort supports Fortran 90
   -- Checking whether /opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/ifort supports Fortran 90 -- yes
   -- HDF5: Using hdf5 compiler wrapper to determine Fortran configuration
   -- Found HDF5: /home/ewinter/ubuntu-20.04/local/hdf5/1.14.1-2/lib/libhdf5_fortran.so;/home/ewinter/ubuntu-20.04/local/hdf5/1.14.1-2/lib/libhdf5.so;/usr/lib/x86_64-linux-gnu/libz.so;/usr/lib/x86_64-linux-gnu/libdl.so;/usr/lib/x86_64-linux-gnu/libm.so (found version "1.14.1-2") found components: Fortran 
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0") 
   -- Found OpenMP: TRUE (found version "5.0") found components: Fortran 
   CMake Warning at cmake/compilers.cmake:61 (message):
     Setting default optimization to O2 to avoid certain Intel compiler bugs
   Call Stack (most recent call first):
     CMakeLists.txt:99 (include)


   -------------------------
   Configuration summary ...
   System: ubuntu-20
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 20.2.9.20230302
             /opt/intel/oneapi/compiler/2023.1.0/linux/bin/intel64/ifort
   HDF5 Wrapper: /home/ewinter/ubuntu-20.04/local/hdf5/1.14.1-2/bin/h5fc
   Version: 3f4e147c / development
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp
   Build Flags: -O2 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo
   -------------------------

   Adding CHIMP module ...
       EB IC file is /home/ewinter/ubuntu-20.04/cgs/kaiju/development/kaiju/src/chimp/ebICs/ebICstd.F90
       TP IC file is /home/ewinter/ubuntu-20.04/cgs/kaiju/development/kaiju/src/chimp/tpICs/tpICstd.F90
       Adding executable project.x
       Adding executable psd.x
       Adding executable push.x
       Adding executable slice.x
       Adding executable chop.x
       Adding executable trace.x
       Adding executable sctrack.x
       Adding executable calcdb.x
       Adding executable wpicheck.x
   Adding Gamera module ...
       Bricksize is 16
       IC file is /home/ewinter/ubuntu-20.04/cgs/kaiju/development/kaiju/src/gamera/ICs/null.F90
       Adding executable gamera.x
   Adding Gamera Helio module ...
       IC file is /home/ewinter/ubuntu-20.04/cgs/kaiju/development/kaiju/src/gamera/ICs/helio/wsa.F90
       Adding executable gamhelio.x
   Adding ReMIX module ...
       Adding executable remix.x
   Adding RCM module ...
       RCM Grid is of size 180 x 361 x 160
       Adding executable rcm.x
   Adding Voltron module ...
       IC file is /home/ewinter/ubuntu-20.04/cgs/kaiju/development/kaiju/src/voltron/ICs/earthcmi.F90
       Adding executable voltron.x
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /home/ewinter/ubuntu-20.04/cgs/kaiju/development/kaiju/build_serial

Step 5: Compile the ``kaiju`` software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now use ``make`` to build the ``kaiju`` software, time-stamping and saving the output:

.. code-block:: shell

   date; time make >& make.out; date


This command should complete in about 9 minutes on a Ubuntu 20.04 system running on an i9 processor. When the command is finished, check the output file ``make.out``. The file is long, but the last few lines should look something like this:

.. code-block::

   [ 99%] Built target voltron
   Scanning dependencies of target remix.x
   [100%] Building Fortran object CMakeFiles/remix.x.dir/src/drivers/remix.F90.o
   [100%] Linking Fortran executable bin/remix.x
   [100%] Built target remix.x
   Scanning dependencies of target remix
   [100%] Built target remix
   Scanning dependencies of target gamera
   [100%] Built target gamera

To verify that all of the ``kaiju`` programs have been built, examine the ``bin`` subdirectory of your ``build_serial`` directory (this list will evolve as more programs are added):

.. code-block:: shell

   ls bin
   calcdb.x  gamhelio.x  psd.x   remix2rcm.x    sctrack.x  voltron.x
   chop.x    kaitoy.x    push.x  remix2remix.x  slice.x    wpicheck.x
   gamera.x  project.x   rcm.x   remix.x        trace.x


Using the ``kaiju`` software
--------------------------------

Once built, you must run the setup script before using the ``kaiju`` software:

.. code-block:: shell

   source $KAIJU_HOME/scripts/setupEnvironment.sh


This script will set environment variables needed by the ``kaiju`` software, including the ``KAIJUHOME`` environment variable (not the ``KAIJU_HOME`` environment variable). However, the path to the compiled programs is not added - you will need to specify the complete path when using compiled programs. For example,  to run the serial version of ``gamera.x``:

.. code-block:: shell

   $KAIJUHOME/build_serial/bin/gamera.x
