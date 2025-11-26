GAMERA-RCM
==========

Voltron is the executable name for the standard magnetosphere simulation. It
consists of GAMERA global MHD model and RCM ring current model. GAMERA-RCM
runs as a coupled magnetosphere-ring current system. Options are available in
Voltron to not enable RCM coupling, i.e. pure MHD mode, while only keep the
field line tracing capacity in RCM.

Checklist
---------

Here is a checklist of SIX necessary ingredients you should go over before
launching a Voltron run:

#. Executable: e.g., voltron.x (serial run) or voltron_mpi.x (parallel run).
   See compilation instructions.
#. Grid file: e.g., lfmQ.h5. See script instructions for genLFM.py.
#. Solar wind file: e.g., bcwind.h5. See script instructions for omni2wind.py
   or gen_SW_kaiju.py.
#. RCM configuration file: e.g., rcmconfig.h5. See script instructions for
   genRCM.py.
#. Configuration file: e.g., cmriQ.xml. See :doc:`XML <./xml/voltronXML>`
   instructions and pay attention to the differences between serial run and
   MPI run.
#. Job submission file: e.g., RunVOLTRON.pbs. See pbs instructions and pay
   attention to the differences between serial run and MPI run.

Compilation
-----------

Assume Kaiju repository has been installed successfully at $KAIJUDIR, and
$KAIJUDIR/scripts has been added to $PATH and $KAIJUDIR added to $PYTHONPATH.

The commands to compile for MPI parallelized Voltron, which is the mostly used
mode, include:

.. code-block:: shell

   cd $KAIJUDIR
   module purge
   module restore kaiju
   rm -r build
   mkdir build
   cd build
   cmake -DENABLE_MPI=ON -DENABLE_MKL=OFF ..
   make voltron_mpi.x

Here are more explanations and exceptions for the operations above.

#. "module restore kaiju" assumes that you have saved the eight necessary
modules under your home directory. The eight necessary modules are listed
below. If you haven't saved them, do "module purge", then
"module load git/2.9.5" and one by one for each of the seven rest modules,
then "module save kaiju".

.. code-block::

      1) git/2.9.5
      2) intel/18.0.5
      3) hdf5/1.10.5
      4) impi/2018.4.274
      5) ncarenv/1.3
      6) ncarcompilers/0.5.0
      7) python/2.7.16
      8) cmake/3.14.4

#. "rm -r build" is to clean up existing build directory to avoid any residual
settings from previous compilations. Skip this if there is no build directory
under $KAIJUDIR.

#. "cmake -DENABLE_MPI=ON -DENABLE_MKL=OFF .." is for MPI parallelized run.
The flag "-DENABLE_MPI" is by default off to compile serial run. The flag
"-DENABLE_MKL" is also by default off and the gmres solver is used for the
Posisson's equation. When "-DENABLE_MKL=ON", the threaded Intel pardiso solver
is used. Recent tests show that Voltron results are not reproducible with MKL
on. Pardiso only brings a few percent improvement of running speed. It is
suggested to keep the default setting of MKL off.

#. "make voltron_mpi.x" is for MPI parallelized run. If compiling serial
Voltron, simply use "make voltron.x" after "cmake ..".

#. Additional notes from history page: The build system uses cmake which will
attempt to auto-detect HDF5/OMP/MPI settings, however optionally you can
provide a file "cmake/user.cmake" to set various variables if the auto-detect
doesn't work.

#. To check if your operations are in the right place, refer to the typical
screen outputs when compiling:

.. code-block:: bash

   [(NPL) ] kaiju/build> cmake ..
   -- The Fortran compiler identification is Intel 18.0.5.20180823
   -- Check for working Fortran compiler: /glade/u/apps/ch/opt/ncarcompilers/0.5.0/intel/18.0.5/ifort
   -- Check for working Fortran compiler: /glade/u/apps/ch/opt/ncarcompilers/0.5.0/intel/18.0.5/ifort  -- works
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Checking whether /glade/u/apps/ch/opt/ncarcompilers/0.5.0/intel/18.0.5/ifort supports Fortran 90
   -- Checking whether /glade/u/apps/ch/opt/ncarcompilers/0.5.0/intel/18.0.5/ifort supports Fortran 90 -- yes
   -- HDF5: Using hdf5 compiler wrapper for all Fortran compiling
   -- Found HDF5: Included by compiler wrappers  found components:  Fortran
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0")
   -- Found OpenMP: TRUE (found version "5.0") found components:  Fortran

   Configuration summary ...
   System: cheyenne4
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 18.0.5.20180823
   HDF5 Wrapper:
   Version: 91c9592 / master
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp

   Build Flags: -O3 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo -march=corei7 -axCORE-AVX2
   -----------------------------------------------------------------------------------------------------------------------

   Adding CHIMP module ...
      EB IC file is /glade/u/home/ldong/aplkaiju/kaiju/src/chimp/ebICs/ebICstd.F90
      TP IC file is /glade/u/home/ldong/aplkaiju/kaiju/src/chimp/tpICs/tpICstd.F90
      Adding executable project.x
      Adding executable psd.x
      Adding executable push.x
      Adding executable slice.x
      Adding executable chop.x
      Adding executable trace.x
   Adding Gamera module ...
      Bricksize is 16
      IC file is /glade/u/home/ldong/aplkaiju/kaiju/src/gamera/ICs/null.F90
      Adding executable gamera.x
   Adding ReMIX module ...
      Adding executable remix.x
   Adding RCM module ...
      Adding executable rcm.x
   Adding Voltron module ...
      IC file is /glade/u/home/ldong/aplkaiju/kaiju/src/voltron/ICs/earthcmi.F90
      Adding executable voltron.x
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /glade/u/home/ldong/aplkaiju/kaiju/build
   [(NPL) ] kaiju/build>

When cmake is ready, you can start compiling the executable with a simple make
command:

.. code-block:: bash

   make voltron

Normal outputs look like this, showing the percentage of completeness:

.. code-block:: shell

   [(NPL) ] kaiju/build> make voltron
   Scanning dependencies of target baselib
   [  1%] Building Fortran object src/base/CMakeFiles/baselib.dir/kdefs.F90.o
   ...... Lines omitted by editor of this wiki page ......
   [ 20%] Linking Fortran static library libbaselib.a
   [ 20%] Built target baselib
   Scanning dependencies of target rcmlib
   ...... Lines omitted by editor of this wiki page ......
   [ 34%] Built target rcmlib
   Scanning dependencies of target chimplib
   [ 34%] Building Fortran object src/chimp/CMakeFiles/chimplib.dir/chmpunits.F90.o
   ...... Lines omitted by editor of this wiki page ......
   [ 57%] Linking Fortran static library libchimplib.a
   [ 57%] Built target chimplib
   Scanning dependencies of target gamlib
   [ 59%] Building Fortran object src/gamera/CMakeFiles/gamlib.dir/gamutils.F90.o
   ...... Lines omitted by editor of this wiki page ......
   [ 78%] Linking Fortran static library libgamlib.a
   [ 78%] Built target gamlib
   Scanning dependencies of target remixlib
   [ 79%] Building Fortran object src/remix/CMakeFiles/remixlib.dir/mixconductance.F90.o
   ...... Lines omitted by editor of this wiki page ......
   [ 85%] Linking Fortran static library libremixlib.a
   [ 85%] Built target remixlib
   Scanning dependencies of target voltlib
   [ 85%] Building Fortran object src/voltron/CMakeFiles/voltlib.dir/ICs/earthcmi.F90.o
   ...... Lines omitted by editor of this wiki page ......
   [ 97%] Linking Fortran static library libvoltlib.a
   [ 97%] Built target voltlib
   Scanning dependencies of target voltron.x
   [ 98%] Building Fortran object CMakeFiles/voltron.x.dir/src/drivers/voltronx.F90.o
   [100%] Linking Fortran executable bin/voltron.x

and followed by hundreds of lines like this with "remark #15009", which are
good messages telling that the compiler is making the code faster (-Kareem):

.. code-block:: shell

   /glade/u/home/ldong/aplkaiju/kaiju/src/drivers/voltronx.F90(3): remark #15009: MAIN__ has been targeted for automatic cpu dispatch

The compilation for voltron.x is successful when you see this at the end:

.. code-block:: shell

   [100%] Built target voltron.x
   Scanning dependencies of target voltron
   [100%] Built target voltron

Check the running status (may move this part to another page)
-------------------------------------------------------------

Check the status of submitted job with

.. code-block:: shell

   qstat -u username

Or real time output with

.. code-block:: shell

   tail -f *.out

Quick lookup of run outcomes
----------------------------

A few python tools for diagnostics are available under $KAIJUDIR/scripts. The
msphpic.py can be used to make some combined RCM/remix/Gamera figures like
below:

.. image:: https://bitbucket.org/repo/kMoBzBp/images/1000851377-qkpic.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/1000851377-qkpic.png
   :alt: qkpic.png


The gamsphVid.py can be used to generate multiple plots for making animations.
Use "gamsphVid.py -h" for instructions on usage.

MPI Differences
---------------

Running a coupled Gamera-RCM case with MPI support requires three things:

#. Building the MPI version of the coupled executable.
#. Modifying case XML to supply additional MPI decomposition information
#. Modifying the submission script to request multiple nodes and use mpirun

Modifying Case XML
~~~~~~~~~~~~~~~~~~

Modifying a case XML for a coupled case is very similar to modifying one for
an MHD-only case. The same modifications are required to **only the Gamera**
section of the XML to define what the MPI decomposition is. Only Gamera
currently supports decomposition, so no other sections of the XML require
modification.

**Note**: currently coupled Gamera-RCM only supports MPI decomposition in the
I and J dimensions. Decomposition along the K dimension will result in errors
or bad results.

Three additional lines are required in the case XML file when running with MPI
decomposition. And these lines are ignored by the non-MPI version of coupled
Gamera, so you can safely leave them in the XML (if you want to) when not
using MPI.

In the Gamera section of the XML, one line is required
for each dimension that tells how many regions that dimension is decomposed
into, and whether that dimension is periodic. Here is an example where the
case is decomposed into 4 regions along the I and J axes, and not decomposed
along the K axis. The I and J axes are not periodic, but the K axis is.

.. code-block:: shell

   <Gamera>
   ...
      <iPdir N="4" bcPeriodic="F"/>
      <jPdir N="4" bcPeriodic="F"/>
      <kPdir N="1" bcPeriodic="T"/>
   ...
   </Gamera>

Modifying Job Submission Script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For information about creating job submission scripts, check the
:doc:`Gamerasphere <Gamerasphere>` article. Assuming that you have a job
submission script suitable for running serial jobs, here is how to modify it
to run coupled MPI gamera. These examples will be designed to work with
Cheyenne, but can be adapted to most clusters.

First, you need to request an appropriate number of nodes and tell the job
submission system how many MPI ranks should be created per node. Continuing
the example above, this case is decomposed 4 times along the I dimension, and
not at all along the J or K dimensions. So this case will need a total of
4\ *1*\ 1=4 MPI ranks **for Gamera**, and one additional MPI rank for
coupled RCM.

For this example we will assign one MPI rank to each physical processor/socket
for Gamera, which provides a reasonable balance between performance and cost.
Each of Cheyenne's compute nodes has two processors/sockets, so this means
that each compute node will receive two MPI Gamera ranks. The coupled RCM MPI
rank is more resource intensive, and so will get an entire compute node to
itself. The original, serial, resource request line from the job submission
script looked like this:

.. code-block:: shell

   #PBS -l select=1:ncpus=72:ompthreads=72

That line requests all 72 cpus on a single compute node, which is perfect for
a single process. We want to create a total of 5 processes spread across three
compute nodes. We want 2 compute nodes with 2 MPI ranks each which will give
us 4 ranks for Gamera, and then 1 compute node with 1 MPI rank by itself for
coupled RCM. That looks like this:

.. code-block:: shell

   #PBS -l select=2:ncpus=36:mpiprocs=2:ompthreads=36+1:ncpus=36:mpiprocs=1:ompthreads=72

**Note** that in larger cases you will want to add helper ranks that share some
of the workload assigned here to the coupled RCM rank. Those helper ranks will
also require entire nodes, and should be added to the latter portion of the
PBS select command. A command that adds 2 more nodes for helpers compared to
the one above (5 nodes total) would look like this:

.. code-block:: shell

   #PBS -l select=2:ncpus=36:mpiprocs=2:ompthreads=36+3:ncpus=36:mpiprocs=1:ompthreads=72

The line that controls the number of OMP threads should also be cut in half
since we now have 2 Gamera processes per node:

.. code-block:: shell

   export OMP_NUM_THREADS=36

The only other line we need to change is the one that calls the executable and
starts the simulation. In the serial case that looked like this:

.. code-block:: shell

   ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out

That command literally calls the executable and passes it the input XML file.
Instead, we now need to use a helper application called mpirun, which will
call our executable for us:

.. code-block:: shell

   mpirun ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out

Using MPT
~~~~~~~~~

If you are running with the MPT mpi library, the submission script will
require some additional modifications, described in a
:doc:`dedicated page <runningWithMPT>`.
