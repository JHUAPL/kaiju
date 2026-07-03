Heliosphere Quick Start Guide
=============================

** Notes for Updating ** Want to replace this page the instructions on
steps to follow to do your own event simulation. Basic points, generate
solar wind file, decide on what components Gamera, Mix, RCM, T*GCM
you'll be using, decide on resolution, generate grid files, generate
batch script files, run the code, refer them to the post processing
steps.

**Note:** This quick start assumes you have completed the build instructions.

**Note:** Throughout the descriptions ``$KAIJUHOME`` refers to the
base directory of the `kaiju <https://bitbucket.org/aplkaiju/kaiju>`_
repository.

Initial Setup
-------------

Simple build instructions
~~~~~~~~~~~~~~~~~~~~~~~~~

To compile GAMERA Helio do the following in your ``~/kaiju/build directory``:

.. code-block:: bash

    cmake  -DENABLE_MPI=ON  ..
    make gamhelio_mpi.x

Using ``makeitso`` to create input and batch submission files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note: Supported HPC systems only - NCAR Derecho and NASA
Pleiades/Electra/Aitken (default Pleiades Broadwell nodes)

To run ``gamhelio-makeitso.py`` you must first ensure you have a proper
python environment with the prerequisite packages, add the kaipy scripts to
your path by running:

.. code-block:: bash

    . ~/kaiju/scripts/setupEnvironment.sh

you may then run the ``gamhelio-makeitso.py`` script which will guide you
through setting up the desired model and job parameters. This script operates
as a command-line input program, where basic parameter prompts will be
displayed with default settings, which may be modified.
``gamhelio-makeitso.py`` may be run in ``BASIC``, ``INTERMEDIATE``,
or ``EXPERT`` mode, with varying levels of customization of run parameters.
The command modes are available as an argument to the script. The script run
with the ``--verbose`` option will display additional status output while the
script runs. Upon completion, the script will have generated appropriate
input files and batch scripts, and finally print instructions on how to
submit the generated run to your system's batch scheduler.

.. code-block:: shell

   gamhelio-makeitso.py --verbose --mode=(BASIC|INTERMEDIATE|EXPERT)

Manual Run Setup
----------------

Set a grid and boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create files with a spherical grid and inner boundary conditions for the
inner heliosphere simulation run

.. code-block:: bash

   wsa2gamera.py ~/kaiju/kaipy/gamhelio/ConfigScripts/startup.config

If needed edit a config file
~/kaiju/kaipy/gamhelio/ConfigScripts/startup.config. Ni, Nj, Nk set a number
of cells in r, theta and phi directions, respectively. tMin and tMax set a
range in theta counting from the North (+Z) direction corresponding to
theta=0. Rin and Rout set a range in radius (distance unit is in solar radii).

.. code-block:: shell

   ;Comments and definitions:
   ;Modify if needed paths to a grid file, output innerbc file and WSA fits file
   ;tMin and tMax set a range for theta [tMin, tMax]*pi
   ;Rin and Rout are inner and outer boundaries in a radial direction
   ;Ni, Nj, Nk set number of cells in r, theta, phi directions
   ;Nghost is a number of ghost cells

   [Gamera]
   gameraGridFile = heliogrid.h5
   GridDir = ./
   gameraIbcFile = innerbc.h5
   IbcDir = ./

   [Grid]
   tMin = 0.1
   tMax = 0.9
   Rin = 21.5
   Rout = 215.
   Ni = 128
   Nj = 64
   Nk = 128

   [WSA]
   ;wsafile is the path to the WSA fits file relative to $KAIJUHOME
   ;Helio test uses WSA file for Carrington Rotation 2193
   wsafile = examples/helio/vel_201708132000R002_ahmi.fits
   density_temperature_infile = no
   gauss_smooth_width = 0 ; 8
   normalized = no

   [Constants]
   gamma = 1.5
   Nghost   = 4
   Tsolar = 27.27

   [Normalization]
   B0 = 1.e-3  ; in [Gs] equal to 100 [nT]
   n0 = 200.   ; in [cm-3]
   T0 = 1.e6 ; in [K]

By default a spherical grid for inner heliosphere simulation is uniform.
Other grid options are in $KAIJUHOME/kaipy/gamera/gamGrids.py. GenKSphNonU
creates a non-uniform grid in r-direction changing smoothly from finer grid
near the inner boundary to coarser grid near the outer boundary; GenKSphNonUG
creates a custom grid for a CME simulation with a fine uniform grid in the
region 0.1-0.3 AU and a non-uniform coarser grid further out to 1 AU. If
needed modify wsa2gamera.py to use any of these options or create your own
grid function in $KAIJUHOME/kaipy/gamera/gamGrids.py.

Check that you successfully generated heliogrid.h5 with the grid and
innerbc.h5 with boundary conditions in a run directory.

.. code-block:: shell

   HDF5 "innerbc.h5" {
   FILE_CONTENTS {
    group      /
    dataset    /br
    dataset    /br_kface
    dataset    /et_kedge
    dataset    /rho
    dataset    /temp
    dataset    /vr
    dataset    /vr_kface
    }
   }

XML input file
~~~~~~~~~~~~~~

An example wsa.xml input file for mpi gamera helio run is as follows (for
resolution NixNjxNk = 128x64x128 or 256x128x256):

.. code-block:: xml

    <?xml version="1.0"?>
    <Kaiju>
        <Gamera>
            <sim runid="wsa" doH5g="T" H5Grid="heliogrid.h5" icType="user" pdmb="1.0" rmeth="7UP"/>
            <time tFin="200."/>
            <output dtOut="50." tsOut="50" timer="F"/>
            <physics doMHD="T" gamma="1.5"/>
            <prob Tsolar = "25.38"/>
            <restart resFile = "wsa.Res.00008.h5" dtRes="1000." doRes="F"/>
            <iPdir N="4" bcPeriodic="F"/>
            <jPdir N="2" bcPeriodic="F"/>
            <kPdir N="4" bcPeriodic="T"/>
        </Gamera>
    </Kaiju>

For high-resolution run 1024x512x1024 use the following de-composition

.. code-block:: XML

   <iPdir N="8" bcPeriodic="F"/>
   <jPdir N="4" bcPeriodic="F"/>
   <kPdir N="8" bcPeriodic="T"/>

Have wsa.xml in a run directory.

PBS script
~~~~~~~~~~

Here is an example pbs script to run mpi gamera (for resolution helio run
128x64x128)

.. code-block:: bash

   #!/bin/bash
   #PBS -A UJHB0015
   #PBS -N heliompi
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=02:00:00
   #PBS -l select=16:ncpus=36:mpiprocs=2:ompthreads=36
   #PBS -m abe
   #PBS -M your_email_address

   #Example usage

   export EXE="./gamera_mpi.x"
   export RUNID="wsa"

   source ~/.bashrc

   module list
   hostname
   date
   #export OMP_NUM_THREADS=36
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   mpirun ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date

The example above uses 16 computer nodes (2 MPI ranks per node) creating 32
processes for 32 MPI ranks (4x2x4 = 32 in decomposition for low resolution
run above).

For high resolution run 1024x512x1024 we have 8x4x8 = 256 MPI ranks so we
select 128 nodes (with 2 MPI ranks/node).

.. code-block:: shell

   #PBS -l walltime=11:59:00
   #PBS -l select=128:ncpus=36:mpiprocs=2:ompthreads=36:mem=109GB

See PBS job basics `here <https://arc.ucar.edu/docs>`_.

Submitting a run
----------------

Copy or link gamera executable ~/kaiju/build/bin/gamera_mpi.x.

.. code-block:: bash

   ln -s ~/kaiju/build/bin/gamhelio_mpi.x gamhelio_mpi.x

Have in a run directory grid file heliogrid.h5, boundary conditions file
innerbc.h5, pbs script gamera.pbs and input xml file wsa.xml.

.. code-block:: shell

   user@cheyenne5:/glade/work/user/helioRun> ls
   gamera_mpi.x  gamera.pbs  heliogrid.h5  innerbc.h5  wsa.xml

Run the job

.. code-block:: bash

   qsub gamera.pbs

Check a status of your job in a queue

.. code-block:: bash

    qstat -u username

Normalization in Gamera-Helio (!Move to Model Description!)
-----------------------------------------------------------

The three main normalization parameters are

#. Length L = 1R_S = 6.955e10 cm
#. Magnetic field magnitude B0 = 100 nT = 1.e-3 Gs
#. Number density n0 = 200 cm-3

Velocity is normalized to the Alfven velocity V0 = B0/sqrt(4 pi rho0)
~ 150 km/s. Time is normalized to t = L/V0 = 4637 s ~ 1.29 h ~ 1 hr 17 min.
Pressure is normalized to the magnetic pressure B0^2/(4*pi).

Helio Test Run
--------------

From Google Doc:

* Take 128x64x128 as the baseline resolution. It should run at
  50min/CR on one Cheyenne node.
* Take CR2193 because it compares well (solar minimum). Ask Nick for
  permission to use this file as a test. Put this file somewhere where
  everything else binary is and use bitbucket LFS.
* Run the script that generates the grid and BC files.
* Make xml from ini (we provide ini file for the test).
* Compile and run the code.
* Compare the results with the quick-look plots.
