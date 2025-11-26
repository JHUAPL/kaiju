CHIMP XML
=========

``push.x``
----------

.. code-block:: xml

   <?xml version="1.0" ?>
   <Chimp>
       <sim runid="erc"/>
       <time T0="1290.0" dt="0.05" tFin="1310.0"/>
       <fields ebfile="ebBBFs.h5" grType="LFM"/>
       <pusher epsht="0.05" imeth="DYN"/>
       <tps Np="500000" species="e"/>
       <units uid="EARTH"/>
       <output T0Out="1300.0" doEQProj="T" dtOut="0.1" tsOut="100"/>
       <radius max="20.0" min="10.0"/>
       <phi max="150.0" min="105.0"/>
       <alpha max="180.0" min="0.0"/>
       <energy doLog="F" max="50.0" min="0.1"/>
       <tracer epsds="0.025"/>
       <stream doStream="T" max="1300.0" min="1290.0"/>
       <parallel Ri="4" Rj="4" Rk="1"/>
   </Chimp>

The XML elements are:

``<Chimp>``
    This element is required.

``<sim>``
    ``runid`` The name of the simulation.

``<time>``
    ``T0`` Start time for the particle trajectory (seconds).

    ``dt`` Time step for the particle trajectory (seconds).

    ``tFin`` End time for the particle trajectory (seconds).

``<fields>``
    ``ebfile`` The name of the file containing the electromagnetic fields.

    ``grType`` The type of grid used for the electromagnetic fields.

``<pusher>``
    ``epsht`` Small parameter for the integrator.

    ``imeth`` Integration method, FO (full orbit), GC (guiding center), DYN
    (switch between both).

``<tps>``
    ``Np`` Number of test particles.

    ``species`` Species of the test particles.

``<units>``
    ``uid`` Which units to use for reading input data. The default value is
    "Earth" for using GAMERA outputs of magnetosphere electromagnetic fields.
    CHIMP also works with LFM outputs by changing the uid to "LFM". Don't
    misuse or you will get absurdly wrong results and crazily slow
    performance. See more details in kaiju/src/chimp/chmpunits.F90.

``<output>``
    ``T0Out`` When to start outputting.

    ``doEQProj`` Whether to project the test particle positions to the
    equatorial plane.

    ``dtOut`` Time interval for output files.

    ``tsOut`` Timestep interval for console output.

``<radius>``
    ``max`` Maximum radius for the test particles.

    ``min`` Minimum radius for the test particles.

``<phi>``
    ``max`` Maximum azimuth for the test particles.

    ``min`` Minimum azimuth for the test particles.

``<alpha>``
    ``max`` Maximum pitch angle for the test particles.

    ``min`` Minimum pitch angle for the test particles.

``<energy>``
    ``doLog`` Whether to use a logarithmic distribution for the test particle
    energies.

    ``max`` Maximum energy for the test particles.

    ``min`` Minimum energy for the test particles.

``<tracer>``
    ``epsds`` Small parameter for the field line tracer.

``<stream>``
    ``doStream`` Whether to continuously create test particles.

    ``max`` Maximum time for the test particle creation.

    ``min`` Minimum time for the test particle creation.

``<parallel>``
    ``Ri`` Number of processors in the i-direction.

    ``Rj`` Number of processors in the j-direction.

    ``Rk`` Number of processors in the k-direction.

Initialization of test particles
--------------------------------

Standard initialization creates a number of test particles (tps/Np) of species
(tps/species; see kaiju/src/chimp/chmpunits.F90@getSpecies) in the (Z=0)
plane.  Radius (in Rx), phi (azimuth in degrees), alpha (pitch angle in
degrees), energy (in keV) specify the bounds over which test particle
parameters are randomly chosen.

By default particles are all created at T0, however stream/doStream can
specify the continuous creation of test particles.  If doStream is specified
then each test particle will be given a randomly assigned birth day between
time stream/min and stream/max.

Gamera output slicer (slice.x)
------------------------------

Example XML file

.. code-block:: xml

   <?xml version="1.0" ?>
   <Chimp>
       <sim runid="ebXY"/>
       <time T0="36000.0" dt="5.0" tFin="43200.0"/>
       <fields doMHD="T" ebfile="msphere" grType="LFM" isMPI="T"/>
       <parallel Ri="6" Rj="12" Rk="1"/>
       <domain dtype="LFMCYL"/>
       <units uid="EARTH"/>
       <slice Npow="1" doXY="T" grType="LFM2D" xSun="25.0"/>
       <tracer epsds="0.05"/>
       <output doSlim="T" doTrc="T"/>
   </Chimp>

An example pbs script to submit slice job:

.. code-block:: bash

   #!/bin/bash
   #PBS -A P28100045
   #PBS -N sliceXZ
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=4:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72

   #Example usage
   #qsub -v KAIJUEXE="./pusher.x" -J 1-5 -N RC_ep RunK.pbs
   #Module savelist "kaiju"
   #Currently Loaded Modules:
   #  1) git/2.9.5    (H)   4) impi/2018.4.274       7) python/2.7.16
   #  2) intel/18.0.5       5) ncarenv/1.3           8) cmake/3.14.4
   #  3) hdf5/1.10.5        6) ncarcompilers/0.5.0

   export EXE=${slice}
   export RUNID=${PBS_JOBNAME}

   ###source ~/.bashrc
   module restore kaiju

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   echo "Running $EXE"
   ./slicer.x ${RUNID}.xml > ${RUNID}.out
   date

Gamera field line tracer (trace.x)
----------------------------------

Example XML file

.. code-block:: xml

   #!XML
   <?xml version="1.0" ?>
   <KAIJU>
       <Chimp>
           <sim runid="IonFL"/>
           <time T0="36000.0" dt="5.0" tFin="43200.0"/>
           <fields doMHD="T" ebfile="msphere" grType="LFM" isMPI="T"/>
           <parallel Ri="6" Rj="12" Rk="1"/>
           <domain dtype="LFMCYL"/>
           <units uid="EARTH"/>
           <interp wgt="TSC"/>
           <tracer epsds="0.05"/>
           <output dtOut="5.0"/>
           <points Nx1="1" Nx2="72" Nx3="16" grType="SPHERICAL"/>
           <radius min="2.05" max="2.05"/>
           <phi min="180.0" max="360.0"/>
           <theta min="5.0" max="45.0"/>
       </Chimp>
   </KAIJU>

This will create an HDF5 file with a similar group structure for time
slices (Step#0,Step#1,xxx) but with each time slice also containing a group
for each field line seed, i.e.

*
  Step#0


  * Line#0
  * Line#1

*
  Step#1


  * Line#0
  * Line#1

In the above example, the field line seed points are set in a spherical
coordinate system with Nx1,Nx2,Nx3 points in r,phi,theta.

The script kaiju/scripts/genXLine.py can be run on the output HDF5 data to
create an XDMF file that can be read by VisIt or Paraview.
