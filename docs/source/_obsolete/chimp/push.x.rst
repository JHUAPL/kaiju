Test particle pusher (push.x)
=============================

Files needed in run directory
-----------------------------

``ebfile`` - MHD data (such as msphere.gam.h5) through which to integrate test
particle trajectories through

``eRB.xml`` - xml file detailing test particle initialization and simulation
parameters

``RunCHIMP.pbs`` - job script to submit run

``push.x`` - compiled executable

No grid file is needed as CHIMP pulls this information form the MHD data file.
CHIMP is not currently set up to run simultaneously with GAMERA, therefore MHD
files are used as input to one-way drive test particle simulations with no
feedback included on the large scale electromagnetic fields.

Can further parallelize the the code by running multiple jobs (one job per
node) simultaneously. Each job uses a randomized seed based on the job number
to seed the test particles enabling the use of the test particles from each to
increase statistics.

Initialization of test particles
--------------------------------

Standard initialization creates a number of test particles (tps/Np) of species
(tps/species; see kaiju/src/chimp/chmpunits.F90@getSpecies) in the (Z=0)
plane.  Radius (in Rx), phi (azimuth in degrees), alpha (pitch angle in
degrees), energy (in keV) specify the bounds over which test particle
parameters are randomly chosen.

By default particles are all created at ``T0``, however stream/doStream can
specify the continuous creation of test particles.  If doStream is specified
then each test particle will be given a randomly assigned birth day between
time stream/min and stream/max.

Example XML file
----------------

.. code-block:: xml

   <?xml version="1.0" ?>
   <Kaiju>
     <Chimp>
       <sim runid="eRB"/>
       <time T0="36000.0" dt="10.0" tFin="86400.5"/>
       <restart doRes="T" dtRes="1800.0" resID="eRB" nRes="-1"/>
       <fields doEBFix="T" doMHD="F" ebfile="msphere" grType="LFM" isMPI="T"/>
       <parallel Ri="8" Rj="8" Rk="1"/>
       <pusher epsht="0.05" imeth="GC"/>
       <tps Np="25000" species="e"/>
       <units uid="EARTH"/>
       <output doTrc="T" doEQProj="T" dtOut="60.0" tsOut="100"/>
       <radius max="8.0" min="2.5"/>
       <phi max="360.0" min="0.0"/>
       <alpha max="177.5" min="2.5"/>
       <energy doEBInit="F" doLog="T" max="5000.0" min="50.0"/>
       <domain dtype="SPH" gtype="LFM" rmax="20.0" rgap="1.6" rmin="1.05"/>
       <tracer epsds="0.05"/>
       <stream doStream="F" max="72000.0" min="71100.0"/>
     </Chimp>
   </Kaiju>

Parameter Descriptions
----------------------

``<sim>`` (optional): Specify identifying information for this computation.

``runid`` (optional, default ``"Sim"``): String specifying an identifier for
this run of ``push.x``. A best practice is to use the ``runid`` in the name of
the XML file.

``<time>`` (optional): Specify time range and interval for magnetic field
calculation.

``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) for
ground magnetic field calculation, relative to start of simulation results
used as input.

``dt`` (optional, default ``"1.0"``): Sets subcycling time of test particles,
sets the frequency when checks particle diagnostics, i.e. if particle on open
or closed field line, etc

``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds) for
ground magnetic field calculation, relative to start of simulation results
used as input.

``<fields>`` (required): Describes the input data from a MAGE model run.

``ebfile`` (optional, default ``"ebdata.h5"``): Path to HDF5 file containing
the electric and magnetic fields computed by a MAGE model run.

``grType`` (optional, default ``"EGG"``): String specifying grid type used by
the MAGE output file. Valid values are ``"EGG"``, ``"LFM"``, ``"SPH"``. If the
string is not one of the supported grid types, the default value (``"EGG"``)
is used, and a warning message is printed.

``doEBFix`` (optional, default ``"false"``): Set to ``"true"`` to "clean" the
electric field E so that the dot product of the electric and magnetic fields
is 0. See ``ebinterp.F90``.

``doMHD`` (optional, default ``"false"``): Set to ``"true"`` to pass the full
set of magnetohydrodynamic variables to CHIMP, rather than just electric and
magnetic fields. See ``ebtypes.F90``.

``isMPI`` (optional, default ``"false"``): Set to ``"true"`` is the MAGE
results file was generated with an MPI version of the model. See
``eblCstd.F90``.

``doNumB0`` (optional, default ``"false"``): Set to ``"true"`` to numerically
compute the background magnetic field. See ``starter.F90``.

``doPureB0`` (optional, default ``"false"``): Set to ``"true"`` use the
analytical form of the dipole field, sets electric field to 0 See
``starter.F90``.

``doEBOut`` (optional, default ``"false"``): Set to ``"true"`` to output
slices of the electric and magnetic fields.See ``starter.F90``.

``<domain>`` (optional): Options for the problem domain

``dtype`` (optional, default ``"SPH"``): Domain over which to perform CHIMP
calculations, separate from grid, enables the user to perform calculation on a
subset of the grid to reduce computation where it is not needed - See
``gridloc.F90``, line 172. Valid values are ``"SPH"``, ``"LFM"``,
``"LFMCYL"``, ``"MAGE"``, ``"EGG"``, ``"ELL"``.

``gtype`` (optional, default ``"NONE"``): If set to "SPH" or "LFM", extends
domain from inner boundary of MAGE file to value set by rmin, this assumes a
dipole magnetic field in this gap region, useful for test particle runs to
allow for particles to bounce in the gap region that would be considered lost
otherwise - See ``gridloc.F90``, line 207. Valid values are ``"NONE"``,
``"SPH"``, ``"LFM"``.

``rClosed`` (optional, default set by choice of ``units/uid``): Radial value
for field line endpoint to reach to be considered closed - See
``chmpunits.F90``.

``rmax`` (optional, default computed): Maximum radius of Domain region - See
``gridloc.F90``.

``rmin`` (optional, default computed): Minimum radius of Domain region - See
``gridloc.F90``.

``rgap`` (optional, default computed): Assume Dipole field from
rmin < r < rgap, usually set to inner boundary radius of MAGE grid, used if
gtype is ``"SPH"`` or ``"LFM"`` - See ``gridloc.F90``.

``xSun`` (optional,default 20.0): if dType is "LFM" or "MAGE", the Domain
region includes all i-shells which have distances along the Earth-Sun line is
less than this value (in Re).

``xTail`` (optional,default -100.0): if dType is "LFM" or "MAGE", the Domain
region includes cells in the magnetotail up until this value (in Re).

``yzMax`` (optional,default 40.0): if dType is "LFM" or "MAGE", the Domain
region includes cells with Y and Z coordinates between +/- yzMax (in Re).

``<interp>`` (optional): Options related to interpolation

``wgt`` (optional, default ``"TSC"``): Sets 1D interpolation type. Valid
values are ``"TSC"`` (1D triangular shaped cloud), ``"LIN"`` (linear),
``"QUAD"`` (parabolic). See ``starter.F90``.

``<gciter>`` (optional): Options for guiding center computation

``maxIter`` (optional, default ``"20"``): Maximum number of iterations for
guiding center computation.

``TolGC`` (optional, default ``"1.0e-3"``): Tolerance for guiding center
computations.

``<output>`` (optional): Options related to driver output

``T0Out`` (optional, default value of ``T0`` attribute of ``<time>`` element):
Time (in seconds from ``T0``) to begin output of data to h5 files See
``starter.F90``.

``dtOut`` (optional, default ``"10.0"``): Output cadence.

``timer`` (optional, default ``"false"``): Set to ``"true"`` to turn time
flags on See ``starter.F90``, line 139.

``tsOut`` (optional, default ``"10"``): Cadence to output diagnostics to
run-log file See ``starter.F90``.

``doEQProj`` (optional, default ``"false"``): Set to ``.true.`` to include
equatorial variables, projected down to magnetic equator along field line
from particle position (i.e. Xeq,Yeq, equatorial pitch angle etc) See
``chmpdefs.F90``.

``doLL`` (optional, default ``"false"``): Set to ``.true.`` to include
latitude and longitude values of the particles position projected to the
northern hemisphere along field line.

``doTrc`` (optional, default ``"false"``): Similar to doEQProj, used in
slice.x See ``chmpdefs.F90``.

``<parallel>`` (optional): Options if ebfile was generated using an MPI
version of the code (read if fields/doMPI is set to ``"true"``, file name in
form of ebfile_Ri_Rj_Rk_i_j_k.gam.h5)

``Ri`` (optional, default ``"1"``): Number of ranks used in decomposition  of
``"i"`` dimension See iotable.F90.

``Rj`` (optional, default ``"1"``): Number of ranks used in decomposition  of
``"j"`` dimension See iotable.F90.

``Rk`` (optional, default ``"1"``): Number of ranks used in decomposition  of
``"k"`` dimension See iotable.F90.

``doOldNaming`` (optional, default ``"false"``): Allow for backward
compatibility for MHD files generated with the now deprecated naming
convention See ``chmpdefs.F90``.

``<plasmapause>`` (optional): Options for calculation to determine plasmapause
location in MHD file

``doPP`` (optional, default ``"false"``): Set to ``.true.`` to calculate
plasmapause location in the equator and include it in output file for slice.x
- See chmpfields.F90.

``Lin`` (optional, default ``"2.0"``): Minimum L-shell value to begin
plasmapause calculation - See plasmaputils.F90.

``Lout`` (optional, default ``"6.0"``):  Maximum L-shell value to end
plasmapause calculation - See plasmaputils.F90.

``Nl`` (optional, default ``"80"``): Number of cells/steps in L-shell See
plasmaputils.F90.

``Nphi`` (optional, default ``"72"``): Number of cells/steps in longitude See
plasmaputils.F90.

``phi0`` (optional, default ``"0.0"``): Minimum longitude to perform
plasmapause calculation between (in radians, 0 is in the +X direction) -
See plasmaputils.F90.

``phi1`` (optional, default ``"2*PI"``): Maximum longitude to perform
plasmapause calculation between See plasmaputils.F90.

``<pusher>`` (optional): Options related to test particle integrator

``imeth`` (optional, default ``"FO"``): Select the integrator used to evolve
test particles. Valid values are ``"FO"``, ``"GC"``, or ``"DYN"`` for full
orbit/relativistic Lorentz equations, guiding center, or dynamic switching
between the two when appropriate.

``epsht`` (optional, default ``"0.05"``): Parameter used in timestep
calculation.

``epsgc`` (optional, default ``"0.05"``): Maximum value of adiabaticity
parameter used to determine when to switch between GC and FO integration.

``<restart>`` (optional): Simulation restart options.

``doRes`` (optional, default ``"false"``): Set to ``"true"`` to read a
restart file.

``resID`` (optional, default value of ``runid`` attribute of ``<sim>``
element): String specifying an identifier for the runID of the files to be
used to restart the test particle run See particleio.F90.

``nRes`` (optional, default ``"-1"``): Restart file number to read in, a value
of -1 reads the symlinked file resID.Res.XXXXX.h5 See chmpio.F90.

``dtRes`` (optional, default ``"-1"``): Cadence to produce restart files

``<tps>`` (optional): Options related to test particle

``Np`` (optional, default ``"100"``): Number of test particles to simulate per
job

``species`` (optional, default ``"X"``): Desired species to simulate. For full
list of available species, see the getSpecies() function in chmpunits.F90.

``<energy>`` (optional): Options for initialization of test particle energies

``min`` (optional, default ``"1.0"``): Minimum energy (in keV).

``max`` (optional, default ``"100.0"``): Maximum energy (in keV).

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute particles in
uniformly in log-space between min/max values.

``doEBInit`` (optional, default ``"false"``): Set to ``"true"`` to have
specified energy in ExB frame instead of lab frame.

``<alpha>`` (optional): Options for initialization of test particle pitch
angles

``min`` (optional, default ``"0.0"``): Minimum pitch angle (in degrees).

``max`` (optional, default ``"360.0"``): Maximum pitch angle (in degrees).

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute particles in
uniformly in log-space between min/max values.

``<psi>`` (optional): Options for initialization for gyrophase of test
particle

``min`` (optional, default ``"0.0"``): Minimum gyrophase angle (in degrees).

``max`` (optional, default ``"360.0"``): Maximum gyrophase angle (in degrees).

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute particles in
uniformly in log-space between min/max values.

``<radius>`` (optional): Options for initialization of test particle locations

``min`` (optional, default ``"5.0"``): Minimum radius particles can be
initialized at.

``max`` (optional, default ``"25.0"``): Maximum radius particles can be
initialized at.

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute particles in
uniformly in log-space between min/max values.

``<phi>`` (optional): Options for initialization of test particle locations

``min`` (optional, default ``"0.0"``): Minimum longitude (in degrees)
particles can be initialized at. A value of 0 corresponds to the +X direction.

``max`` (optional, default ``"360.0"``): Maximum longitude (in degrees)
particles can be initialized at.

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute particles in
uniformly in log-space between min/max values.

``<height>`` (optional): Options for initialization of test particle locations

``min`` (optional, default ``"0.0"``): Minimum value for "height" particles
can be initialized at. Default is height above Z=0 plane.

``max`` (optional, default ``"0.0"``): Maximum value for "height" particles
can be initialized at. Default is height above Z=0 plane.

``doLog`` (optional, default ``"false"``): Set to ``"true"`` to distribute
particles in uniformly in log-space between min/max values.

``doOutflow`` (optional, default ``"false"``): Set to ``"true"`` for height to
be treated as a latitude (in degrees), useful for tracing particle outflow
from inner MHD boundary - See ``tpICstd.F90``.

``doWind`` (optional, default ``"false"``): Set to ``"true"`` to keep
particles between Z = +/- zWind. Only valid for full orbit particles. See
``tpICstd.F90``.

``zWind`` (optional, default ``"0.0"``): Value to constrain test particles
between to keep particles between Z = +/- zWind See ``tpICstd.F90``.

``<stream>`` (optional): Options related to time dependent initialization of
test particles

``doStream`` (optional, default ``"false"``): Set to ``"true"`` to inject test
particles over time instead of all being injected at ``T0`` See chmpdefs.F90.

``min`` (optional, default ``"0.0"``): Minimum time (in seconds from ``T0``)
to begin streaming particles.

``max`` (optional, default ``"0.0"``): Maximum time (in seconds from ``T0``)
to stop streaming particles. Particles will be uniformly released between
min/max times.

``<tracer>`` (optional): Options related to field line tracing performed by
CHIMP

``epsds`` (optional, default ``"1.0e-2"``): Tolerance for field line tracing
computations See chmpdefs.F90.

``<units>`` (optional): Name of units system used in the model run.

``uID`` (optional, default ``"Earth"``): See chmpunits.F90 line 148. Valid
values are ``"EARTH"``, ``"EARTHCODE"``, ``"JUPITER"``, ``"JUPITERCODE"``
``"SATURN"`` ``"SATURNCODE"``\ , ``"HELIO"``\ , ``"LFM"``\ , ``"LFMJUPITER"``.

Run Script example
------------------

An example pbs script on ``cheyenne``\ , RunCHIMP.pbs to submit test particle
simulations:

.. code-block:: bash

   #!/bin/bash
   #PBS -A P28100045
   #PBS -N CHIMP
   #PBS -j oe
   #PBS -q economy
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72

   #Example usage
   #qsub -v CHIMPEXE="pusher.x" -J 1-5 -N RC_ep RunCHIMP.pbs

   export EXE=${CHIMPEXE:-"./push.x"}
   export RUNID=${PBS_JOBNAME}

   source ~/.bashrc
   module restore kaiju

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date

This command submits a batch job of 20 test particle runs

.. code-block:: bash

   qsub -J 1-20 -N eRB RunCHIMP.pbs

Output
------

Creates runID.job#.h5part files (one for each job, holding the particles it
simulated). These files can be read in directly into Visit or Paraview. It is
beneficial to also extract a 2D plane from the equator of the MHD data (see
slice.x documentation) and visualize particle location and MHD solution
together for context.
