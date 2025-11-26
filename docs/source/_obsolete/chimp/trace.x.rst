Executable to perform field line tracing through GAMERA output data (trace.x)
=============================================================================

Introduction
------------

This executable takes output from an MHD simulation and extracts magnetic
field lines. Seed points are set uniformly in the in the equator (Z=0), within
a region selected by the user. See below for more details.

Example XML file

.. code-block:: xml

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
           <points Nx1="1" Nx2="72" Nx3="16" grType="SPHERICAL"/>
           <radius min="2.05" max="2.05"/>
           <phi min="180.0" max="360.0"/>
           <theta min="5.0" max="45.0"/>
       </Chimp>
   </KAIJU>

Parameter Descriptions
----------------------

``<sim>`` (optional): Specify identifying information for this computation.

``runid`` (optional, default ``"Sim"``): String specifying an identifier for
this run of ``psd.x``. A best practice is to use the ``runid`` in the name of
the XML file.

``<time>`` (optional): Specify time range and interval for calculation.

``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) relative
to start of simulation results used as input.

``dt`` (optional, default ``"1.0"``): Time interval and output cadence
(simulated seconds) for calculation.

``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds)
relative to start of simulation results used as input.

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
magnetic fields. Includes velocity vector, density and pressure in the output
file. See ``ebtypes.F90``.

``isMPI`` (optional, default ``"false"``): Set to ``"true"`` is the MAGE
results file was generated with an MPI version of the model. See
``eblCstd.F90``.

``<domain>`` (optional): Options for the problem domain

``dtype`` (optional, default ``"SPH"``): Domain over which to perform CHIMP
calculations, separate from grid, enables the user to perform calculation on
a subset of the grid to reduce computation where it is not needed - See
``gridloc.F90``. Valid values are ``"SPH"``, ``"LFM"``, ``"LFMCYL"``,
``"MAGE"``, ``"EGG"``, ``"ELL"``.

``rClosed`` (optional, default set by choice of ``units/uid``): Radial value
for field line endpoint to reach to be considered closed - See
``chmpunits.F90``.

``rmax`` (optional, default computed): Maximum radius of Domain region - See
``gridloc.F90``.

``rmin`` (optional, default computed): Minimum radius of Domain region - See
``gridloc.F90``.

``xSun`` (optional,default 20.0): if dType is "LFM" or "MAGE", the Domain
region includes all i-shells which have distances along the Earth-Sun line is
less than this value (in Re)

``xTail`` (optional,default -100.0): if dType is "LFM" or "MAGE", the Domain
region includes cells in the magnetotail up until this value (in Re)

``yzMax`` (optional,default 40.0): if dType is "LFM" or "MAGE", the Domain
region includes cells with Y and Z coordinates between +/- yzMax (in Re)

``<output>`` (optional): Options related to driver output

``timer`` (optional, default ``"false"``): Set to ``"true"`` to turn time
flags on See ``starter.F90``.

``tsOut`` (optional, default ``"10"``): Cadence to output diagnostics to
run-log file See ``starter.F90``.

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

``<units>`` (optional): Name of units system used in the model run.

``uID`` (optional, default ``"Earth"``): See chmpunits.F90. Valid values are
``"EARTH"``, ``"EARTHCODE"``, ``"JUPITER"``, ``"JUPITERCODE"``, ``"SATURN"``,
``"SATURNCODE"``, ``"HELIO"``, ``"LFM"``, ``"LFMJUPITER"``.

``<interp>`` (optional): Options related to interpolation

``wgt`` (optional, default ``"TSC"``): Sets 1D interpolation type. Valid
values are ``"TSC"`` (1D triangular shaped cloud), ``"LIN"`` (linear),
``"QUAD"`` (parabolic). See ``starter.F90``.

Parameters specific to field line tracing

``<tracer>`` (optional): Options related to field line tracing performed by
CHIMP

``epsds`` (optional, default ``"1.0e-2"``): Tolerance for field line tracing
computations See chmpdefs.F90.

``<points>`` (optional): Options related to initialization for field line seed
points. All points are set uniformly between min/max bounds. 2D grids are
seeded in the magnetic equator (Z=0).

``grType`` (optional, default ``"RP"``): Sets seed point grid type. Valid
values are ``"RP"`` (2D polar grid), ``"XY"`` (2D cartesian grid), ``"RTP"``
or ``"SPHERICAL"`` (3D spherical grid). See ``tracex.F90``.

``Nx1`` (optional, default ``"100"``): Number of cells in X or R depending on
grid specified.

``Nx2`` (optional, default ``"100"``): Number of cells in Y or Phi depending
on grid specified.

``Nx3`` (optional, default ``"1"``): Number of cells in Theta, set only if
grType is ``"RTP"`` or ``"SPHERICAL"``.

``x1Max`` (optional, default ``"10.0"``): Maximum value of the X dimension
used to initialize the seed points. Used if grType is ``"XY"``.

``x1Min`` (optional, default ``"-10.0"``): Minimum value of the X dimension
used to initialize the seed points.  Used if grType is ``"XY"``.

``x2Max`` (optional, default ``"10.0"``): Maximum value of the Y dimension
used to initialize the seed points. Used if grType is ``"XY"``.

``x2Min`` (optional, default ``"-10.0"``): Minimum value of the Y dimension
used to initialize the seed points. Used if grType is ``"XY"``.

``<radius>`` (optional): Radial range for initialization of seed point
locations. Used if grType is ``"RP"``, ``"RTP"`` or ``"SPHERICAL"``.

``min`` (optional, default ``"5.0"``): Minimum seed points are initialized at.

``max`` (optional, default ``"25.0"``): Maximum radius seed points are
initialized at.

``<phi>`` (optional): Azimuthal angle range for initialization of seed point
locations. Used if grType is ``"RP"``, ``"RTP"`` or ``"SPHERICAL"``.

``min`` (optional, default ``"0.0"``): Minimum longitude (in degrees) seed
points are initialized at. A value of 0 corresponds to the +X direction.

``max`` (optional, default ``"360.0"``\): Maximum longitude (in degrees) seed
points are initialized at.

``<theta>`` (optional): Polar angle range for initialization of seed point
locations. Used if grType is ``"RTP"`` or ``"SPHERICAL"``.

``min`` (optional, default ``"0.0"``): Minimum longitude (in degrees) seed
points are initialized at. A value of 0 corresponds to the +X direction.

``max`` (optional, default ``"90.0"``): Maximum longitude (in degrees) seed
points are initialized at.

Run Script example
------------------

An example pbs script on ``cheyenne``, RunTrace.pbs to submit field line
tracing job:

.. code-block:: bash

   #!/bin/bash
   #PBS -A P28100045
   #PBS -N traceFL
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=4:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72
   export EXE=${slice}
   export RUNID=${PBS_JOBNAME}

   #Replace this with your module set
   module purge
   module restore kaiju

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   echo "Running $EXE"
   ./trace.x ${RUNID}.xml > ${RUNID}.out
   date

This job can be submitted with the command

.. code-block:: bash

   qsub RunTrace.pbs

This will create an HDF5 file with a similar group structure for time slices
(Step#0,Step#1,xxx) but with each time slice also containing a group for each
field line seed, i.e.

Step#0

Line#0

Line#1

Step#1

Line#0

Line#1

In the above example, the field line seed points are set in a spherical
coordinate system with Nx1,Nx2,Nx3 points in r,phi,theta.

The script kaiju/scripts/genXLine.py can be run on the output HDF5 data to
create an XDMF file that can be read by VisIt or Paraview.
