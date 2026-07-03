Calculate phase space density from test particle simulations (psd.x)
====================================================================

Introduction
------------

This executable takes output from an MHD and test particle simulation and
calculates the time evolution of the resulting phase space density. The test
particles are weighted assuming an initial phase space density distribution,
taken from either the MHD solution at the particles location or from an
external hdf5 file provided by the user. The weight corresponds to the number
of real particles each test particle represents. For a detailed description of
how the calculations are performed, see appendix A2 of `Sorathia et al (2018)
<https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JA025506>`_.

Example XML file
----------------

.. code-block:: xml

   <?xml version="1.0" ?>
   <Kaiju>
     <Chimp>
       <sim runid="eRBpsdH5" doShape="T"/>
       <tps species="e"/>
       <time T0="36600.0" dt="60.0" tFin="36610.0"/>
       <fields doMHD="T" ebfile="msphere" grType="LFM" isMPI="T"/>
       <parallel Ri="8" Rj="8" Rk="1"/>
       <units uid="EARTH"/>
       <output dtOut="60.0" tsOut="100" doFat="T"/>
       <radius N="30" doLog="F" max="10" min="2.5"/>
       <phi N="24" max="360.0" min="0.0"/>
       <alpha N="9" max="90.0" min="0.0"/>
       <energy N="30" doLog="T" max="7000.0" min="50.0"/>
       <domain rmax="20.0" rmin="1.05"/>
       <tracer epsds="0.05"/>
       <population f0="HDF5IN" f0data="psdInit.h5" ns="1" ne="20" popid="eRB"/>
     </Chimp>
   </Kaiju>

Parameter Descriptions
----------------------

``<sim>`` (optional): Specify identifying information for this computation.

``runid`` (optional, default ``"Sim"``): String specifying an identifier for
this run of ``psd.x``. A best practice is to use the ``runid`` in the name of
the XML file.

``doShape`` (optional, default ``"true"``): Set to ``"true"`` to use a
shaping function when calculating phase space density.

``<tps>`` (optional): Options related to test particles

``species`` (optional, default ``"X"``): Species simulated by push.x that
you would like to weight. For full list of available species, see the
getSpecies() function in chmpunits.F90.

``<time>`` (optional): Specify time range and interval for magnetic field
calculation.

``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) for
ground magnetic field calculation, relative to start of simulation results
used as input.

``dt`` (optional, default ``"1.0"``): Time interval and output cadence
(simulated seconds) for ground magnetic field calculation.

``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds) for
ground magnetic field calculation, relative to start of simulation results
used as input.

``<fields>`` (required): Describes the input data from a MAGE model run.

``ebfile`` (optional, default ``"ebdata.h5"``): Path to HDF5 file containing
the electric and magnetic fields computed by a MAGE model run.

``grType`` (optional, default ``"EGG"``): String specifying grid type used
by the MAGE output file. Valid values are ``"EGG"``, ``"LFM"``, ``"SPH"``. If
the string is not one of the supported grid types, the default value
(``"EGG"``) is used, and a warning message is printed.

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

``rho0`` (optional, default ``"1.0"``): Default density used if not using MHD
values to determine distribution function used to weight particles.

``kT0`` (optional, default ``"1.0"``): Default temperature used if not using
MHD values to determine distribution function used to weight particles.

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

``doFat`` (optional, default ``"false"``): Set to ``"true"`` to include 4D
variable information in output files. See psdio.F90.

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

``<tracer>`` (optional): Options related to field line tracing performed by
CHIMP

``epsds`` (optional, default ``"1.0e-2"``): Tolerance for field line tracing
computations See chmpdefs.F90.

``<units>`` (optional): Name of units system used in the model run.

``uID`` (optional, default ``"Earth"``): See chmpunits.F90. Valid values are
``"EARTH"``, ``"EARTHCODE"``, ``"JUPITER"``, ``"JUPITERCODE"``, ``"SATURN"``,
``"SATURNCODE"``, ``"HELIO"``, ``"LFM"``, ``"LFMJUPITER"``.

``<energy>`` (optional): Options for initialization of phase space density
grid

``min`` (optional, default ``"1.0"``): Minimum energy (in keV) of grid.

``max`` (optional, default ``"100.0"``): Maximum energy (in keV) of grid.

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute cells uniformly
in log-space between min/max values.

``N`` (optional, default ``"15"``): Number of cells to use in this dimension.

``<alpha>`` (optional): Options for initialization of phase space density
grid

``min`` (optional, default ``"0.0"``): Minimum pitch angle (in degrees) of
grid.

``max`` (optional, default ``"360.0"``): Maximum pitch angle (in degrees) of
grid.

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute cells uniformly
in log-space between min/max values.

``N`` (optional, default ``"10"``): Number of cells to use in this dimension.

``<radius>`` (optional): Options for initialization of phase space density
grid

``min`` (optional, default ``"5.0"``): Minimum radius of grid.

``max`` (optional, default ``"25.0"``): Maximum radius of grid.

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute cells uniformly
in log-space between min/max values.

``N`` (optional, default ``"20"``): Number of cells to use in this dimension.

``<phi>`` (optional): Options for initialization of phase space density grid

``min`` (optional, default ``"0.0"``): Minimum longitude (in degrees) of grid.
A value of 0 corresponds to the +X direction.

``max`` (optional, default ``"360.0"``): Maximum longitude (in degrees) of
grid.

``doLog`` (optional, default ``"false"``): Default behavior is uniform
distribution between min/max. Set to ``"true"`` to distribute cells uniformly
in log-space between min/max values.

``N`` (optional, default ``"8"``): Number of cells to use in this dimension.

``<population>`` (optional): Options specific to phase space density
calculation, see psdinit.F90 and pdfuns.F90

``popid`` (optional, default ``"chimp"``): String specifying an identifier for
runID of test particle files being read in for PSD calculation.

``ns`` (optional, default ``"1"``): Starting job id/number of test particle
files that should be included in PSD calculation.

``ne`` (optional, default ``"4"``): Ending job id/number of test particle
files that should be included in PSD calculation.

``kTScl`` (optional, default ``"1.0"``): Factor used to scale temperature
used in calculation.

``f0`` (optional, default ``"Max"``): String specifying an identifier for
type of distribution function to use to weight the test particles. Valid
values are ``"MAX"`` (maxwellian based on local MHD parameters),
``"KAP"`` (kappa based on local MHD parameters), ``"RBSP"`` (function fit to
Van Allen Probes data),  ``"HDF5IN"`` (provided by an external hdf5 file).

``k0`` (optional, default ``"3.0"``): kappa value to use if f0=``"KAP"``

``f0data`` (optional, default ``"psd.h5"``: String specifying an name of hdf5
file containing PSD distribution to weight test particles with. Read if
f0=``"HDF5IN"``

``<stream>`` (optional): Options related to time dependent initialization of
test particles

``dShell`` (optional, default ``"1.0"``): Width of streaming region, used to
rescale PSD initial condition if test particles were injected over time.

Run Script Example
------------------

An example pbs script on ``cheyenne``, RunPSD.pbs to submit a phase space
density run:

.. code-block:: bash

   #!/bin/bash
   #PBS -A P28100045
   #PBS -N eRBpsdH5
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72

   export EXE=${CHIMPEXE:-"psd.x"}
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
   ${EXE} ${RUNID}.xml > ${RUNID}.out
   date

This job can be submitted with the command

.. code-block:: bash

   qsub RunPSD.pbs

Description of Output
---------------------

The run outputs three files:


``runID.wgt.h5``  - contains all the weights of the test particles, in order
of their particle ID

``runID.ps.h5``   - phase space density and intensity as a function of energy,
radius, and longitude, integrated over equatorial pitch angle.

``fPSD`` - phase space density [(keV s)^-3]

``jPSD`` - intensity [cm^-2 sr^-1 s^-1 keV^-1]

``Ntp``  - number of test particles in each cell

``dG``   - phase space density cell volume element

``runID.pseq.h5`` - outputs the moments of the phase space density
distribution function in the equatorial plane.

A fourth file is output if the ``doFat`` flag is set to ``"true"`` in the xml
file.

``runID.ps4.h5`` - 4D phase space density (fPSD) as a function of energy,
equatorial pitch angle, radius, and longitude.

``fPSD`` - phase space density [(keV s)^-3]

``Ntp``  - number of test particles in each cell

``dG``   - phase space density cell volume element

``dGp``  - phase space density cell momentum volume element

``dGx``  - phase space density cell spatial volume element
