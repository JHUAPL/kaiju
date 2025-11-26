``slice.x`` - Extract 2D subsets of ``kaiju`` results
=====================================================


Introduction
------------

``slice.x`` extracts a 2D (``Z=z0`` or ``Y=y0``) slice from GAMERA output on
either the MAGE grid or interpolated onto a cartesian or a polar grid. Output
is saved to an HDF5 file that can be visualized and analyzed similarly to
GAMERA output files.


Before you start
----------------

In order to use ``slice.x``, you need the following:

* A set of output files from a ``kaiju`` magnetosphere run.

* An XML file describing the magnetosphere results and the settings for
  ``slice.x``.


``kaiju`` magnetosphere result files
------------------------------------

Any set of ``kaiju`` magnetosphere results can be used. This documentation
assumes you used the MPI version the the ``kaiju`` software, e.g.
``gamera_mpi.x``. This should generate a number of HDF5 files which contain
the model results from each MPI rank.

Assume all result files are in the current directory.


The XML file for ``slice.x``
----------------------------

The XML file read by ``slice.x`` should look something like this (file
``slice.xml``):

.. code-block:: xml

    <?xml version="1.0" ?>
    <Kaiju>
        <Chimp>
            <sim runid="ebXY"/>
            <time T0="0.0" dt="15.0" tFin="3600.0"/>
            <fields doMHD="T" ebfile="geospace" grType="LFM" isMPI="T"/>
            <parallel Ri="2" Rj="2" Rk="1"/>
            <domain dtype="LFM" xSun="30.0" xTail="-50.0" yzMax="30.0"/>
            <units uid="EARTH"/>
            <slice doXY="T" grType="XY" z0="1.0" xSun="20" xTail="-20" yM="20" Nx1="400" Nx2="400"/>
            <interp wgt="TSC"/>
            <tracer epsds="0.20"/>
            <output doSlim="T" doTrc="T"/>
            <plasmapause doPP="F"/>
        </Chimp>
    </Kaiju>

The elements and attributes of the XML file are described below.

.. note::

    Standard XML is *case-sensitive*. However, the XML parser used by the
    Fortran code in ``kaiju`` is *case-insensitive*. Therefore, all of the
    XML tags and attributes described below can be written with any combintion
    of upper and lower case.

* ``<?xml version="1.0"?>`` (required): Specifies XML version.

* ``<Kaiju>`` (required): Outer element for all ``kaiju`` elements.

* ``<Chimp>`` (required): Inner element for Chimp-specific elements.

* ``<sim>`` (optional): Specify identifying information for this computation.

  * ``runid`` (optional, default ``"Sim"``): String specifying an identifier
    for this run of ``slice.x``. A best practice is to use the ``runid`` in
    the name of the XML file.

* ``<time>`` (optional): Specify time range and interval for magnetic field
  calculation.

  * ``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) for
    ground magnetic field calculation, relative to start of simulation results
    used as input.

  * ``dt`` (optional, default ``"1.0"``): Time interval and output cadence
    (simulated seconds) for ground magnetic field calculation.

  * ``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds) for
    ground magnetic field calculation, relative to start of simulation results
    used as input.

* ``<fields>`` (required): Describes the input data from a MAGE model run.

  * ``ebfile`` (optional, default ``"ebdata.h5"``): Path to HDF5 file
    containing the electric and magnetic fields computed by a MAGE model run.

  * ``grType`` (optional, default ``"EGG"``): String specifying grid type used
    by the MAGE output file. Valid values are ``"EGG"``, ``"LFM"``, ``"SPH"``.
    If the string is not one of the supported grid types, the default value
    (``"EGG"``) is used, and a warning message is printed.

  * ``doEBFix`` (optional, default ``"false"``): Set to ``"true"`` to "clean"
    the electric field E so that the dot product of the electric and magnetic
    fields is 0. See ``ebinterp.F90``.

  * ``doMHD`` (optional, default ``"false"``): Set to ``"true"`` to pass the
    full set of magnetohydrodynamic variables to CHIMP, rather than just
    electric and magnetic fields. See ``ebtypes.F90``.

  * ``isMPI`` (optional, default ``"false"``): Set to ``"true"`` is the MAGE
    results file was generated with an MPI version of the model. See
    ``eblCstd.F90``.

* ``<domain>`` (optional): Options for the problem domain

  * ``dtype`` (optional, default ``"SPH"``): Domain over which to perform
    CHIMP calculations, separate from grid, enables the user to perform
    calculation on a subset of the grid to reduce computation where it is not
    needed - See ``gridloc.F90``. Valid values are ``"SPH"``, ``"LFM"``,
    ``"LFMCYL"``, ``"MAGE"``, ``"EGG"``, and ``"ELL"``.

  * ``rClosed`` (optional, default set by choice of ``units/uid``): Radial
    value for field line endpoint to reach to be considered closed - See
    ``chmpunits.F90``.

  * ``rmax`` (optional, default computed): Maximum radius of Domain region.
    See ``gridloc.F90``.

  * ``rmin`` (optional, default computed): Minimum radius of Domain region.
    See ``gridloc.F90``.

  * ``xSun`` (optional,default 20.0): if dType is "LFM" or "MAGE", the Domain
    region includes all i-shells which have distances along the Earth-Sun line
    is less than this value (in Re)

  * ``xTail`` (optional,default -100.0): if dType is "LFM" or "MAGE", the
    domain region includes cells in the magnetotail up until this value (in
    Re).

  * ``yzMax`` (optional,default 40.0): if dType is "LFM" or "MAGE", the Domain
    region includes cells with Y and Z coordinates between +/- yzMax (in Re)

* ``<output>`` (optional): Options related to driver output

  * ``timer`` (optional, default ``"false"``): Set to ``"true"`` to turn time
    flags on See ``starter.F90``.

  * ``tsOut`` (optional, default ``"10"``): Cadence to output diagnostics to
    run-log file See ``starter.F90``.

  * ``doEQProj`` (optional, default ``"false"``): Set to ``true`` to include
    equatorial variables, projected down to magnetic equator along field line
    from cell location (i.e. Xeq,Yeq, if field line is open or closed etc) See
    ``chmpdefs.F90``.

  * ``doSlim`` (optional, default ``"false"``):  Set to ``"true"`` to remove
    vector electric field and current data from slice.x output See
    ``chmpdefs.F90``.

  * ``doTrc`` (optional, default ``"false"``): Similar to doEQProj, used iZzn
    slice.x See ``chmpdefs.F90``.

* ``<parallel>`` (optional): Options if ebfile was generated using an MPI
  version of the code (read if fields/doMPI is set to ``"true"``, file name in
  form of ebfile_Ri_Rj_Rk_i_j_k.gam.h5)

  * ``Ri`` (optional, default ``"1"``): Number of ranks used in decomposition
    of ``"i"`` dimension See ``iotable.F90``.

  * ``Rj`` (optional, default ``"1"``): Number of ranks used in decomposition
    of ``"j"``    zdimension See ``iotable.F90``.

  * ``Rk`` (optional, default ``"1"``): Number of ranks used in decomposition
    of ``"k"`` dimension See ``iotable.F90``.

  * ``doOldNaming`` (optional, default ``"false"``): Allow for backward
    compatibility for MHD files generated with the now deprecated naming
    convention See ``chmpdefs.F90``.

*  ``<tracer>`` (optional): Options related to field line tracing performed by
   CHIMP

  * ``epsds`` (optional, default ``"1.0e-2"``): Tolerance for field line
    tracing computations See ``chmpdefs.F90``.

* ``<units>`` (optional): Name of units system used in the model run.

  * ``uID`` (optional, default ``"EARTH"``): See ``chmpunits.F90``. Valid
    values are ``"EARTH"``, ``"EARTHCODE"``, ``"JUPITER"``, ``"JUPITERCODE"``,
    ``"SATURN"``, ``"SATURNCODE"``, ``"HELIO"``, ``"LFM"``, and
    ``"LFMJUPITER"``.

* ``<interp>`` (optional): Options related to interpolation

  * ``wgt`` (optional, default ``"TSC"``): Sets 1D interpolation type. Valid
    values are ``"TSC"`` (1D triangular shaped cloud), ``"LIN"`` (linear),
    ``"QUAD"`` (parabolic). See ``starter.F90``.

* ``<slice>`` (optional): Options specific to slice.x driver, see
  ``sliceio.F90``.

  * ``doXY`` (optional, default ``"true"``): Perform a slice in the XY plane,
    if ``"false"``, perform a slice in the XZ plane.

  * ``z0`` (optional, default ``"0.0"``): Take a slice at Z=z0, if doXY
    ``"true"``.

  * ``y0`` (optional, default ``"0.0"``): Take a slice at Y=y0, if doXY
    ``"false"``.

  * ``grType`` (optional, default ``"XY"``): String specifying an identifier
    for the grid to be used in sliced plane. Valid Values are ``"XY"``
    (Cartesian), ``"RP"`` (spherical), and ``"LFM2D"`` (MAGE grid).

  * ``Nx1`` (optional, default ``"128"``): Number of cells in X or R depending
    on grid specified. Not used if grType is set to ``"LFM2D"``.

  * ``Nx2`` (optional, default ``"256"``): Number of cells in Y/Z or R
    depending on grid specified. Not used if grType is set to ``"LFM2D"``.

  * ``Npow`` (optional, default ``"0"``): Number of times to increase the grid
    used in the sliced plane by a factor of 2.

  * ``xSun`` (optional, default ``"12.5"``): If ``grType="XY"``, sets the
    maximum value of X used to initialize the slice grid. If
    ``grType="LFM2D"``, used similarly to ``domain/xSun`` where all i-shells
    which have distances along the Earth-Sun line is less than this value (in
    Re) are included in the sliced grid.

  * ``xTail`` (optional, default ``"-20.0"``): If ``grType="XY"``, sets the
    minimum value of X used to initialize the slice grid.

  * ``yM`` (optional, default ``"20.0"``): If ``grType="XY"``, the
    minimum/maximum values of the grid of the second axis (either Y or Z
    depending on the value of doXY) are set to +/-yM.

  * ``Rin`` (optional, default ``"Rin of ebFile grid"``): If ``grType="RP"``,
    sets the minimum value of R used to initialize the slice grid.

  * ``Rout`` (optional, default ``"Max R in the +X direction of ebFile
    grid"``): If ``grType="RP"``\ , sets the maximum value of R used to
    initialize the slice grid.

  * ``Pin`` (optional, default ``"0.0"``): If ``grType="RP"``, sets the
    minimum value of longitude (in degrees) used to set the slice grid. A
    value of 0 corresponds to the +X direction.

  * ``Pout`` (optional, default ``"360.0"``): If ``grType="RP"``, sets the
    maximum value of longitude (in degrees) used to set the slice grid.

  * ``doLog`` (optional, default ``"false"``): Set to ``"true"`` to distribute
    cells uniformly in log-space between Rin/Rout values.

  * ``dSlc`` (optional, default ``"0.05"``): Spacing used to average slice over
    y0+/-dSlc or z0+/-dSlc, depending on the value of doXY.

* ``<plasmapause>`` (optional): Options for calculation to determine
  plasmapause location in MHD file

  * ``doPP`` (optional, default ``"false"``): Set to ``.true.`` to calculate
    plasmapause location in the equator and include it in output file for
    ``slice.x``. See ``chmpfields.F90``.

  * ``Lin`` (optional, default ``"2.0"``): Minimum L-shell value to begin
    plasmapause calculation. See ``plasmaputils.F90``.

  * ``Lout`` (optional, default ``"6.0"``):  Maximum L-shell value to end
    plasmapause calculation. See ``plasmaputils.F90``.

  * ``Nl`` (optional, default ``"80"``): Number of cells/steps in L-shell.
    See ``plasmaputils.F90``.

  * ``Nphi`` (optional, default ``"72"``): Number of cells/steps in longitude.
    See ``plasmaputils.F90``.

  * ``phi0`` (optional, default ``"0.0"``): Minimum longitude to perform
    plasmapause calculation between (in radians, 0 is in the +X direction) -
    See ``plasmaputils.F90``.

  * ``phi1`` (optional, default ``"2*PI"``): Maximum longitude to perform
    plasmapause calculation between See ``plasmaputils.F90``.


Running ``slice.x``
-------------------

You can now run ``slice.x`` with the following command:

.. code-block:: bash

  slice.x slice.xml

.. note:: A ``slice.x`` run can be quite long...

The extracted data will be in a new file called ``ebXY.eb.h5``. The file
can be quite large, depending on the parameters chosen. The XML shown above
resulted in a file size of 4 GB.
