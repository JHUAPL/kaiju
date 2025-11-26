``chop.x`` - Extract a subset of ``kaiju`` results
==================================================


Introduction
------------

``chop.x`` is a compiled tool that extracts a 3D portion of the domain from
``kaiju`` output on either the MAGE grid or interpolated onto a Cartesian or a
spherical grid. ``chop.x`` also performs additional calculations such as field
line tracing. Output is to an HDF5 file that can be visualized and analyzed
similarly to ``kaiju`` output files.


Before you start
----------------

In order to use ``chop.x``, you need the following:

* A set of output files from a ``kaiju`` magnetosphere run.

* An XML file describing the magnetosphere results and the settings for
  ``chop.x``.


``kaiju`` magnetosphere result files
------------------------------------

Any set of ``kaiju`` magnetosphere results can be used. This documentation
assumes you used the MPI version the the ``kaiju`` software, e.g.
``gamera_mpi.x``. This should generate a number of HDF5 files which contain
the model results from each MPI rank.

Assume all result files are in the current directory.


The XML file for ``chop.x``
---------------------------

The XML file read by ``chop.x`` should look something like this (file
``chop.xml``):

.. code-block:: xml

    <?xml version="1.0"?>
    <Kaiju>
        <Chimp>
            <sim runid="geospace"/>
            <time T0="0.0" dt="30.0" tFin="3600.0"/>
            <fields doMHD="T" ebfile="geospace" grType="LFM" isMPI="T"/>
            <parallel Ri="2" Rj="2" Rk="1"/>
            <domain dtype="LFMCYL" xSun="30.0" xTail="-50.0" yzMax="30.0"/>
            <units uid="EARTH"/>
            <chop grType="RTP" x1Max="15.0" x1Min="2.25" x2Max="180" x2Min="0"
             x3Max="360" x3Min="0" Nx1="1" Nx2="360" Nx3="720"/>
            <interp wgt="TSC"/>
            <tracer epsds="0.20"/>
            <output doTrc="T"/>
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

* ``<sim>`` (optional): Specify identifying information for the ``kaiju``
  results to use,

  * ``runid`` (optional, default ``"Sim"``): String specifying the runid for
    the ``kaiju`` results to use as the basis for the extraction.

* ``<time>`` (optional): Specify time range and interval for extraction.

  * ``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) for
    extraction, relative to start of simulation results used as input.

  * ``dt`` (optional, default ``"1.0"``): Time interval and output cadence
    (simulated seconds) for extraction.

  * ``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds) for
    extraction, relative to start of simulation results used as input.

* ``<fields>`` (required): Describes the results from a MAGE model run.

  * ``doMHD`` (optional, default ``"F"``): Set to ``"T"`` to pass the full set
    of magnetohydrodynamic variables to ``chop.x``, rather than just the
    electric and magnetic fields. Includes velocity vector, density and
    pressure in the output file. See ``ebtypes.F90``.

  * ``ebfile`` (optional, default ``"ebdata"``): Root name for HDF5 files
    containing the results produced by a MAGE model run. This is usually the
    same as the runid.

  * ``grType`` (optional, default ``"EGG"``): String specifying grid type used
    by the MAGE output files. Valid values are ``"EGG"``, ``"LFM"``,
    and ``"SPH"``. If the string is not one of the supported grid types, the
    default value (``"EGG"``) is used, and a warning message is printed. Note
    that for a magnetosphere simulation, the ``grType`` is nearly always
    ``"LFM"``.

  * ``isMPI`` (optional): If ``"T"``, the MAGE results are from an MPI run.

* ``<parallel>`` (required): Describes the MPI decomposition of MAGE model
  run.

  * ``Ri`` (optional, default ``"1"``): Number of ranks used in MPI
    decomposition  of ``i`` dimension.

  * ``Rj`` (optional, default ``"1"``): Number of ranks used in MPI
    decomposition of ``j`` dimension.

  * ``Rk`` (optional, default ``"1"``): Number of ranks used in MPI
    decomposition of ``k`` dimension.

* ``<domain>`` (optional): Options for the problem domain

  * ``dtype`` (optional, default ``"SPH"``): Domain over which to perform
    CHIMP calculations, separate from grid. This enables the user to perform
    calculation on a subset of the grid to reduce computation where it is not
    needed - See ``gridloc.F90``. Valid values are ``"SPH"``, ``"LFM"``,
    ``"LFMCYL"``, ``"MAGE"``, ``"EGG"``, ``"ELL"``.

  * ``xSun`` (optional,default ``"20.0"``: If ``dType`` is ``"LFM"`` or
    ``"MAGE"``, the domain region includes all i-shells which have distances
    along the Earth-Sun line less than this value (in Re).

  * ``xTail`` (optional,default ``"-100.0"``): if ``dType`` is ``"LFM"`` or ``"MAGE"``,
    the domain region includes cells in the magnetotail up until this value
    (in Re).

  * ``yzMax`` (optional,default ``"40.0"``): if ``dType`` is ``"LFM"`` or ``"MAGE"``, the
    domain region includes cells with Y and Z coordinates between +/- ``yzMax``
    (in Re).

* ``<units>`` (optional): Data on units system used in the model run.

  * ``uid`` (optional, default ``"Earth"``): Name of units system used in the
    model run. See ``chmpunits.F90``. Valid values are ``"EARTH"``,
    ``"EARTHCODE"``, ``"JUPITER"``, ``"JUPITERCODE"``, ``"SATURN"``,
    ``"SATURNCODE"``, ``"HELIO"``, ``"LFM"``, ``"LFMJUPITER"``..

* ``<chop>`` (optional): Options specific to ``chop.x``.

  * ``grType`` (optional, default ``"XYZ"``): String specifying an identifier
    for the grid type for data extraction. Valid Values are ``"XYZ"``
    (Cartesian), ``"RTP"`` (r-theta-phi, i.e. spherical), ``"LFM"`` (MAGE grid).

  * ``x1Max`` (optional, default ``"10.0"``): Maximum value of the first
    dimension used to initialize the extracted grid. If ``grTyp="LFM"``,
    ``x1Max`` is used similar to ``domain/xSun``.

  * ``x1Min`` (optional, default ``"-x1Max"``): Minimum value of the first
    dimension used to initialize the extracted grid.

  * ``x2Max`` (optional, default ``"10.0"``): Maximum value of the second
    dimension used to initialize the extracted grid.

  * ``x2Min`` (optional, default ``"-x2Max"``): Minimum value of the second
    dimension used to initialize the extracted grid.

  * ``x3Max`` (optional, default ``"10.0"``): Maximum value of the third
    dimension used to initialize the extracted grid.

  * ``x3Min`` (optional, default ``"-x3Max"``): Minimum value of the second
    dimension used to initialize the extracted grid.

  * ``Nx1`` (optional, default ``"64"``): Number of cells in X or R depending
    on grid specified.

  * ``Nx2`` (optional, default ``"64"``): Number of cells in Y or Theta
    depending on grid specified.

  * ``Nx3`` (optional, default ``"64"``): Number of cells in Z or Phi
    depending on grid specified.

  * ``rClosed`` (optional, default set by choice of ``units/uid``): Radial
    value for field line endpoint to reach to be considered closed. See
    ``chmpunits.F90``.

  * ``rmax`` (optional, default computed): Maximum radius of domain region.
    See ``gridloc.F90``.

  * ``rmin`` (optional, default computed): Minimum radius of domain region.
    See ``gridloc.F90``.

* ``<output>`` (optional): Options related to ``chop.x```` output

  * ``timer`` (optional, default ``"false"``): Set to ``"true"`` to turn time
    flags on See ``starter.F90``.

  * ``tsOut`` (optional, default ``"10"``): Cadence to output diagnostics to
    run-log file. See ``starter.F90``.

  * ``doEQProj`` (optional, default ``"false"``): Set to ``"true"`` to include
    equatorial variables, projected down to magnetic equator along field line
    from cell location (i.e. Xeq, Yeq, if field line is open or closed etc).
    See ``chmpdefs.F90``.

  * ``doSlim`` (optional, default ``"false"``):  Set to ``"true"`` to remove
    vector electric field and current data from ``chop.x`` output. See
    ``chmpdefs.F90``.

  * ``doTrc`` (optional, default ``"false"``): Similar to ``doEQProj``, used in
    ``chop.x``. See ``chmpdefs.F90``.

* ``<parallel>`` (optional): Use if ``ebfile`` was generated using an MPI
  version of the code.

  * ``Ri`` (optional, default ``"1"``): Number of ranks used in decomposition
    of the first dimension See ``iotable.F90``.

  * ``Rj`` (optional, default ``"1"``): Number of ranks used in decomposition
    of the second dimension See ``iotable.F90``.

  * ``Rk`` (optional, default ``"1"``): Number of ranks used in decomposition
    of the third dimension See ``iotable.F90``.

  * ``doOldNaming`` (optional, default ``"false"``): Allow for backward
    compatibility for MHD files generated with the now-deprecated naming
    convention. See ``chmpdefs.F90``.

* ``<tracer>`` (optional): Options related to field line tracing.

  * ``epsds`` (optional, default ``"1.0e-2"``): Tolerance for field line
    tracing computations. See ``chmpdefs.F90``.

* ``<interp>`` (optional): Options related to interpolation.

  * ``wgt`` (optional, default ``"TSC"``): Sets 1D interpolation type. Valid
    values are ``"TSC"`` (1D triangular shaped cloud), ``"LIN"`` (linear),
    ``"QUAD"`` (parabolic). See ``starter.F90``.


Running ``chop.x``
------------------

You can now run ``chop.x`` with the following command:

.. code-block:: bash

  chop.x chop.xml

.. note:: A ``chop.x`` run can be quite long...

The extracted data will be in a new file called ``geospace.eb3.h5``. The file
can be quite large, depending on the parameters chosen. The XML shown above
resulted in a file size of 3 GB.
