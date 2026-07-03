``sctrack.x`` - Interpolate ``kaiju`` results along a spacecraft trajectory
===========================================================================


Introduction
------------

``sctrack.x`` is a compiled tool that interpolates ``kaiju`` results onto a
user-specified trajectory, typically representing a spacecraft.  This is used
when comparing ``kaiju`` results to measured data.


Before you start
----------------

In order to use ``sctrack.x``, you need the following:

* A set of output files from a ``kaiju`` magnetosphere run.

* An HDF5 file containing the spacecraft trajectory, with variables ``T``,
  ``X``, ``Y``, and ``Z``.

* An XML file describing the magnetosphere results and the settings for
  ``sctrack.x``.


``kaiju`` magnetosphere result files
------------------------------------

Any set of ``kaiju`` magnetosphere results can be used. This documentation
assumes you used the MPI version the the ``kaiju`` software, e.g.
``gamera_mpi.x``. This should generate a number of HDF5 files which contain
the model results from each MPI rank.

Assume all result files are in the current directory.


The XML file for ``sctrack.x``
------------------------------

The XML file read by ``chop.x`` should look something like this (file
``chop.xml``):

.. code-block:: xml

    <?xml version="1.0" ?>
    <Kaiju>
        <Chimp>
            <sim runid="CLUSTER1"/>
            <fields doMHD="T" grType="LFM" ebfile="geospace" isMPI="T"/>
            <parallel Ri="2" Rj="2" Rk="1"/>
            <units uid="EARTH"/>
            <trajectory H5Traj="CLUSTER1.sc.h5" doSmooth="F"/>
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

* ``<Chimp>`` (required): Inner element for CHIMP-specific elements.

* ``<sim>`` (optional): Specify identifying information for the ``kaiju``
  results to use,

  * ``runid`` (optional, default ``"Sim"``): String specifying the runid for
    the ``kaiju`` results to use as the basis for the extraction.

* ``<fields>`` (required): Describes the results from a MAGE model run.

  * ``doMHD`` (optional, default ``"F"``): Set to ``"T"`` to pass the full set
    of magnetohydrodynamic variables to ``sctrack.x``, rather than just the
    electric and magnetic fields. Includes velocity vector, density and
    pressure in the output file. See ``ebtypes.F90``.

  * ``grType`` (optional, default ``"EGG"``): String specifying grid type used
    by the MAGE output files. Valid values are ``"EGG"``, ``"LFM"``,
    ``"SPH"``. If the string is not one of the supported grid types, the
    default value (``"EGG"``) is used, and a warning message is printed. Note
    that for a magnetosphere simulation, the ``grType`` is nearly always
    ``"LFM"``.

  * ``ebfile`` (optional, default ``"ebdata"``): Root name for HDF5 files
    containing the results produced by a MAGE model run. This is usually the
    same as the runid.

  * ``isMPI`` (optional): If ``"T"``, the MAGE results are from an MPI run.

* ``<units>`` (optional): Data on units system used in the model run.

  * ``uid`` (optional, default ``"EARTH"``): Name of units system used in the
    model run. See ``chmpunits.F90``. Valid values are ``"EARTH"``,
    ``"EARTHCODE"``, ``"JUPITER"``, ``"JUPITERCODE"``, ``"SATURN"``,
    ``"SATURNCODE"``, ``"HELIO"``, ``"LFM"``, ``"LFMJUPITER"``..

* ``<trajectory>`` (required): Options specific to the spacecraft trajectory.

  ``H5Traj`` (required): Name of HDF5 file containing spacecraft trajectory.
  The interpolated data will be written back to this file.

  ``doSmooth`` (required): Set to ``"T"`` to smooth the trajectory, ``F``
  otherwise.


Running ``sctrack.x``
---------------------

You can now run ``sctrack.x`` with the following command:

.. code-block:: bash

  sctrack.x sctrack.xml

The extracted data will be in the file defined in the ``H5Traj`` attribute.
