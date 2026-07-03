``calcdb.x`` - Computing ground delta-B values from ``kaiju`` results
=====================================================================


Introduction
------------

The tool ``calcdb.x`` is a compiled Fortran program which uses the output from
a MAGE simulation to compute the values of the ground
magnetic field disturbance.


Before you start
----------------

In order to use ``calcdb.x``, you need the following:

* A set of output files from a ``kaiju`` magnetosphere run.

* An XML file describing the magnetosphere results and the settings for the
  ground delta B calculation.

* (optional) A PBS job file to run ``calcdb.x``.


``kaiju`` magnetosphere result files
------------------------------------

Any set of ``kaiju`` magnetosphere results can be used. This documentation
assumes you used the MPI version the the ``kaiju`` software, e.g.
``gamera_mpi.x``. This should generate a number of HDF5 files which contain
the model results from each MPI rank.

Assume all result files are in the current directory.


The XML file for ``calcdb.x``
-----------------------------

The XML file read by ``calcdb.x`` should look something like this (file
``calcdb.xml``):

.. code-block:: xml

    <?xml version="1.0"?>
    <Kaiju>
        <Chimp>
            <sim runid="geospace"/>
            <time T0="0.0" dt="60.0" tFin="14400.0"/>
            <fields ebfile="geospace" grType="LFM" doJ="T" isMPI="true"/>
            <parallel Ri="4" Rj="4" Rk="1"/>
            <grid Nlat="360" Nlon="720" Nz="1"/>
            <parintime NumB="4"/>
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
    the ``kaiju`` results to use in the ground delta B calculation.

* ``<time>`` (optional): Specify time range and interval for magnetic field
  calculation.

  * ``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) for
    ground magnetic field calculation, relative to start of simulation
    results used as input.

  * ``dt`` (optional, default ``"1.0"``): Time interval and output cadence
    (simulated seconds) for ground magnetic field calculation.

  * ``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds) for
    ground magnetic field calculation, relative to start of simulation
    results used as input.

* ``<fields>`` (required): Describes the MAGE model results to use.

  * ``ebfile`` (optional, default ``"ebdata"``): Root name for HDF5 files
    containing the results produced by a MAGE model run. This is usually the
    same as the runid.

  * ``grType`` (optional, default ``"EGG"``): String specifying grid type used
    by the MAGE output files. Valid values are ``"EGG"``, ``"LFM"``, and
    ``"SPH"``. If the string is not one of the supported grid types, the
    default value (``"EGG"``) is used, and a warning message is printed. Note
    that for a magnetosphere simulation, the ``grType`` is always ``"LFM"``.

  * ``doJ`` (required, must be ``"T"``): If ``"T"``, compute currents from the
    ``kaiju`` model results.

  * ``isMPI`` (optional): If ``"true"``, the MAGE results are from an MPI run.

* ``<parallel>`` (required): Describes the MPI decomposition of the MAGE model
  run.

  * ``Ri`` (optional, default ``"1"``): Number of ranks used in MPI
    decomposition  of ``i`` dimension.

  * ``Rj`` (optional, default ``"1"``): Number of ranks used in MPI
    decomposition  of ``j`` dimension.

  * ``Rk`` (optional, default ``"1"``): Number of ranks used in MPI
    decomposition of ``k`` dimension.

* ``<grid>`` (optional): Options to specify the grid on the ground used in
  ``calcdb.x``.

  * ``Nlat`` (optional, default ``"45"``): Number of latitude cells.

  * ``Nlon`` (optional, default ``"90"``): Number of longitude cells.

  * ``Nz`` (optional, default ``"2"``): Number of altitude cells.


Running ``calcdb.x``
--------------------

We can now run ``calcdb.x`` with the following command:

.. code-block:: bash

  calcdb.x calcdb.xml

For large simulations, this process can take a long time. A more efficient
approach is to use PBS to run multiple ``calcdb.x`` jobs in parallel, and
combine the results. To do this, add the following element to the XML file. in
the ``<Chimp>`` element:

.. code-block:: xml

  <parintime NumB="4"/>

where:

* ``<parintime>`` (optional): Options to run a PBS job array of ``calcdb.x``
  to increase calculation speed.

  * ``NumB`` (optional, default ``"0"``): Number of jobs into which the
    calculation will be split for parallel computation.

We also need a PBS script, as shown below (the example runs on ``derecho``):

.. code-block:: bash

  #PBS -a YOUR_PBS_ACCOUNT
  #PBS -q main
  #PBS -l walltime=01:00:00
  #PBS -l select=1:ncpus=128
  #PBS -l job_priority=economy

  module purge
  # Load your modules here.
  module list

  export KAIJU_INSTALL_DIR=/path/to/your/kaiju/code
  source $KAIJU_INSTALL_DIR/scripts/setupEnvironment.sh
  export KAIJU_BUILD_DIR=/path/to/your/kaiju/build
  cp $KAIJU_BUILD_DIR/bin/calcdb.x ./calcdb.x

  export JNUM=${PBS_ARRAY_INDEX:-0}
  ./calcdb.x calcdb.xml $JNUM

Now submit the job to the PBS system:

.. code-block:: bash

  $ qsub -J 1-4 calcdb.pbs

.. note::

  The second digit in ``-J 1-4`` must match the value of ``NumB`` in the XML
  file.

When the calculation has finished, your directory will contain a set of files
named like this:

.. code-block:: bash

  ls -1 *.deltab.h5
  geospace.0001.deltab.h5
  geospace.0002.deltab.h5
  geospace.0003.deltab.h5
  geospace.0004.deltab.h5

Now combine the results into a single file with the ``pitmerge.py`` command
(``pitmerge.py`` is provided in the ``kaipy`` python package):

.. code-block:: bash

  pitmerge.py -runid geospace

where ``runid`` is the runid for the ``kaiju`` run being used for the
calculation. This command will create the file ``geospace.deltab.h5``, which
contains the ground delta B values for the simulation.
