``create_gamhelio_ensemble.py``
===============================

The Python script ``create_gamhelio_ensemble.py`` was developed to simplify
the process of configuring and running an ensemble of heliosphere models with
the ``kaiju`` software. It provides an interactive, prompt-driven interface to
specify all of the parameters needed for an ensemble of model runs.


Running the ``create_gamhelio_ensemble.py`` script
--------------------------------------------------

The ``create_gamhelio_ensemble.py`` script is provided as part of the
``kaiju`` software. It is found at
``$KAIJUHOME/scripts/makeitso-gamhelio/create_gamhelio_ensemble.py``, where
``$KAIJUHOME`` is the location of your ``kaiju`` software tree. After
configuring your ``kaiju`` software, you can get help text for the script
like this:

.. code-block:: bash
    
    create_gamhelio_ensemble.py --help
    usage: create_gamhelio_ensemble.py [-h] [--clobber] [--debug] [--verbose] ensemble_description_path

    Create an ensemble of gamhelio runs.

    positional arguments:
    ensemble_description_path
                            Path to .ini file describing ensemble (default: None)

    optional arguments:
      -h, --help            show this help message and exit
      --clobber             Overwrite existing files and directories (default: False).
      --debug, -d           Print debugging output (default: False).
      --verbose, -v         Print verbose output (default: False).

The ``ensemble_description_path`` option allows the user to specify an
``.ini``-format file which provides the parameter ranges for the ensemble.


The ensemble description file
-----------------------------

The ensemble description file is an ``.ini``-format file which provides the
parameter ranges for the ensemble, as well as the system-specific information
required to run the ensemble. In this way, ``create_gamhelio_ensemble.py`` is a
stripped-down, non-interactive version of ``makeitso-gamhelio.py``.

The ensemble description file starts with a section called ``[glparams]``
which defines the ranges for each parameter of the Gibson & Low CME model to
use in the ensemble. Next is the section ``[paths]``, which defines several
important directory and file locations. This is followed by a ``[norm]``
section which defines a few normalization parameters for the model. Last is
the ``[pbs]`` section, which defines all of the parameters needed to create
the ensemble as a set of PBS jobs.

.. note:: This script works on both ``derecho`` and ``pleiades``, but does not
  contain any system-specific information. All system-specific information
  must be included in the ``[pbs]`` section of the ensemble description file.

An example of an ensemble description file is provided below.

.. code-block:: ini

  [glparams]
  crot = 2095
  den_cme = 1500.
  dores = T
  dtres = 10.
  gl_bpar = 0.001, 0.002
  gl_lat = 7.
  gl_legsang = 10.
  gl_lon = 45.
  gl_orientation = 0.
  gl_topmorph = 1.75
  gl_vel_fh = 400.
  nres = 00000
  resid = helio
  t_cme = 500000.

  [paths]
  rundir = /glade/u/home/ewinter/cgs/runs/test/create_gamhelio_ensemble
  execpath = /glade/u/home/ewinter/cgs/runs/test/create_gamhelio_ensemble/gamhelio_mpi.x
  gridpath = /glade/u/home/ewinter/cgs/runs/test/create_gamhelio_ensemble/heliogrid.h5
  innerbcpath = /glade/u/home/ewinter/cgs/runs/test/create_gamhelio_ensemble/innerbc.h5
  restartpath = /glade/u/home/ewinter/cgs/runs/test/create_gamhelio_ensemble/helio_0002_0002_0002_0000_0000_0000.gam.Res.00000.h5

  [norm]
  gn0 = 200.
  gB0 = 1.e-3
  x0 = 6.956e10

  [pbs]
  hpc_system = derecho
  account_name = P28100045
  queue = main
  job_priority = economy
  select = 4
  ncpus = 128
  mpiprocs = 2
  ompthreads = 64
  modules = ncarenv/23.06, craype/2.7.20, intel/2023.0.0, ncarcompilers/1.0.0, cray-mpich/8.1.25, hdf5/1.12.2
  kaiju_install_dir = /glade/u/home/ewinter/cgs/aplkaiju/kaiju-private/ewinter-gamhelio_ensembles/kaiju-private
  kaipy_install_dir = /glade/u/home/ewinter/cgs/aplkaiju/kaipy-private/development/kaipy-private

The most important of these parameters are:

:crot: Number of Carrington rotation represented in the WSA FITS file

:dores: Set to ``T`` to launch the simulation from a restart file.

:dtres: Cadence for creation of restart files (simulated hours).

:gl_bpar: Maximum magnitude (Gauss) of the magnetic field in the Gibson-Low
  CME model.

:gl_lat: Latitude (degrees) of the CME emergence in the initial WSA map at
  21.5 R\ :sub:`S`.

:gl_lon: Longitude (degrees) of the CME emergence in the initial WSA map at
  21.5 R\ :sub:`S`.

:gl_orientation: Orientation angle (degrees) of the flux rope in the Gibson-
  Low CME model. 0° orientation means that the flux rope is in the
  equatorial plane.

:gl_topmorph: Sets the topology parameter of the Gibson-Low CME model.

:gl_vel_fh: Radial velocity of CME front (km/s).

:nres: Specify index of restart file to use for restart run.

:resid: Run ID of the restart file to use.

:t_cme: Time of CME emergence (simulated hours relative to start of
  simulation).


Creating the ensemble
---------------------

The ensemble definition file is passed to ``create_gamhelio_ensemble.py`` on
the command line:

.. code-block:: bash

  create_gamhelio_ensemble.py --verbose ensemble.ini

  Loading ensemble description from ensemble.ini.
  Creating ensemble parameter grid.
  Computing additional parameters.
  Creating run directories.
  Creating gamhelio XML files.
  Creating PBS files.
  Creating bash script for ensemble.

These messages illustrate the following steps in the ensemble creation
process:

* Read the ensemble description file.

* Create a grid of all possible parameter combinations.

* Compute additional parameters from the ensemble parameters.

* Create a directory for each ensemble member simulation.

* Create the input XML files for each ensemble member for use by
  ``gamhelio_mpi.x``.

* Create the PBS script for each ensemble member.

* Create a bash script that submits each ensemble member as a separate PBS
  job.

Once this step is complete, you can submit the entire set of ensemble members
for execution by running the script:

.. code-block:: bash

  bash ensemble.sh
