Makeitso -- Use for MAGE without TIEGCM
=============================================


Introduction
------------

The Python script ``makeitso.py`` was developed to simplify the process of
configuring and running GR MAGE (that is, the geospace application  of
the ``kaiju`` software without TIEGCM.) It
provides an interactive, prompt-driven interface to specify all of the
parameters needed for a model run.

The ``makeitso.py`` script can operate in one of three different modes:
``BASIC``, ``INTERMEDIATE``, or ``EXPERT``. Each mode provides access to
a subset of the ``kaiju`` parameters. The ``BASIC`` mode requires
the user to provide the minimum set of parameters needed to specify a model
run, such as the run ID, and the simulation time period. The ``INTERMEDIATE``
mode allows the user to specify all of the parameters from the ``BASIC`` mode,
as well as a wider set of run parameters, such as non-standard file
locations and some MHD parameters. The ``EXPERT`` mode provides access to all
of the user-adjustable parameters from the ``kaiju`` software. When finished,
the script generates the files needed to run a magnetosphere model, and saves
all options in a convenient JSON file so that the run can be repeated at a
later date.

The entry point for MAGE runs is the ``engage.py`` script, which starts
``makeitso.py`` for a GR run (``engage.py --component gr``) and asks the same
questions. More information on ``engage.py`` is available
:doc:`here </makeitso/engage>`. This page describes ``makeitso.py`` used
directly.


Running the ``makeitso.py`` script
----------------------------------

The ``makeitso.py`` script is provided as part of the ``kaiju`` software. It
is found at ``$KAIJUHOME/scripts/makeitso/makeitso.py``, where ``$KAIJUHOME``
is the location of your ``kaiju`` software tree. After configuring your
``kaiju`` software, you can get help text for the script like this:

.. code-block:: bash
    :force:

    makeitso.py --help
    usage: makeitso.py [-h] [--clobber] [--debug] [--mode MODE] [--options_path OPTIONS_PATH]
                       [--rederive GROUP[,GROUP]] [--verbose] [--no-tui] [--tui]

    Interactive script to prepare a MAGE magnetosphere model run

    options:
      -h, --help            show this help message and exit
      --clobber             Overwrite existing options file (default: False).
      --debug, -d           Print debugging output (default: False).
      --mode MODE           User mode: BASIC|INTERMEDIATE|EXPERT, case-insensitive (default: BASIC).
      --options_path OPTIONS_PATH, -o OPTIONS_PATH
                            Path to JSON file of options (default: None).
      --rederive GROUP[,GROUP]
                            With -o: re-derive these groups; '--rederive list' names them.
      --verbose, -v         Print verbose output (default: False).
      --no-tui              Ask the options as linear prompts instead of the wizard.
      --tui                 The wizard (the default; kept for old command lines).

    The run is written to pbs.run_directory; a segmented run's GR spin-up runs in gr_spinup/. Layout
    and replay: docs/source/makeitso/makeitso.rst.

The ``--options_path`` option allows the user to specify an existing JSON file
from a previous run of ``makeitso.py`` so that the entire process of model
generation can be automated. The ``--mode`` option specifies the user mode to
run in, with ``BASIC`` being the default.


An example in ``BASIC`` mode
----------------------------

This section provides an annotated example session of ``makeitso.py`` running
in the default ``BASIC`` mode on the ``aitken`` supercomputer. Run:

.. code-block:: bash

  makeitso.py

This opens an arrow-key wizard that asks the questions of this session, one at a time: type a
value or press :kbd:`Enter` to accept the default and go to the next one.
The keys are listed at the bottom of the screen: ``↑ back   ↓/Enter next   Ctrl-S finish   Ctrl-R reset   Ctrl-C cancel``.
This is the recommended way to run ``makeitso.py``.
The wizard asks the following questions. They are shown here, for reference,
as the plain prompts that ``--no-tui`` prints (the fallback when your terminal
cannot show the wizard; the tool also falls back by itself):

.. code-block:: bash

  Name to use for PBS job(s) [geospace]:

Enter an identifying string to use for your model run. This name will be used
as the basis for most of the files created by ``makeitso.py`` and the
``kaiju`` software. The default name is ``geospace``.

.. code-block:: bash

  Do you have an existing solar wind file to use? (Y|y|N|n) [N]:

If you already have a file containing solar wind data to use for the boundary
conditions of your simulation, enter ``Y``, and you will then be prompted for
the path to the file. If you don't have the file, enter ``N`` and you will be
prompted for the date range to use.

.. code-block:: bash

  Start date for simulation (yyyy-mm-ddThh:mm:ss) [2001-06-01T23:00:00]: 2016-08-09T09:00:00
  Stop date for simulation (yyyy-mm-ddThh:mm:ss) [2001-06-02T01:00:00]: 2016-08-09T11:00:00

Enter the start and stop date and time for the solar wind data you want to
use. The required data will be fetched from CDAWeb, and converted into a
format usable by the ``kaiju`` software.

.. code-block:: bash

  GAMERA grid type (D|Q|O|H) [Q]:

The codes represent double- (``D``), quad- (``Q``), oct- (``O``) and
hex- (``H``) resolutions in the LFM grid used in the ``kaiju`` software.

.. code-block:: bash

  Name of HPC system (derecho|aitken|linux|perlmutter) [derecho]: aitken

Enter ``aitken`` to prepare the run for ``aitken``. The selection you make
here will customize the remaining prompts for the selected system.

.. code-block:: bash

  PBS account name [<YOUR_LOGIN_NAME_HERE>]:

Enter the PBS account to charge your jobs to. On ``aitken`` the default is
your login name. On ``derecho`` there is no default, and the prompt lists your
projects.

.. code-block:: bash

  Run directory [.]:

Specify the directory that you wish to perform the simulation in. The
directory will contain all of the files generated by ``makeitso.py``. A
relative path is taken relative to the directory you run ``makeitso.py``
from, and the directory is created if it does not exist.

.. code-block:: bash

  Path to kaiju installation [<YOUR_KAIJUHOME_HERE>]:
  Path to kaiju build directory [<YOUR_KAIJUHOME_HERE>/build_mpi]:

Enter the paths to the location of your ``kaiju`` code, and the location of
your ``kaiju`` build directory. The defaults are ``$KAIJUHOME`` and its
``build_mpi`` subdirectory.

.. code-block:: bash

  PBS queue name (low|normal|long|debug|devel) [normal]:

Select a PBS queue to use on the selected supercomputer.

.. code-block:: bash

  Requested wall time for each PBS job segment (HH:MM:SS) [08:00:00]:
  NOTE: no wall-clock rate is recorded for a gr segment on grid Q: the walltime rule uses 3 h per simulated hour (estimate), conservatively. (aitken: config/machines.yaml rates)

Specify the wall clock time to request for your job (or each segment, if you
split your job into multiple segments). On ``aitken`` the script then prints
a ``NOTE`` about the wall time estimate for the selected grid.

.. code-block:: bash

  PBS group list:

On ``aitken`` you can enter a PBS group list for your jobs. This prompt has no
default, and can be left empty.

.. code-block:: bash
  :force:

  Do you want to split your job into multiple segments? (Y|y|N|n) [N]:
  WARNING: at grid Q on aitken (rate 3 h/h estimate x margin 2) walltime 08:00:00 is expected to hold a segment of at most 1 h: not refused (no measured rate to refuse on); watch the first job's wall time.

Enter ``Y`` here if you want to split your simulation into multiple PBS jobs.
This will allow you to run long simulations that are chained together, with
each using the results of the previous job as a starting point. If you enter
``Y``, you will be prompted for segment length (in simulated time). The
``WARNING`` compares the wall time with the expected run time of a segment.

.. code-block:: bash

  (GAMERA) Relative path to HDF5 file containing solar wind data [bcwind.h5]:

This is the name of the solar wind file in the run directory. When the data
are fetched from CDAWeb, ``makeitso.py`` writes them to this file.

.. code-block:: bash

  (VOLTRON) File output cadence in simulated seconds [60.0]:

How often (in simulated seconds) the ``kaiju`` software should output results
during the course of the simulation.

The script then runs several additional tools to prepare the files needed for
your simulation. This includes input files, as well as PBS job scripts for
your simulation, and a bash shell script to submit the PBS jobs.

.. code-block:: bash

  Generating Quad LFM-style grid ...

  Output: lfmQ.h5
  Size: (96,96,128)
  Inner Radius: 2.000000
  Sunward Outer Radius: 30.000000
  Tail Outer Radius: 322.511578
  Low-lat BC: 45.000000
  Ring params:
  <ring gid="lfm" doRing="T" Nr="8" Nc1="8" Nc2="16" Nc3="32" Nc4="32" Nc5="64" Nc6="64" Nc7="64" Nc8="64"/>

  Writing to lfmQ.h5
  05-Oct-26 02:38:02: /path/to/your/conda/envs/kaiju-3.12/lib/python3.12/site-packages/spacepy/time.py:2448: UserWarning: Leapseconds may be out of date. Use spacepy.toolbox.update(leapsecs=True)
    _read_leaps()

  Retrieving f10.7 data from CDAWeb
  Retrieving solar wind data from CDAWeb
          Using Bx fields
  Bx Fit Coefficients are  [-3.78792744 -0.77915822 -1.0774984 ]
  Saving "OMNI_HRO_1MIN.txt_bxFit.png"
  Converting to Gamera solar wind file
          Found 21 variables and 120 lines
          Offsetting from LFM start ( 0.00 min) to Gamera start ( 0.00 min)
  Saving "OMNI_HRO_1MIN.txt.png"
  Writing Gamera solar wind to bcwind.h5
  Making new raijuconfig.h5, destroying pre-existing file if there
  Stamping file with git hash and branch, and script args
  Adding waveModel to raijuconfig.h5
  Reading /path/to/your/kaipy/kaipy/raiju/waveModel/chorus_polynomial.txt
  Added name: 0_Plasmasphere with flavor: 0
  Added name: Hot Electrons with flavor: 1
  Added name: Hot Protons with flavor: 2
  Adding Species to raijuconfig.h5
  Adding params used to generate lambda distribution as root attribute


  Template creation complete!


  The PBS scripts ['./geospace.pbs'] have been created, each with a corresponding XML file. To submit the jobs with the proper dependency (to ensure each segment runs in order), please run the script geospace_pbs.sh like this:
  bash geospace_pbs.sh
  The options are saved in makeitso_parameters.json; replay them with makeitso.py -o makeitso_parameters.json.

When finished, the script creates the file ``makeitso_parameters.json`` in
the run directory. This file contains a record of all of the parameters used
in your simulation. This file can be passed back to ``makeitso.py`` with the
``-o`` option in a subsequent session to repeat the simulation, and also
provides a convenient starting point for minor tweaks to your simulation
parameters. In the run directory, add ``--clobber`` so that the file can be
overwritten:

.. code-block:: bash

  makeitso.py --clobber -o makeitso_parameters.json

A file passed back with ``-o`` is replayed as it was saved. If you edit the
dates, the grid, the job name or the machine settings in it, add the
``--rederive`` option to re-derive the values that depend on them
(``makeitso.py --rederive list`` prints the groups).


Additional parameters in ``INTERMEDIATE`` and ``EXPERT`` mode
-------------------------------------------------------------

Many more parameters are available in ``INTERMEDIATE`` and ``EXPERT`` modes.
These parameters are documented in the file ``option_descriptions.json``,
which is stored in the same directory as the ``makeitso.py`` script.
