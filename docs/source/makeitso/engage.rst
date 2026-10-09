Engage -- Use for MAGE runs (GR, GRT and GRX)
=============================================


Introduction
------------

The Python script ``engage.py`` was developed to simplify the process of
configuring and running MAGE (that is, the geospace application of the
``kaiju`` software). It is the entry point for three types of model run:

  * GR MAGE, the geospace application without TIEGCM ("GAMERA-RAIJU").

  * GRT MAGE, the geospace application coupled with TIEGCM
    ("GAMERA-TIEGCM-RAIJU").

  * GRX MAGE, the geospace application coupled with WACCM-X
    ("GAMERA-RAIJU-WACCM-X"), available on ``derecho`` only.

It provides an interactive, prompt-driven interface to specify all of the
parameters needed for a model run.

``engage.py``, ``makeitso.py`` and ``tiegcmrun.py`` ask their options in an
arrow-key wizard by default, which is the recommended way to use them. The
walkthrough pages show the same questions as the plain prompts that
``--no-tui`` prints, for reference. The first screen of the GRT wizard looks
like this:

.. code-block:: none

    engage (GRT) wizard — BASIC    step 1/16    [simulation]
     job_name   [BASIC]                                       │  ▶ job_name               geospace
     Name to use for PBS job(s)                               │    start_date             2001-06-01T2
     default (default):  geospace                             │    stop_date              2001-06-02T0
                                                              │     … 13 more
     › geospace                                               │

     (1 field(s) need fixing)   ↑ back   ↓/Enter next   Ctrl-S finish   Ctrl-R reset   Ctrl-C cancel

The ``engage.py`` script is a wrapper around the ``makeitso.py`` script,
TIE-GCM's ``tiegcmrun`` script and the WACCM-X ``waxrun`` script, which are
used to prepare the necessary files for a MAGE model run.

  - For more details on the ``makeitso.py`` script, see the :doc:`makeitso </makeitso/makeitso>` documentation.

  - For more details on the TIE-GCM ``tiegcmrun`` script, see the `tiegcmrun <https://tiegcm-docs.readthedocs.io/en/latest/tiegcm/quickstart>`_ documentation.

  - For more details on the ``waxrun`` parameters, see the :doc:`GRX </makeitso/engage_grx>` documentation.

The ``engage.py`` script can operate in one of three different modes:
``BASIC``, ``INTERMEDIATE``, or ``EXPERT``. Each mode provides access to
a subset of the ``kaiju``, TIE-GCM and WACCM-X parameters.

  * The ``BASIC`` mode
    Requires the user to provide the minimum set of parameters needed to specify a model
    run, such as the run ID, and the simulation time periods.

  * The ``INTERMEDIATE`` mode
    Allows the user to specify all of the
    parameters from the ``BASIC`` mode, as well as a wider set of run parameters,
    such as non-standard file locations and some MHD parameters and TIE-GCM parameters.

  * The ``EXPERT`` mode
    Provides access to all of the user-adjustable
    parameters from the ``kaiju``, TIE-GCM and WACCM-X software.

When finished, the script generates the files needed to run a magnetosphere model, and saves
all options in convenient JSON files so that the run can be repeated at a
later date.


Running the ``engage.py`` script
----------------------------------

The ``engage.py`` script is provided as part of the ``kaiju`` software. It
is found at ``$KAIJUHOME/scripts/engage/engage.py``, where ``$KAIJUHOME``
is the location of your ``kaiju`` software tree. After configuring your
``kaiju`` software, you can get help text for the script like this:

.. code-block:: bash

    engage.py --help
    usage: engage.py [-h] [--clobber] [--recreate-cases] [--debug] [--mode MODE]
                     [--component {grt,grx,tiegcm,gr}] [--engage_options_path ENGAGE_OPTIONS_PATH]
                     [--makeitso_options_path MAKEITSO_OPTIONS_PATH]
                     [--tiegcm_options_path TIEGCM_OPTIONS_PATH] [--rederive GROUP[,GROUP]]
                     [--verbose] [--waxrun_options_path WAXRUN_OPTIONS_PATH] [--dry-run]
                     [--skip-build] [--phase {spinup,coupled,both}] [--no-tui] [--tui]

    Configure and generate a MAGE run (GRT, GRX or GR) or a standalone TIE-GCM run.

    options:
      -h, --help            show this help message and exit
      --clobber             Overwrite the run's existing replay files (default: False).
      --recreate-cases      GRX: delete and recreate existing WACCM-X cases (default: False).
      --debug, -d           Print debugging output (default: False).
      --mode MODE           User mode: BASIC|INTERMEDIATE|EXPERT, case-insensitive (default: BASIC).
      --component {grt,grx,tiegcm,gr}, -m {grt,grx,tiegcm,gr}
                            Run type: coupled grt or grx, standalone tiegcm or gr (default: grt).
      --engage_options_path ENGAGE_OPTIONS_PATH, -eo ENGAGE_OPTIONS_PATH
                            Path to engage JSON file of options (default: None).
      --makeitso_options_path MAKEITSO_OPTIONS_PATH, -mo MAKEITSO_OPTIONS_PATH
                            Path to makeitso JSON file of options (default: None).
      --tiegcm_options_path TIEGCM_OPTIONS_PATH, -to TIEGCM_OPTIONS_PATH
                            Path to tiegcm JSON file of options (default: None).
      --rederive GROUP[,GROUP]
                            With -eo/-mo/-to: re-derive these groups; '--rederive list' names them.
      --verbose, -v         Print verbose output (default: False).
      --waxrun_options_path WAXRUN_OPTIONS_PATH, -wo WAXRUN_OPTIONS_PATH
                            GRX: path to waxrun JSON file of options (default: None).
      --dry-run             GRX: write the case specs as plain files; run no CIME script.
      --skip-build          GRX: skip case.build; every cesm.exe must already exist.
      --phase {spinup,coupled,both}
                            GRX: which WACCM-X case(s) to generate (default: both).
      --no-tui              Ask the options as linear prompts instead of the wizards.
      --tui                 The wizards (the default; kept for old command lines).

    The run is written to pbs.run_directory; a segmented run's GR lead-in runs in gr_spinup/. Layout
    and replay: docs/source/makeitso/engage.rst and the GR, GRT and GRX pages.

The ``--component`` option selects the type of model run. Each type has a
page describing its prompts, and a quick start guide that walks through a
complete run:

.. list-table::
   :header-rows: 1

   * - ``--component``
     - Run type
     - Prompts
     - Quick start
   * - ``gr``
     - MAGE - Without TIEGCM (GR)
     - :doc:`here </makeitso/engage_gr>`
     - :doc:`here </running/geoGRQuickStart>`
   * - ``grt`` (the default)
     - MAGE - With TIEGCM (GRT)
     - :doc:`here </makeitso/engage_grt>`
     - :doc:`here </running/geoGRTQuickStart>`
   * - ``grx``
     - MAGE - With WACCM-X (GRX), ``derecho`` only
     - :doc:`here </makeitso/engage_grx>`
     - :doc:`here </running/geoGRXQuickStart>`

For a GRX run, the prompts offer GRX defaults (grid ``Q``, 4-hour segments
and dates in March 2013), and the ``waxrun`` parameters are requested after
the ``engage`` and ``makeitso`` parameters. The GRX-only options are described
:doc:`here </makeitso/engage_grx>`.

The ``--mode`` option specifies the user mode to run in, with ``BASIC`` being
the default. The ``-eo``, ``-mo``,
``-to`` and ``-wo`` options allow the user to specify existing JSON files so
that the entire process of model generation can be automated. More info on
this given below.


Additional parameters in ``INTERMEDIATE`` and ``EXPERT`` mode
-------------------------------------------------------------

Many more parameters are available in ``INTERMEDIATE`` and ``EXPERT`` modes.
These parameters are documented in the files
``$KAIJUHOME/scripts/engage/option_engage_descriptions.json`` and
``$KAIJUHOME/scripts/makeitso/option_descriptions.json``, and, for the
TIE-GCM and WACCM-X parameters, in
``$TIEGCMHOME/tiegcmrun/options_description.json`` and
``$KAIJUHOME/scripts/waxrun/config/options_description.json``.


Using JSON files for ``engage.py``
-------------------------------------------------------------

The ``engage.py`` script can also be run in a non-interactive mode, where it
reads JSON files containing the parameters for the simulation. This allows
you to automate the process of running the simulation, and to easily repeat
the simulation with the same parameters.

Each type of run takes its own set of JSON files:

  * GR: ``-mo``
  * GRT: ``-eo``, ``-mo`` and ``-to``
  * GRX: ``-eo``, ``-mo`` and ``-wo``

The files of a run type must be given together; ``engage.py`` then asks no
questions. If only some of them are given, ``engage.py`` stops and names the
missing ones.

The four files contain:

  * ``engage_options_path`` (``-eo``)
    This is the path to the JSON file containing the parameters for the
    ``engage.py`` script. It contains the parameters that are shared by the
    coupled run, such as the run ID, start and stop dates, segment length,
    GAMERA grid type and PBS account.
  * ``makeitso_options_path`` (``-mo``)
    This is the path to the JSON file containing the parameters for the
    ``makeitso.py`` script. It contains the parameters that are specific to the
    ``makeitso.py`` script, such as the solar wind file and the output cadence.
    For a GR run it is the only file, so it also contains the run ID, start and
    stop dates, GAMERA grid type and PBS account.
  * ``tiegcm_options_path`` (``-to``)
    This is the path to the JSON file containing the parameters for the
    ``tiegcmrun`` script. It contains the parameters that are specific to the
    TIE-GCM simulation, such as the TIE-GCM executables, F10.7 flux level, and so on.
  * ``waxrun_options_path`` (``-wo``)
    This is the path to the JSON file containing the parameters for the
    ``waxrun`` script. It contains the parameters that are specific to the
    WACCM-X simulation, such as the CAM source tree and the start file
    directory.

A file only needs to contain the parameters you want to set. Every parameter
it leaves out is set to the default value of its prompt, and ``engage.py``
prints a ``NOTE`` saying which parameters it filled in (or, when there are
many, how many); the full set is saved in the run directory, as described
below. The start and stop dates, the GAMERA grid type (GR and GRT), the PBS
account (on ``derecho``) and, for GRX, the CAM source tree (``srcroot``,
unless ``WAXHOME`` is set) have no default, and must be in the files.

To run the ``engage.py`` script in non-interactive mode, copy the template
files for your run type (below) into your run directory, edit them, and run
``engage.py`` from that directory. For a GR run:

.. code-block:: bash

    engage.py --component gr -mo makeitso_input.json

For a GRT run:

.. code-block:: bash

    engage.py -eo engage_input.json -mo makeitso_input.json -to tiegcmrun_input.json

For a GRX run:

.. code-block:: bash

    engage.py --component grx -eo engage_input.json -mo makeitso_input.json -wo waxrun_input.json

Here are templates for the JSON files:

  - GR:

    - derecho: :download:`makeitso_input.json <engage_template/derecho/gr/makeitso_input.json>`
    - aitken: :download:`makeitso_input.json <engage_template/aitken/gr/makeitso_input.json>`

  - GRT on derecho:

    - :download:`engage_input.json <engage_template/derecho/engage_input.json>`
    - :download:`makeitso_input.json <engage_template/derecho/makeitso_input.json>`
    - :download:`tiegcmrun_input.json <engage_template/derecho/tiegcmrun_input.json>`

  - GRT on aitken:

    - :download:`engage_input.json <engage_template/aitken/engage_input.json>`
    - :download:`makeitso_input.json <engage_template/aitken/makeitso_input.json>`
    - :download:`tiegcmrun_input.json <engage_template/aitken/tiegcmrun_input.json>`

  - GRX on derecho:

    - :download:`engage_input.json <engage_template/derecho/grx/engage_input.json>`
    - :download:`makeitso_input.json <engage_template/derecho/grx/makeitso_input.json>`
    - :download:`waxrun_input.json <engage_template/derecho/grx/waxrun_input.json>`

These JSON files can be used as a starting point for your own simulations. You will
need to modify certain parameters in them:

  - engage_input.json (GRT and GRX), and makeitso_input.json for a GR run:

    - account_name: Your PBS account. On ``aitken`` this is not needed: your
      login name is used.
    - start_date: The start date of your simulation.
    - stop_date: The stop date of your simulation.
    - segment_duration (GRT and GRX): The duration of each segment in simulated seconds.
    - gamera_grid_type (GR and GRT): The GAMERA grid type to use.

  - makeitso_input.json (GRT and GRX):

    - Nothing needs to be changed. It sets the solar wind file option and the
      output cadence, which you can modify if needed.

  - tiegcmrun_input.json:

    - modelexe: The path to the TIE-GCM standalone executable (``tiegcm.exe``).
    - coupled_modelexe: The path to the TIE-GCM coupled executable (``tiegcm.x``).
    - A relative path is read relative to the directory of the JSON file.

  - waxrun_input.json:

    - srcroot: The path to your CAM source tree.

When finished, ``engage.py`` saves all of the parameters used in your
simulation in the run directory: ``makeitso_parameters.json`` for a GR run;
``engage_parameters.json``, ``makeitso_parameters.json`` and
``tiegcmrun_parameters.json`` for a GRT run; and ``engage_parameters.json``,
``makeitso_parameters.json`` and ``waxrun_parameters.json`` for a GRX run.
These files can be passed back to ``engage.py`` in the same way to repeat the
simulation, and also provide a convenient starting point for minor tweaks to
your simulation parameters. The command to repeat each type of run is given on
its page: :doc:`GR </makeitso/engage_gr>`, :doc:`GRT </makeitso/engage_grt>`
and :doc:`GRX </makeitso/engage_grx>`.
