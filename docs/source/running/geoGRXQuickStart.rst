MAGE - With WACCM-X (GRX) Quick Start
===================================================

These instructions illustrate the process of running a geospace
simulation using the MAGE model coupled with WACCM-X. We call
this version of the model "GRX" ("GAMERA-RAIJU-WACCM-X") for brevity.
GRX runs are supported on ``derecho`` only.


Before you begin
------------------------------------------------

Setting environment variables
************************************************

*Source* (not *run*) the environment setup scripts for the ``kaiju`` software,
set ``WAXHOME`` to the path to your CAM clone, then activate the ``conda``
environment that holds ``kaipy`` (see :doc:`/python/index`). For example:

.. code-block:: bash

    source /path/to/your/kaiju-clone/scripts/setupEnvironment.sh
    export WAXHOME=/path/to/your/cam
    module load conda
    conda activate kaiju-3.12

.. note::

    The ``WAXHOME`` environment variable points to the CAM source code
    directory, which ``engage.py`` uses to create and build the WACCM-X
    cases. ``engage.py`` also needs ``ncks``, which is loaded with the
    module set for GRX runs on the :doc:`build page </building/buildDerecho_GRX>`.

The ``kaiju`` software (in the ``build_grx`` subdirectory) must be built and
the CAM source code cloned before you continue. Instructions are available
:doc:`here </building/buildDerecho_GRX>`.

Build ``voltron_mpi.x`` from this branch. It includes the ``xml-engage`` change: DragonKing reads
its conductance and precipitation settings from a ``<DRAGONKING>`` section of the deck, and
``makeitso``/``engage`` write only that section. They refuse a ``bin/voltron_mpi.x`` built without
the ``Kaiju/DRAGONKING`` reader, and every generated job re-checks the binary at launch and stops
with ``GATE`` exit 91 rather than run with default conductance.

Running a geospace simulation with MAGE
--------------------------------------------

The MAGE software needs several files in order to run. The detailed steps
for creating these files have been combined into a script called
``engage.py``. The script is provided in the ``kaiju`` code repository. More
information on ``engage.py`` is available
:doc:`here </makeitso/engage_grx>`.

You can see the options supported by ``engage.py`` by running it with the
``--help`` or ``-h`` command-line option. The full help text is shown
:doc:`here </makeitso/engage>`.

.. code-block:: bash

    engage.py --help

Running ``engage.py``
************************************************

For this example, we will run the code on ``derecho``, and use the default
``BASIC`` mode, which requires the minimum amount of input from the user. At
each prompt, you can either type in a value, or hit the :kbd:`Return` key to
accept the default value (shown in square brackets at the end of the prompt).

1. First we will create a directory for this run in your scratch space.

.. code-block:: bash

    mkdir quickstart_grx
    cd quickstart_grx

2. Run ``engage.py`` with ``--component grx``.

.. code-block:: bash

    engage.py --component grx

This opens an arrow-key wizard that asks the questions below, one at a time: type a
value or press :kbd:`Enter` to accept the default and go to the next one.
The keys are listed at the bottom of the screen: ``↑ back   ↓/Enter next   Ctrl-S finish   Ctrl-R reset   Ctrl-C cancel``.
This is the recommended way to run ``engage.py``.
The wizard asks the following questions. They are shown here, for reference,
as the plain prompts that ``--no-tui`` prints (the fallback when your terminal
cannot show the wizard; the tool also falls back by itself):

.. code-block:: bash

    Name to use for PBS job(s) [geospace]:
    Start date for simulation (yyyy-mm-ddThh:mm:ss) [2013-03-16T16:00:00]:
    Stop date for simulation (yyyy-mm-ddThh:mm:ss) [2013-03-18T12:00:00]:
    GAMERA grid type (Q) [Q]:
    Name of HPC system (derecho) [derecho]:
    PBS account name (your projects: <YOUR_ACCOUNT_HERE>): <YOUR_ACCOUNT_HERE>
    Run directory [.]:
    Path to kaiju installation [<YOUR_KAIJUHOME_HERE>]:
    Path to kaiju build directory [<YOUR_KAIJUHOME_HERE>/build_grx]:
    PBS queue name (develop|main) [main]:
    Job priority (premium|regular|economy) [economy]:
    Requested wall time for each PBS job segment (HH:MM:SS) [12:00:00]:
    Suggested segment for walltime 12:00:00 on main: 14400 s (4 h); valid values (seconds): [14400 (4 h)]
    Segment length in simulated seconds [14400]:
    Extend TFIN by dtCouple - 1 seconds (T|F) [T]:

.. warning::

    Make sure to set ``Path to kaiju build directory`` to the directory where you built the
    ``voltron_mpi.x`` executable with the module set for GRX runs. If you followed the build
    instructions, this should be the ``build_grx`` subdirectory. This is required for running
    the model in GRX mode.

``engage.py`` will then prompt you for the following additional information from ``makeitso``:

.. code-block:: bash

    Instructions:
    -> Default selected input parameter is given in GREEN
    -> Warnings and information are given in YELLOW
    -> Errors are given in RED
    -> Valid values (if any) are given in brackets eg. (value1 | value2 | value3)


    Run Options:
    Do you have an existing solar wind file to use? (Y|y|N|n) [N]:
    (GAMERA) Relative path to HDF5 file containing solar wind data [bcwind.h5]:
    (VOLTRON) File output cadence in simulated seconds [60.0]:

After these inputs, the script fetches data from CDAWeb for the specified time
range to use in the solar wind boundary condition file, and generates the XML
files for the run, as well as a grid file for use in the model.

You should see output similar to this:

.. code-block:: bash

    Generating Quad LFM-style grid ...
    ..
    ..
    Template creation complete!


    makeitso: wrote 13 GRX voltron config files to <YOUR_RUN_DIRECTORY_HERE> (2 GR lead-in decks in gr_spinup/)

``engage.py`` will then prompt you for the following additional information from ``waxrun``:

.. code-block:: bash

    Instructions:
    -> Default Selected input parameter is given in GREEN
    -> Warnings and Information are given in YELLOW
    -> Errors are given in RED
    -> Valid values (if any) are given in brackets eg. (value1 | value2 | value3)
    -> Enter '?' for any input parameter to get a detailed description
    -> 17 option(s) already set by engage are not asked again


    Run Options:
    CAM source tree (your clone) [<YOUR_CAM_TREE_HERE>]:
    The default series covers 2003-12-11 to 2004-05-01 (winter and spring storms only).
    Start file directory: a CESM rest/ directory, a series root or one cam.i file (native ne30pg3 L145) [/glade/campaign/hao/itmodel/liuh/csm/waccmx_ne30_2004_cmip6/rest]:
    Generate nudging data? (yes/no) [Y]:
    nearest restart: 2004-03-06T00:00:00, re-dated to 2013-03-06 (256 h of spinup)
    available: 2003-12-11T00:00:00, 2003-12-21T00:00:00, 2004-01-06T00:00:00, 2004-01-11T00:00:00, 2004-01-21T00:00:00, 2004-02-05T00:00:00, 2004-02-15T00:00:00, 2004-02-25T00:00:00, 2004-03-06T00:00:00, 2004-03-16T00:00:00, 2004-03-26T00:00:00, 2004-04-05T00:00:00, 2004-04-15T00:00:00, 2004-04-25T00:00:00, 2004-05-01T00:00:00
    Start file date (yyyy-mm-ddThh:mm:ss) [2004-03-06T00:00:00]:
    Spinup start date (yyyy-mm-ddThh:mm:ss) [2013-03-06T00:00:00]:

After these inputs, the script creates and builds the two WACCM-X cases, and
generates the PBS files for the run. Building WACCM-X compiles CESM, which
takes a while.

You should see output similar to this (the output on this page is from a
short run named ``ltgrx``, with two one-hour coupled segments):

.. code-block:: bash

    ..
    ..
    waxrun: parameters saved to <YOUR_RUN_DIRECTORY_HERE>/waxrun_parameters.json
    Creating 4 GRX segment PBS scripts (ic=1, maxhours=3)
    GRX voltron configs = ['<YOUR_RUN_DIRECTORY_HERE>/gr_spinup/ltgrx-SPINUP.xml', '<YOUR_RUN_DIRECTORY_HERE>/gr_spinup/ltgrx-WARMUP-01.xml', '<YOUR_RUN_DIRECTORY_HERE>/ltgrx-01.xml', '<YOUR_RUN_DIRECTORY_HERE>/ltgrx-02.xml']
    GRX spinup submits  = ['<YOUR_RUN_DIRECTORY_HERE>/ltgrx_spinup_submit_01.sh']
    GRX coupled case    = ['<YOUR_RUN_DIRECTORY_HERE>/wax_case/ltgrx', '<YOUR_RUN_DIRECTORY_HERE>/wax_case/ltgrx/user_nl_cam']
    GRX segment PBS     = ['<YOUR_RUN_DIRECTORY_HERE>/gr_spinup/ltgrx-SPINUP.pbs', '<YOUR_RUN_DIRECTORY_HERE>/gr_spinup/ltgrx-WARMUP-01.pbs', '<YOUR_RUN_DIRECTORY_HERE>/ltgrx-01.pbs', '<YOUR_RUN_DIRECTORY_HERE>/ltgrx-02.pbs']
    GRX submit script   = <YOUR_RUN_DIRECTORY_HERE>/ltgrx_grx.sh

Looking at files generated by ``engage.py``
**************************************************************

Your run directory now contains the following files and directories:

.. code-block:: bash

    ls
    atm_in                  init_generated_files       nudging_gen
    bcwind.h5               lfmQ.h5                    nuopc.runconfig
    CASEROOT                lnd_in                     nuopc.runseq
    cesm                    ltgrx-01.pbs               OMNI_HRO_1MIN.txt_bxFit.png
    docn_in                 ltgrx-01.xml               OMNI_HRO_1MIN.txt.png
    docn.streams.xml        ltgrx-02.pbs               raijuconfig.h5
    drv_flds_in             ltgrx-02.xml               shared
    drv_in                  ltgrx_grx.sh               voltron_mpi.x
    engage_parameters.json  ltgrx_spinup_submit_01.sh  wax_case
    fd.yaml                 makeitso_parameters.json   waxrun_parameters.json
    gr_spinup               met_data_ltgrx
    ice_in                  mosart_in

Among them:

* ``engage_parameters.json``, ``makeitso_parameters.json``, ``waxrun_parameters.json``
    These are the JSON files that contain the parameters for the simulation. They
    are generated by the ``engage.py`` script with all the parameters required to run the
    simulation.
* ``geospace-**.pbs`` and ``geospace-**.xml``
    These are the PBS scripts and the XML files for the segments of the GRX
    coupled model. The ``-01``, ``-02``, etc. suffixes indicate the segment
    number, and the segments are run in order.
* ``geospace_grx.sh``
    This is the script that submits all the jobs of the run.
* ``gr_spinup``
    This directory holds the ``geospace-SPINUP.*`` and ``geospace-WARMUP-01.*``
    segments, which run the GAMERA RAIJU model before the coupled model execution.
* ``wax_case`` and ``cesm``
    These directories hold the two WACCM-X cases (the WACCM-X spin-up and the
    coupled case) and their build.
* ``nudging_gen`` (when nudging data are missing)
    This directory holds the job that generates the missing nudging data for WACCM-X.

The ``gr_spinup`` directory contains:

.. code-block:: bash

    ls gr_spinup
    bcwind.h5         ltgrx-SPINUP.xml     raijuconfig.h5
    lfmQ.h5           ltgrx-WARMUP-01.pbs  voltron_mpi.x
    ltgrx-SPINUP.pbs  ltgrx-WARMUP-01.xml

Submitting the GRX model run
*************************************

Finally, submit the model run using the script generated by ``engage.py``.
You will see the resulting PBS job ID (your job ID will differ from what is
shown below). The job IDs are also saved in ``geospace_grx.jobids``.

.. code-block:: bash

    bash ltgrx_grx.sh
    nud=7716687.desched1
    submit_jobs case.run
    Submit job case.run
    nudging=7716687.desched1
    spinup=7716689
    gr_spinup=7716690.desched1
    gr_warmup=7716691.desched1
    coupled=7716692.desched1
    To resubmit the coupled chain after a fix (remove CESM leftovers first, or rerun with CLOBBER=1):
        cd <YOUR_RUN_DIRECTORY_HERE> && qsub ltgrx-01.pbs

Replaying the run
*************************************

The options you entered are saved in ``engage_parameters.json``,
``makeitso_parameters.json`` and ``waxrun_parameters.json``. To generate the
same run again, run this command in the run directory. It also recreates and
rebuilds the two WACCM-X cases:

.. code-block:: bash

    engage.py --component grx --clobber --recreate-cases -eo engage_parameters.json -mo makeitso_parameters.json -wo waxrun_parameters.json

Understanding the output files of the GRX model run
**************************************************************************

When the jobs are complete, you will see many new HDF5 files in your run
directory, along with PBS housekeeping files and logs. The most important
files are (repeated upper-case letters in the names represent integer
strings):

* ``geospace_LLLLL_MMMMM_NNNNN_IIIII_JJJJJ_KKKKK.gam.h5``

    These files contain the core MHD variables from the simulation, computed
    by the GAMERA portion of the MAGE model. The strings ``LLLLL``, ``MMMMM``,
    and ``NNNNN`` contain the number of subsections of the ``X``, ``Y``, and
    ``Z`` dimensions used to divide the domain among MPI ranks. The strings
    ``IIIII``, ``JJJJJ``, and ``KKKKK`` represent the MPI rank index along
    each dimension.

* ``geospace.mix.h5``

    This file contains the results from the
    `REMIX <https://cgs.jhuapl.edu/Models/remix.php>`_ portion of the
    `MAGE <https://cgs.jhuapl.edu/Models>`_ model.

* ``geospace.raiju.h5``

    This file contains the results from the
    RAIJU portion of the `MAGE <https://cgs.jhuapl.edu/Models>`_ model.

* ``geospace.cam.h0i.*.nc``, ``geospace.cam.h1i.*.nc``,
  ``geospace.cam.h2a.*.nc``, ``geospace.cam.h2i.*.nc``,
  ``geospace.cam.h3i.*.nc`` and ``geospace.cam.h4i.*.nc``

    These are the output files from the WACCM-X portion of the
    `MAGE <https://cgs.jhuapl.edu/Models>`_ model. The ``ltgrx`` run wrote,
    for example, ``ltgrx.cam.h0i.2013-03-16-21600.nc`` (at the start of the
    coupled run) and ``ltgrx.cam.h4i.2013-03-16-21660.nc`` (every 60
    seconds).

* ``geospace.cam.r.*.nc``, ``geospace.cam.rs.*.nc``,
  ``geospace.cam.rh0.*.nc`` and ``geospace.cam.rh2.*.nc``

    These are the WACCM-X checkpoint files, written at the end of each
    segment, for example ``ltgrx.cam.r.2013-03-16-25200.nc``.

* ``geospace_*.gam.Res.RRRRR.h5``

    These are checkpoint files generated during the simulation which can be
    used as restart points for future simulations.
