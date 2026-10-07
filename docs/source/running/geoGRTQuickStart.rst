MAGE - With TIEGCM (GRT) Quick Start
===================================================

These instructions illustrate the process of running a geospace
simulation using the MAGE model coupled with TIEGCM. We call
this version of the model "GRT" ("GAMERA-TIEGCM-RAIJU") for brevity.


Before you begin
------------------------------------------------

Setting environment variables
************************************************

*Source* (not *run*) the environment setup scripts for the ``kaiju`` software,
add paths to ``TIEGCMHOME`` and ``TIEGCMDATA``, then activate the ``conda``
environment that holds ``kaipy`` (see :doc:`/python/index`). For example:

.. code-block:: bash

    source /path/to/your/kaiju-clone/scripts/setupEnvironment.sh
    export TIEGCMHOME=/path/to/your/tiegcm
    export TIEGCMDATA=/path/to/your/tiegcm/data
    module load conda
    conda activate kaiju-3.12

.. note::

    The ``TIEGCMHOME`` and ``TIEGCMDATA`` environment variables are required
    for running the GRT model. They should point to the TIEGCM source code
    directory and the TIEGCM data directory, respectively.

    The TIEGCMDATA directory is located in the following locations:
        - On ``derecho``: ``/glade/campaign/hao/itmodel/tiegcm3.0/new_data``
        - On ``aitken``: ``/nobackup/nrao3/tiegcm/tiegcm3.0/data``
        - The required data files can be downloaded from the NCAR Globus endpoint using the following link: `TIEGCM Data Files <https://app.globus.org/file-manager?origin_id=b2502c58-c3eb-470f-86d4-cbdcd0aeb6c8&origin_path=%2F>`_

The ``kaiju`` software (in the ``build_grt`` subdirectory) and the TIEGCM
executables must be built before you continue. Instructions for both builds
are available :doc:`here </building/buildDerecho_GRT>`.

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
information on running GRT with ``engage.py`` is available
:doc:`here </makeitso/engage_grt>`.

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

    mkdir quickstart_grt
    cd quickstart_grt

2. Copy the TIEGCM executables you built in the previous steps to your run
   directory. The ``Standalone Executable`` and ``Coupled Executable`` prompts
   default to these copies; you can also type the paths to the executables at
   those prompts instead.

.. note::
    The TIEGCM executables for a ``Q`` run were built in the ``tiegcm_build_Q`` directory in
    the build guide.

.. code-block:: bash

    cp $TIEGCMHOME/tiegcm_build_Q/stdout/tiegcm.exe .
    cp $TIEGCMHOME/tiegcm_build_Q/stdout/tiegcm.x .


3. Run ``engage.py``. GRT is the default run type of ``engage.py``, so no
   ``--component`` option is needed.

.. code-block:: bash

    engage.py

This opens an arrow-key wizard that asks the questions below, one at a time: type a
value or press :kbd:`Enter` to accept the default and go to the next one.
The keys are listed at the bottom of the screen: ``↑ back   ↓/Enter next   Ctrl-S finish   Ctrl-R reset   Ctrl-C cancel``.
This is the recommended way to run ``engage.py``.
The wizard asks the following questions. They are shown here, for reference,
as the plain prompts that ``--no-tui`` prints (the fallback when your terminal
cannot show the wizard; the tool also falls back by itself):

.. code-block:: bash

    Name to use for PBS job(s) [geospace]:
    Start date for simulation (yyyy-mm-ddThh:mm:ss) [2001-06-01T23:00:00]:
    Stop date for simulation (yyyy-mm-ddThh:mm:ss) [2001-06-02T01:00:00]:
    GAMERA grid type (D|Q|O) [Q]:
    Name of HPC system (derecho|aitken|linux|perlmutter) [derecho]:
    PBS account name (your projects: <YOUR_ACCOUNT_HERE>): <YOUR_ACCOUNT_HERE>
    Run directory [.]:
    Path to kaiju installation [<YOUR_KAIJUHOME_HERE>]:
    Path to kaiju build directory [<YOUR_KAIJUHOME_HERE>/build_grt]:
    PBS queue name (develop|main) [main]:
    Job priority (premium|regular|economy) [economy]:
    Requested wall time for each PBS job segment (HH:MM:SS) [12:00:00]:
    Do you want to split your job into multiple segments? (Y|y|N|n) [Y]:
    Suggested segment for walltime 12:00:00 on main: 7200 s (2 h); valid values (seconds): [3600 (1 h), 7200 (2 h)]
    Segment length in simulated seconds [7200]: 3600
    Conda environment to use for the simulation [<YOUR_CONDA_ENVIRONMENT_HERE>]:
    Extend TFIN by dtCouple - 1 seconds (T|F) [T]:

.. warning::

    Make sure to set ``Path to kaiju build directory`` to the directory where you built the
    ``voltron_mpi.x`` executable with the module set for GRT runs. If you followed the build
    instructions, this should be the ``build_grt`` subdirectory. This is required for running
    the model in GRT mode.

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
range to use in the solar wind boundary condition file.

You should see output similar to this:

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
    05-Oct-26 02:12:50: /path/to/your/conda/envs/kaiju-3.12/lib/python3.12/site-packages/spacepy/time.py:2448: UserWarning: Leapseconds may be out of date. Use spacepy.toolbox.update(leapsecs=True)
      _read_leaps()

    Retrieving f10.7 data from CDAWeb
    Retrieving solar wind data from CDAWeb
            Using Bx fields
    Bx Fit Coefficients are  [ 1.14301522 -0.42752137 -0.36840806]
    Saving "OMNI_HRO_1MIN.txt_bxFit.png"
    Converting to Gamera solar wind file
            Found 21 variables and 360 lines
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


    ..
    ..
    Template creation complete!


    The PBS scripts ['./gr_spinup/geospace-SPINUP.pbs', './gr_spinup/geospace-WARMUP-01.pbs', './gr_spinup/geospace-WARMUP-02.pbs', './gr_spinup/geospace-WARMUP-03.pbs', './gr_spinup/geospace-WARMUP-04.pbs', './geospace-01.pbs', './geospace-02.pbs'] have been created, each with a corresponding XML file. To submit the jobs with the proper dependency (to ensure each segment runs in order), please run the script geospace_pbs.sh like this:
    bash geospace_pbs.sh

``engage.py`` will then prompt you for the following additional information from ``tiegcmrun``:

.. code-block:: bash
    :force:

    Instructions:
    -> Default Selected input parameter is given in GREEN
    -> Warnings and Information are given in YELLOW
    -> Errors are given in RED
    -> Valid values (if any) are given in brackets eg. (value1 | value2 | value3)
    -> Enter '?' for any input parameter to get a detailed description


    Run Options:
    User Mode = BASIC
    Compile = False
    Execute = False
    Coupling = True
    Engage = True


    Directory of model [<YOUR_TIEGCMHOME_HERE>]:
    Directory of Tiegcm Data Files [/glade/campaign/hao/itmodel/tiegcm3.0/new_data]:
    Standalone Executable [<YOUR_RUN_DIRECTORY_HERE>/tiegcm.exe]:
    Coupled Executable [<YOUR_RUN_DIRECTORY_HERE>/tiegcm.x]:
    Low = 70, Medium = 140 , High = 200
    F107 flux level for TIEGCM spin up (low|medium|high) [low]:
    SOURCE file location [/glade/campaign/hao/itmodel/tiegcm3.0/new_data/source/junsol_f70.nc]:
    STEP number [30]:
    NSTEP_SUB number [10]:
    Secondary Output Fields [['TN', 'UN', 'VN', 'O2', 'O1', 'HE', 'N2', 'N2D', 'N4S', 'NO', 'NE', 'Z', 'WN', 'DEN', 'TEC']] / ENTER to go next:
    High-latitude potential model that is going to be used (HEELIS|WEIMER) [HEELIS]:
    If GPI_NCFILE is specified, then KP and POWER/CTPOTEN are skipped. If further POTENTIAL_MODEL is WEIMER and IMF_NCFILE is specified, then the Weimer model and aurora will be driven by the IMF data, and only F107 and F107A will be read from the GPI data file.
    Enter 'gen' to generate a GPI file for this run's dates with gcmprocpy (27-day trailing F10.7 average).
    GPI file [/glade/campaign/hao/itmodel/tiegcm3.0/new_data/boundary_files/GPI/gpi_1960001-2026180.nc]:

After these inputs, the script interpolates the source file for TIEGCM, and
generates the namelist and PBS files for the run.

You should see output similar to this:

.. code-block:: bash

    ..
    ..
    Interpolating primary file /glade/campaign/hao/itmodel/tiegcm3.0/new_data/source/junsol_f70.nc to create new primary file <YOUR_RUN_DIRECTORY_HERE>/tiegcm_standalone/geospace-tiegcm-standalone_prim.nc at horizontal resolution 2.5 and vertical resolution 0.25 with zitop 7.0.
    Creating new primary file:  <YOUR_RUN_DIRECTORY_HERE>/tiegcm_standalone/geospace-tiegcm-standalone_prim.nc
    Creating 2 PBS scripts
    GR_pbs_scripts = ['./gr_spinup/geospace-SPINUP.pbs', './gr_spinup/geospace-WARMUP-01.pbs', './gr_spinup/geospace-WARMUP-02.pbs', './gr_spinup/geospace-WARMUP-03.pbs', './gr_spinup/geospace-WARMUP-04.pbs']
    Tiegcm_pbs_scripts = ['tiegcm_standalone/geospace-tiegcm-standalone-01.pbs', 'tiegcm_standalone/geospace-tiegcm-standalone-02.pbs']
    GRT_pbs_scripts = ['./geospace-01.pbs', './geospace-02.pbs']
    submit_all_jobs_script = geospace_pbs.sh

Looking at files generated by ``engage.py``
**************************************************************

You should now see the following files in your run directory:

.. code-block:: bash

    ls
    bcwind.h5               geospace-02.xml              raijuconfig.h5
    engage_parameters.json  geospace_pbs.sh              tiegcm.exe
    geospace-01.inp         gr_spinup                    tiegcmrun_parameters.json
    geospace-01.pbs         lfmQ.h5                      tiegcm_standalone
    geospace-01.xml         makeitso_parameters.json     tiegcm.x
    geospace-02.inp         OMNI_HRO_1MIN.txt_bxFit.png  voltron_mpi.x
    geospace-02.pbs         OMNI_HRO_1MIN.txt.png

The GAMERA spin-up and warm-up segments are in the ``gr_spinup``
subdirectory:

.. code-block:: bash

    ls gr_spinup
    bcwind.h5               geospace-WARMUP-02.pbs  geospace-WARMUP-04.xml
    geospace-SPINUP.pbs     geospace-WARMUP-02.xml  lfmQ.h5
    geospace-SPINUP.xml     geospace-WARMUP-03.pbs  raijuconfig.h5
    geospace-WARMUP-01.pbs  geospace-WARMUP-03.xml  voltron_mpi.x
    geospace-WARMUP-01.xml  geospace-WARMUP-04.pbs

The TIEGCM spin-up segments are in the ``tiegcm_standalone`` subdirectory:

.. code-block:: bash

    ls tiegcm_standalone
    geospace-tiegcm-standalone-01.inp  geospace-tiegcm-standalone-02.pbs
    geospace-tiegcm-standalone-01.pbs  geospace-tiegcm-standalone_prim.nc
    geospace-tiegcm-standalone-02.inp  tiegcm_resolution_upscale.py

There are several types of files created for each of the jobs, including:

* ``*.pbs``
    These are the PBS scripts that will be submitted to the job scheduler to run
    the segments of the simulation.
* ``*.xml``
    These are the XML files that contain the parameters for GAMERA and RAIJU of the
    segment.
* ``*.inp``
    These are the namelist files that contain parameters for TIEGCM of the segment.
* ``*_parameters.json``
    These are the JSON files that contain the parameters for the simulation. They
    are generated by the ``engage.py`` script with all the parameters required to run the
    simulation.

The run is divided into segments:

* ``gr_spinup/geospace-SPINUP.*``
    This segment runs the GAMERA model to create the initial conditions for the
    simulation. It is run first, and its output is used by the next segment.
* ``gr_spinup/geospace-WARMUP-**.*``
    These segments run the GAMERA RAIJU model to "warm up" for the coupled model execution.
    The ``-01``, ``-02``, etc. suffixes indicate the segment number, and the
    segments are run in order.
* ``tiegcm_standalone/geospace-tiegcm-standalone-**.*``
    These segments run the TIEGCM model to create the initial conditions for the coupled model.
    The ``-01``, ``-02``, etc. suffixes indicate the segment number, and the
    segments are run in order.
* ``geospace-**.*``
    These segments run the GRT coupled model. The ``-01``, ``-02``, etc.
    suffixes indicate the segment number, and the segments are run
    in order.

This image shows how the segments are run in order:

.. image:: GRTSegment.png

The image files are summaries of the CDAWeb data used in the initial condition
file (``bcwind.h5``). Those plots should look similar to this:

.. image:: Bx_fit_grt.png

.. image:: sw_grt.png

Submitting the GRT model run
*************************************

Finally, submit the model run using the script generated by ``engage.py``.
You will see the resulting PBS job ID for each segment (your job IDs will
differ from what is shown below). The job IDs are also saved in the file
``geospace.jobids``.

.. code-block:: bash

    bash geospace_pbs.sh
    geospace-SPINUP.pbs=9770226.desched1
    geospace-WARMUP-01.pbs=9770227.desched1
    geospace-WARMUP-02.pbs=9770228.desched1
    geospace-WARMUP-03.pbs=9770229.desched1
    geospace-WARMUP-04.pbs=9770230.desched1
    geospace-tiegcm-standalone-01.pbs=9770231.desched1
    geospace-tiegcm-standalone-02.pbs=9770232.desched1
    geospace-01.pbs=9770233.desched1
    geospace-02.pbs=9770234.desched1

Replaying the run
*************************************

The options you entered are saved in ``engage_parameters.json``,
``makeitso_parameters.json`` and ``tiegcmrun_parameters.json``. To generate the
same run again, run this command in the run directory:

.. code-block:: bash

    engage.py --clobber -eo engage_parameters.json -mo makeitso_parameters.json -to tiegcmrun_parameters.json

Understanding the output files of the GRT model run
**************************************************************************

Once the job is started in the queue, it should take about 80 minutes to run
(on ``derecho``). When complete, you will see many new HDF5 files in your
run directory, along with PBS housekeeping files and logs. The most important
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

* ``geospace_sech_*.nc``

    These are secondary output files that contain the results from the
    `TIEGCM <https://tiegcm-docs.readthedocs.io/>`_ portion of the
    `MAGE <https://cgs.jhuapl.edu/Models>`_ model.

* ``geospace_*.gam.Res.RRRRR.h5``

    These are checkpoint files generated during the simulation which can be
    used as restart points for future simulations.

* ``geospace_prim_*.nc``

    These are the primary output files from the TIEGCM portion of the model
    that are designed as checkpoint files.

* ``geospace_temp_*.nc``

    These are temporary output files from the TIEGCM portion of the model.

Visualizing the results
-----------------------

Now perform a quick visualization of the results from your model using the
``msphpic.py`` script, provided in the ``kaipy`` package.

.. code-block:: bash

    msphpic.py -id geospace

This script will create a file called ``qkmsphpic.png``, which should look
similar to this:

.. image:: qkmsphpic_grt.png
