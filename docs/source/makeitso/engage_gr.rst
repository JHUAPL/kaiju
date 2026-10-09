Engage for MAGE without TIEGCM (GR)
=============================================


Introduction
------------

The Python script ``engage.py`` can be used to configure and run GR MAGE
(that is, the geospace application of the ``kaiju`` software without
TIEGCM). With the ``--component gr`` option, ``engage.py`` starts the
``makeitso.py`` script, which asks for the parameters needed for the model
run.

The questions, their order and their default values are the same as when you
run ``makeitso.py`` directly, and the same files are created. More information
on ``engage.py`` is available :doc:`here </makeitso/engage>`.


Running the ``engage.py`` script
----------------------------------

To configure a GR run, run ``engage.py`` with ``--component gr`` in your run
directory:

.. code-block:: bash

  engage.py --component gr

This opens an arrow-key wizard that asks the questions below, one at a time: type a
value or press :kbd:`Enter` to accept the default and go to the next one.
The keys are listed at the bottom of the screen: ``↑ back   ↓/Enter next   Ctrl-S finish   Ctrl-R reset   Ctrl-C cancel``.
This is the recommended way to run ``engage.py``.
The wizard asks the following questions. They are shown here, for reference,
as the plain prompts that ``--no-tui`` prints (the fallback when your terminal
cannot show the wizard; the tool also falls back by itself):

.. code-block:: bash

  engage: launching standalone makeitso (GR: GAMERA+RAIJU) ...


  Instructions:
  -> Default selected input parameter is given in GREEN
  -> Warnings and information are given in YELLOW
  -> Errors are given in RED
  -> Valid values (if any) are given in brackets eg. (value1 | value2 | value3)


  Run Options:
  Name to use for PBS job(s) [geospace]:

Each prompt is explained in the section "An example in ``BASIC`` mode" of the
:doc:`makeitso </makeitso/makeitso>` documentation. The options of
``engage.py`` are listed :doc:`here </makeitso/engage>`.

When finished, the script creates the input files for the run, the PBS job
scripts with their XML files, and the script ``geospace_pbs.sh`` that submits
the PBS jobs (here ``geospace`` is the name of the PBS job). A segmented run
also creates the subdirectory ``gr_spinup``.


Using JSON files for ``engage.py``
----------------------------------

When finished, the script creates the file ``makeitso_parameters.json`` in the
run directory. This file contains a record of all of the parameters used in
your simulation. This file can be passed back to ``engage.py`` with the
``-mo`` option in a subsequent session to repeat the simulation, and also
provides a convenient starting point for minor tweaks to your simulation
parameters. In the run directory, add ``--clobber`` so that the file can be
overwritten:

.. code-block:: bash

  engage.py --component gr --clobber -mo makeitso_parameters.json

The closing message of the session suggests ``makeitso.py``; the command
``makeitso.py --clobber -o makeitso_parameters.json`` works as well.

To start a new run from a JSON file, download a template:

  - derecho: :download:`makeitso_input.json <engage_template/derecho/gr/makeitso_input.json>`
  - aitken: :download:`makeitso_input.json <engage_template/aitken/gr/makeitso_input.json>`

Copy it into your run directory, set your PBS account (not needed on
``aitken``), the start and stop dates and the GAMERA grid type, and run:

.. code-block:: bash

  engage.py --component gr -mo makeitso_input.json

More information on JSON files is available
:doc:`here </makeitso/engage>`.

A complete GR run, from the run directory to the results, is described
:doc:`here </running/geoGRQuickStart>`.
