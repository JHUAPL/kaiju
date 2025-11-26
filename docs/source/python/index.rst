Creating a Python environment for ``kaiju``
===========================================

``kaiju`` software is written in modern Fortran but we use Python for
pre- and post-processing. To set up a run, you will need to use a
Python script, provided with ``kaiju``, that is called :doc:`makeitso
</makeitso/index>` or its variants called ``engage`` and
``makeitso-helio``. For some of its functions (e.g., grid
generation or downloading NASA data to set boundary conditions),
``makeitso`` relies on our Python package `kaipy
<https://kaipy-docs.readthedocs.io/en/latest/index.html>`_. You will also need to install ``kaipy`` to be able to analyze
the simulation results.

For this reason, before you can run ``kaiju`` software (e.g., ``MAGE``
or ``GAMERA-helio``), you will need to install ``kaipy`` first. These
pages describe how to set up and configure a Python environment with
the proper prerequisites to use the ``kaiju`` and ``kaipy`` tools.


.. toctree::
    :maxdepth: 1

    derecho
    aitken
