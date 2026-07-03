Running the ``kaiju`` software
==============================


Before you begin
----------------

These instructions use the ``makeitso.py``, ``engage.py`` and
``makeitso-gamhelio.py`` tools, all belonging to the family of scripts
we call ``makeitso``,
to simplify the process of configuring and running the ``kaiju`` code. These
scripts are provided as part of the ``kaiju`` package, but also require the
``kaipy`` Python package, which can be installed into your Python environment
with the command:

.. code-block:: bash
    
    pip install kaipy

The documentation for ``kaipy``, including instructions for building a Python
environment with all prerequisites, is available on
`Read the Docs <https://kaipy-docs.readthedocs.io/en/latest/index.html>`_.

.. important:: **Building the NASA CDF (Common Data Format) library.**
    The ``kaipy`` package requires the NASA CDF library, which
    should be available when you create your Python environment for ``kaipy``.
    This requirement will be removed in a future release of the ``kaiju``
    software. Instructions for installing the CDF library can be found
    :doc:`here <buildCDF>`.
         
.. important::

    The ``kaipy`` code requires python version 3.10 or later. If you try to
    run the ``kaipy`` code using a different version of ``python``, you may
    encounter compatibility issues.


Running Guides
-----------------

.. toctree::
    :maxdepth: 1

    geoGRQuickStart
    geoGTRQuickStart
    helioQuickStart
