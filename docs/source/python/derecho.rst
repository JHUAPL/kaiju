Creating a Python environment for ``kaiju`` on ``derecho``
==========================================================


Introduction
------------

This file describes how to set up a Python ``conda`` environment on
``derecho`` which can be used to run the ``kaiju`` code.

These steps assume that ``pip`` (only) is used for installing additional
packages, and packages are only taken from PyPI.

These instructions assume you are using the ``bash`` shell.



Building the python environment
-------------------------------

To create a Python environment for ``kaiju``, we first install/load the 
Conda software package, then create the environment, then populate it with 
the required Python packages.

.. note:: 
    These instructions are designed for novice users. If you are
    comfortable with building and managing ``conda``-based environments, feel
    free to look into the CISL-provided
    `conda documentation <https://ncar-hpc-docs.readthedocs.io/en/stable/environment-and-software/user-environment/conda>`_
    for futher details.


.. code-block:: bash


    # Load the conda module.
    module load conda

    # Now create the environment for kaiju, specifying only the python
    # version.
    conda create -n kaiju-3.12 python=3.12

    # Activate the new environment.
    conda activate kaiju-3.12

    # Install the kaipy software.
    # If you are using kaipy via pip:
    pip install kaipy
    # OR
    # If you are using a clone of the kaipy repository:
    cd /path/to/kaipy
    pip install -e .
    # -e is used to install the package in editable mode, which allows
    # you to make changes to the code and have them reflected without needing
    # to reinstall the package.

Using the python environment
----------------------------

Once your python environment is created, you must activate it for use with the
``kaiju`` software:

.. code-block:: bash

    conda activate kaiju-3.12