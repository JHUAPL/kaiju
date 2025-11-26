Building a python environment for ``kaiju`` on Ubuntu 20.04
===============================================================

Introduction
------------

This file describes how to set up a python 3.8-based ``conda`` environment on an Intel-based machine running Ubuntu 20.04 which can be used to run the ``kaiju`` code.

These steps assume that ``pip`` (only) is used for installing additional packages, and packages are only taken from PyPI (and the default ``conda`` repository, if needed) - no other repositories are used.

*A NOTE OF EXPLANATION:* These instructions install ``miniconda3`` into a Ubuntu 20.04-specific subdirectory of the home directory, to maintain compatibility with instructions for other systems. Feel free to install ``miniconda3`` wherever is convenient.

Building the python environment
-------------------------------

.. code-block:: shell

   # Download the installer.
   cd $HOME
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

   # Run the installer.
   # Install into $HOME/ubuntu-20.04/miniconda3, use all defaults.
   bash ./Miniconda3-latest-MacOSX-x86_64.sh

   # NOTE: This installation *should* run "conda init", which *should*
   # add the conda setup code to ~/.bashrc or ~/.bash_profile.

   # NOTE: This installation creates ~/.conda/environments.txt, but
   # nothing else. This file just contains the path to the miniconda3
   # installation directory.

   # Turn off auto-activation of base environment.
   # NOTE: This command creates ~/.condarc.
   conda config --set auto_activate_base false

   # Update everything to latest version.
   conda update --all

   # NOTE: This creates the directory ~/.cache/conda.

   # Now create the environment for kaiju, specifying only the python version.
   conda create -n kaiju-3.8 python=3.8

   # NOTE: This adds the path to the new environment in
   # ~/.conda/environments.txt.

   # Activate the new environment.
   conda activate kaiju-3.8

   # Update everything to latest version.
   conda update --all

   # IMPORTANT: Use a pip cache specific to this machine.
   export PIP_CACHE_DIR=$HOME/ubuntu-20.04/pip_cache

   # Make sure CDF, HDF5, and GEOS have been built and installed.

   # Install packages required by kaipy.
   pip install -r requirements.txt

   # Save the environment description.
   conda list >& kaiju-3.8-ubuntu-20.04-YYYYMMDD.env



Where ``requirements.txt`` is a text file containing:

.. code-block::

   ai.cdas
   alive_progress
   cartopy
   cdasws
   cdflib
   configparser
   dataclasses_json
   h5py
   jinja2
   matplotlib
   netCDF4
   progressbar
   pyhdf
   pytest
   scipy
   spacepy
   sunpy
   xarray


Using the python environment
----------------------------

Once your python environment is created, you must activate it for use with the ``kaiju`` software:

.. code-block:: shell

   conda activate kaiju-3.8   # or whatever name you used for the python environment
