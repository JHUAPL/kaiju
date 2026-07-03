Building the ``kaiju`` software on Ubuntu 20.04
===============================================

Introduction
------------

This page provides instructions for building the ``kaiju`` software on an
Intel-based Ubuntu 20.04 system. These instructions were developed for using a
Ubuntu 20.04 virtual machine running under VirtualBox on a Mac, but the
instructions should work for other Linux distributions, although minor changes
may be needed, especially if your Linux distribution does no use the ``apt``
package manager.

Step 0: Prepare your system
---------------------------

A few tools must be installed on your system before you begin building the
``kaiju`` software. These instructions assume you have permission to use
``sudo`` to install software.

.. code-block:: shell

   sudo apt-get install build-essential
   sudo apt-get install cmake
   sudo apt-get install git
   sudo apt-get install git-lfs
   sudo apt-get install libgeos-dev
   sudo apt-get install pkg-config


Install the Intel compiler suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``kaiju`` software is developed using the Intel compiler suite, available
`here <https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.29x43u>`_.
The compiler suite is free, but comes in several optional parts. You will need
to install the `Intel oneAPI Base Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`_
and the `Intel oneAPI HPC Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`_.
Once installed, this software should be available on your system under
``/opt/intel``.

Once the Intel tools are installed, you will need to make them available
during the build process. Set up the tools as follows:

.. code-block:: shell

   source /opt/intel/oneapi/setvars.sh

Step 1: Build prerequisite libraries
------------------------------------

NASA CDF library
~~~~~~~~~~~~~~~~

The `NASA CDF (Common Data Format) library <https://cdf.gsfc.nasa.gov/>`_ is
used in parts of the ``kaiju`` post-processing software when fetching
spacecraft data from `CDAWeb <https://cdaweb.gsfc.nasa.gov/>`_. Prior to
building the ``kaiju`` software, the CDF library must be built and installed
on your system, since it is not available via the ``apt`` package manager.

Instructions for building and installing the CDF library on your system are
available :doc:`here <ubuntu-20.04_build_cdf>`.

HDF5 library
~~~~~~~~~~~~

The `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ file format is used
extensively for data storage in the ``kaiju`` software. The HDF5 software is
available via ``apt``, but this version will not link with the ``kaiju``
software when built with the Intel compiler suite. Prior to building the
``kaiju`` software, the HDF5 library must be built and installed on your
system.

Instructions for building and installing the HDF5 library on your system are
available :doc:`here <ubuntu-20.04_build_hdf5>`.

NOTE: If you successfully build the ``kaiju`` software using the
``apt``-installed version of HDF5 on Ubuntu 20.04, please let us know how you
did it!

Step 2: Create a python environment
-----------------------------------

Most of the ``kaiju`` software for pre-processing, post-processing, and
analysis is written in `Python <https://www.python.org/>`_. Python is
available in many forms (or 'distributions'), but we recommend use of the
`Miniconda distribution <https://docs.conda.io/en/latest/miniconda.html>`_
for simplicity and compactness.

Instructions for installing python and building a python environment on your
system are available :doc:`here <ubuntu-20.04_build_python>`.

Step 3: Compile the ``kaiju`` software
--------------------------------------

The ``kaiju`` software can be built in serial or MPI versions. The serial
version should be used when running the code on a single computer. The MPI
version should be used when running on an HPC system (a supercomputer).

Instructions for building the serial version of the ``kaiju`` software on
Ubuntu 20.04 are available :doc:`here <ubuntu-20.04_build_kaiju>`. The MPI
version of ``kaiju`` is not supported on Ubuntu 20.04 at this time.

.. toctree::

   ubuntu-20.04_build_cdf
   ubuntu-20.04_build_hdf5
   ubuntu-20.04_build_kaiju
   ubuntu-20.04_build_python