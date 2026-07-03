Building the ``kaiju`` software on CentOS-Stream 9
==================================================

Introduction
------------

This page provides instructions for building the ``kaiju`` software on an Intel-based CentOS-Stream 9 system. These instructions were developed for using a CentOS-Stream 9 virtual machine running under VirtualBox on a Mac, but the instructions should work for other Linux distributions, although minor changes may be needed, especially if your Linux distribution does no use the ``yum`` package manager.

NOTE: These instructions should also work, essentially unchanged, on a RHEL (Red Hat Enterprise Linux) system.

Step 0: Prepare your system
---------------------------

A few tools must be installed on your system before you begin building the ``kaiju`` software. These instructions assume you have permission to use ``sudo`` to install software.

.. code-block::

   #!shell
   sudo yum install epel-release
   sudo yum install make
   sudo yum install gcc
   sudo yum install cmake
   sudo yum install git
   sudo yum install git-lfs
   sudo yum install geos-devel


Install the Intel compiler suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``kaiju`` software is developed using the Intel compiler suite, available `here <https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.29x43u>`_. The compiler suite is free, but comes in several optional parts. You will need to install the `Intel oneAPI Base Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`_ and the `Intel oneAPI HPC Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`_. Once installed, this software should be available on your system under ``/opt/intel``.

Once the Intel tools are installed, you will need to make them available during the build process. Set up the tools as follows:

.. code-block::

   #!shell
   source /opt/intel/oneapi/setvars.sh


Step 1: Build prerequisite libraries
------------------------------------

NASA CDF library
^^^^^^^^^^^^^^^^

The `NASA CDF (Common Data Format) library <https://cdf.gsfc.nasa.gov/>`_ is used in parts of the ``kaiju`` post-processing software when fetching spacecraft data from `CDAWeb <https://cdaweb.gsfc.nasa.gov/>`_. Prior to building the ``kaiju`` software, the CDF library must be built and installed on your system, since it is not available via the ``yum`` package manager.

Instructions for building and installing the CDF library on your Mac are available `here <centos-stream-9_build_cdf>`__.

HDF5 library
^^^^^^^^^^^^

The `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ file format is used extensively for data storage in the ``kaiju`` software. The HDF5 software is available via ``yum``\ , but this version will not link with the ``kaiju`` software when built with the Intel compiler suite. Prior to building the ``kaiju`` software, the HDF5 library must be built and installed on your system.

Instructions for building and installing the HDF5 library on your Mac are available `here <centos-stream-9_build_hdf5>`__.

NOTE: If you successfully build the ``kaiju`` software using the ``yum``\ -installed version of HDF5 on CentOS-Stream 9, please let us know how you did it!

Step 2: Create a python environment
-----------------------------------

Most of the ``kaiju`` software for pre-processing, post-processing, and analysis is written in `Python <https://www.python.org/>`_. Python is available in many forms (or 'distributions'), but we recommend use of the `Miniconda] distribution <https://docs.conda.io/en/latest/miniconda.html>`_ for simplicity and compactness.

Instructions for installing python and building a python environment on your Mac are available `here <centos-stream-9_build_python>`__.

Step 3: Compile the ``kaiju`` software
------------------------------------------

The ``kaiju`` software can be built in serial or MPI versions. The serial version should be used when running the code on a single computer. The MPI version should be used when running on an HPC system (a supercomputer).

Instructions for building the serial version of the ``kaiju`` software on CentOS-Stream 9 are available `here <centos-stream-9_build_kaiju>`__. The MPI version of ``kaiju`` is not supported on CentOS-Stream 9 at this time.

.. toctree::
    :maxdepth: 1

    centos-stream-9_build_cdf
    centos-stream-9_build_hdf5
    centos-stream-9_build_kaiju
    centos-stream-9_build_python
