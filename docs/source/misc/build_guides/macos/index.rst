Building the ``kaiju`` software on MacOS
========================================

Introduction
------------
This page provides instructions for building the ``kaiju`` software on an
Intel-based MacOS system. These instructions were developed for MacOS Ventura,
but the instructions should work for other MacOS versions.

Step 0: Prepare your Mac
------------------------

A few tools must be installed on your Mac before you begin building the
``kaiju`` software.

Install Xcode
~~~~~~~~~~~~~

`Xcode <https://developer.apple.com/xcode/>`_ is the standard software
development suite for MacOS systems. You can download and install Xcode from
the MacOS App Store on your Mac.

Install HomeBrew
~~~~~~~~~~~~~~~~

`HomeBrew <https://brew.sh/>`_ is a software package management system for
MacOS. HomeBrew will be used to install tools required for building the
``kaiju`` software.

There are other package management systems available for MacOS. If you choose
to use an alternate package management system, replace these instructions with
the equivalent instructions for your chosen system.

Install ``cmake``
~~~~~~~~~~~~~~~~~

`CMake <https://cmake.org/>`_ is used to control the build process for the
``kaiju`` software. This tool is not provided on MacOS systems. You can
install it with HomeBrew with the command:

.. code-block:: shell

   brew install cmake

This command will install ``cmake`` and its support files under ``/usr/local``
on your Mac.

Install the Intel compiler suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``kaiju`` software is developed using the Intel compiler suite, available
`here <https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.29x43u>`_.
The compiler suite is free, but comes in several optional parts. You will need
to install the `Intel oneAPI Base Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`_
and the `Intel oneAPI HPC Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`_.
Once installed, this software should be available on your Mac under
``/opt/intel``.

Once the Intel tools are installed, you will need to make them available
during the build process. Add the Intel tools to your command path as follows:

.. code-block:: shell

   export INTEL_HOME=/opt/intel   # Or wherever you installed it
   export PATH=$INTEL_HOME/oneapi/compiler/latest/mac/bin/intel64:$PATH

Step 1: Build prerequisite libraries
------------------------------------

NASA CDF library
~~~~~~~~~~~~~~~~

The `NASA CDF (Common Data Format) library <https://cdf.gsfc.nasa.gov/>`_ is
used in parts of the ``kaiju`` post-processing software when fetching
spacecraft data from `CDAWeb <https://cdaweb.gsfc.nasa.gov/>`_. Prior to
building the ``kaiju`` software, the CDF library must be built and installed
on your Mac, since it is not available via the HomeBrew package manager.

Instructions for building and installing the CDF library on your Mac are
available :doc:`here <macos_build_cdf>`.

HDF5 library
~~~~~~~~~~~~

The `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ file format is used
extensively for data storage in the ``kaiju`` software. The HDF5 software is
available via HomeBrew, but this version will not link with the ``kaiju``
software when built with the Intel compiler suite. Prior to building the
``kaiju`` software, the HDF5 library must be built and installed on your Mac.

Instructions for building and installing the HDF5 library on your Mac are
available :doc:`here <macos_build_hdf5>`.

NOTE: If you successfully build the ``kaiju`` software using the
HomeBrew-installed version of HDF5 on MacOS, please let us know how you did
it!

GEOS library
~~~~~~~~~~~~

The `GEOS <https://libgeos.org/>`_ library is needed for certain
geography-related coordinate transformations. It can be installed with the
command:

.. code-block:: shell

   brew install geos


Step 2: Create a python environment
-----------------------------------

Most of the ``kaiju`` software for pre-processing, post-processing, and
analysis is written in `Python <https://www.python.org/>`_. Python is
available in many forms (or 'distributions'), but we recommend use of the
`Miniconda distribution <https://docs.conda.io/en/latest/miniconda.html>`_
for simplicity and compactness.

Instructions for installing python and building a python environment on your
Mac are available :doc:`here <macos_build_python>`.

Step 3: Compile the ``kaiju`` software
------------------------------------------

The ``kaiju`` software can be built in serial or MPI versions. The serial
version should be used when running the code on a single computer, such as a laptop. The MPI version should be used when running on an HPC system
(a supercomputer).

Instructions for building the serial version of the ``kaiju`` software on
MacOS are available :doc:`here <macos_build_kaiju>`. The MPI version of
``kaiju`` is not supported on MacOS at this time.

.. toctree::
   :maxdepth: 1

   macos_build_cdf
   macos_build_hdf5
   macos_build_kaiju
   macos_build_python