Building the ``kaiju`` software on ``derecho`` for MAGE - With WACCM-X (GRX)
==============================================================================================


Introduction
------------

This page provides instructions for building the ``kaiju`` software on the
``derecho`` supercomputer. These instructions assume that you have cloned the
``kaiju`` repository.


Prepare your software environment
---------------------------------

Like most HPC systems, ``derecho`` uses the ``module`` system to manage the
versions of software packages available to the user. When you log in to
``derecho``, the following modules are loaded by default:

.. code-block:: bash

   module list

   Currently Loaded Modules:
     1) ncarenv/23.09  (S)   4) ncarcompilers/1.0.0   7) netcdf/4.9.2
     2) craype/2.7.23        5) cray-mpich/8.1.27
     3) intel/2023.2.1       6) hdf5/1.12.2

     Where:
      S:  Module is Sticky, requires --force to unload or purge

This set of modules **cannot** be used to build the ``kaiju`` code.

Start by purging any currently-loaded modules, then loading the following
module set for MAGE runs coupled with WACCM-X (known as "GRX"):

.. code-block:: bash

    module purge
    module load cesmdev/1.0
    module load ncarenv/23.09
    module load craype
    module load cmake
    module load intel-classic/2023.2.1
    module load mkl
    module load ncarcompilers
    module load cray-mpich/8.1.27
    module load hdf5-mpi/1.12.2
    module load netcdf-mpi/4.9.2
    module load parallel-netcdf/1.12.3
    module load parallelio/2.6.2
    module load esmf/8.6.0
    module load nco

.. important::

    You must use these exact versions of the modules to ensure the software
    compiles properly. If you use different versions of any of these modules,
    a successful build cannot be guaranteed. This module list is current as of
    **4 October 2026**, and is subject to change as the compute environment
    changes.

Build the ``kaiju`` software
----------------------------

These instructions show how to build the MPI version of the ``kaiju``
software. The GRX version is built in the subdirectory ``build_grx``
under the ``kaiju`` source code directory. In practice, you can place the
build directory in any convenient location.

.. code-block:: bash

    # Move to your kaiju clone.
    cd /path/to/kaiju

    # Create the GRX build directory and enter it.
    mkdir build_grx
    cd build_grx

    # Run cmake to create the Makefile, saving output.
    # NOTE: The FC definition is *required* for proper cmake operation.
    FC=`which ifort` cmake -DENABLE_MPI=ON -DALLOW_INVALID_COMPILERS=ON .. >& cmake.out

    # Compile the MAGE model for geospace simulations
    make -j8 voltron_mpi.x >& make-voltron.out

When finished, your build directory will contain a ``bin``
subdirectory which will contain the compiled ``kaiju`` executables.

Get the WACCM-X source code
---------------------------

The GRX model runs WACCM-X from a clone of the CAM source code. You
only need to clone it and tell ``engage.py`` where it is.

Getting the WACCM-X source code
************************************************

The CAM source code can be obtained by cloning the ``CAM`` repository on
GitHub. Clone the branch and commit used for GRX runs, apply the CMEPS patch
provided in the ``kaiju`` repository, then fetch the CAM components:

.. code-block:: bash

    # Clone the CAM source code and enter it.
    git clone --branch waccmx_mage_cplng https://github.com/fvitt/CAM.git /path/to/your/cam
    cd /path/to/your/cam

    # Check out the commit used for GRX runs.
    git checkout da3bc8d52f652b64b12d9978ba97b50fbcb750c7

    # Apply the CMEPS patch from your kaiju clone.
    git apply /path/to/kaiju/scripts/waxrun/config/cmeps_pin.patch

    # Fetch the CAM components.
    ./bin/git-fleximod update
    ./bin/git-fleximod status

Setting environment variables
************************************************

Set the ``WAXHOME`` environment variable to the path to your CAM clone. For
example:

.. code-block:: bash

    export WAXHOME=/path/to/your/cam

.. note::

    The ``WAXHOME`` environment variable is used by ``engage.py`` as the
    default CAM source tree for GRX runs. You can also type the path at the
    ``CAM source tree`` prompt.

    There is no separate build step for WACCM-X. ``engage.py`` creates and
    builds the WACCM-X cases from your CAM clone when it generates a GRX run.

.. note:: Documentation on the analysis tools is found
    :doc:`here </tools/index>`.
