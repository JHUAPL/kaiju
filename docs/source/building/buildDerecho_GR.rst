Building the ``kaiju`` software on ``derecho`` for MAGE - Without TIEGCM (GR) 
=============================================================================================


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
module set:

.. code-block:: bash

    module --force purge
    module load ncarenv/23.06
    module load cmake/3.26.3
    module load craype/2.7.20
    module load intel/2023.0.0
    module load ncarcompilers/1.0.0
    module load cray-mpich/8.1.25
    module load hdf5-mpi/1.12.2
    module load conda/latest

.. important::

    You must use these exact versions of the modules to ensure the software
    compiles properly. If you use different versions of any of these modules,
    a successful build cannot be guaranteed. This module list is current as of
    **11 April 2025**, and is subject to change as the compute environment
    changes.

Build the ``kaiju`` software
----------------------------

These instructions show how to build the MPI version of the ``kaiju``
software. The MPI version is built in the subdirectory ``build_mpi``
under the ``kaiju`` source code directory. In practice, you can place the
build directory in any convenient location.

.. code-block:: bash

    # Move to your kaiju clone.
    cd /path/to/kaiju

    # Create the build directory and enter it.
    mkdir build_mpi
    cd build_mpi

    # Run cmake to create the Makefile, saving output.
    # NOTE: The FC definition is *required* for proper cmake operation.
    FC=`which ifort` cmake -DENABLE_MPI=ON .. >& cmake.out
    # If you are building the GTR version of the code, use the following cmake command instead:

    # You can pick one compile target below or compile all of them, if you'd like

    # Compile the MAGE model for geospace simulations
    make -j4 voltron_mpi.x >& make-voltron.out

    # Compile the GAMERA-helio model for inner heliosphere simulations
    make -j4 gamhelio_mpi.x >& make-gamhelio.out

    # Compile analysis tools
    make -j4 calcdb.x chop.x sctrack.x slice.x >& make-analysis.out

    
When finished, your build directory will contain a ``bin``
subdirectory which will contain the compiled ``kaiju`` executables.

.. note:: Documentation on the analysis tools is found
    :doc:`here </tools/index>`.

