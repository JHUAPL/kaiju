Building the ``kaiju`` software on ``pleiades`` for MAGE - With TIEGCM (GTR)
==============================================================================================

Introduction
------------

This page provides instructions for building the ``kaiju`` software on the
``pleiades`` supercomputer. These instructions assume that you have cloned the
``kaiju`` repository.

Prepare your software environment
---------------------------------

Like most HPC systems, ``pleiades`` uses the ``module`` system to manage the
versions of software packages available to the user. When you log in to
``pleiades``, no modules are loaded by default:

.. code-block:: bash

   module list
   No Modulefiles Currently Loaded.

Start by purging any currently-loaded modules, then loading the following
module set for MAGE runs coupled with TIEGCM (known as "GTR"):

.. warning::

    The GTR currently required custom built NetCDF and ESMF modules on ``pleiades``. If you need to
    run GTR, you will need access to ``/home7/nrao3/local3`` and ``/nobackup/nrao3/tiegcm/tiegcm3.0/data``,
    please reach out to ``nikhilr@ucar.edu`` with the following:
    
        - Your Pleiades username
        - Your Name
        - Your Institution

.. code-block:: bash

    module --force purge
    
    module use -a /nasa/modulefiles/testing
    module use -a /swbuild/analytix/tools/modulefiles
    module load nas
    module load comp-intel/2020.4.304
    module load mpi-hpe/mpt.2.30
    module load szip/2.1.1
    module load hdf5/1.12.3_mpt
    module load miniconda3/v4

    export FC=mpif90
    export CC=mpicc
    export CXX=mpicxx

    export PREFIX=/home7/nrao3/local3
    export LIBRARY_PATH=${LIBRARY_PATH}:$PREFIX/lib
    export LD_LIBRARY_PATH=$LIBRARY_PATH
    export CPATH=$PREFIX/include
    export PATH=${PATH}:$PREFIX/bin

.. important::

    You must use these exact versions of the modules to ensure the software
    compiles properly. If you use different versions of any of these modules,
    a successful build cannot be guaranteed. This module list is current as of
    **11 April 2025**, and is subject to change as the compute environment
    changes.

Build the ``kaiju`` software
----------------------------

These instructions show how to build the MPI version of the ``kaiju``
software. The GTR version is built in the subdirectory ``build_gtr``
under the ``kaiju`` source code directory. In practice, you can place the
build directory in any convenient location.

.. code-block:: bash

    # Move to your kaiju clone.
    cd /path/to/kaiju

    # Create the GTR build directory and enter it.
    mkdir build_gtr
    cd build_gtr

    # Run cmake to create the Makefile, saving output.
    # NOTE: The FC definition is *required* for proper cmake operation.
    FC=`which ifort` cmake -DENABLE_MPI=ON .. >& cmake.out

    # You can pick one compile target below or compile all of them, if you'd like

    # Compile the MAGE model for geospace simulations
    make -j4 voltron_mpi.x >& make-voltron.out

    # Compile the GAMERA-helio model for inner heliosphere simulations
    make -j4 gamhelio_mpi.x >& make-gamhelio.out

    # Compile analysis tools
    make -j4 calcdb.x chop.x sctrack.x slice.x >& make-analysis.out

.. warning::

    Pleiades faces issues while loading the ``hdf5/1.12.3_mpt`` module at times.
    If you aren't able to build ``kaiju`` with the above module set, unload 
    the ``hdf5/1.12.3_mpt`` module and load it again:

    .. code-block:: bash

        module unload hdf5/1.12.3_mpt
        module load hdf5/1.12.3_mpt

    
When finished, your build directory will contain a ``bin``
subdirectory which will contain the compiled ``kaiju`` executables.


Build the ``tiegcm`` software
-----------------------------

`TIEGCM <https://tiegcm-docs.readthedocs.io/>`_ is a comprehensive, first-principles, three-dimensional, 
non-linear representation of the coupled thermosphere and ionosphere system that includes a self-consistent solution 
of the middle and low-latitude dynamo field. 

Getting the TIE-GCM source code
************************************************

The ``TIE-GCM`` source code can be obtained by cloning the ``TIE-GCM`` repository
on GitHub:

.. code-block:: bash

    git clone https://github.com/NCAR/tiegcm.git

Setting environment variables
************************************************

Add paths to ``TIEGCMHOME`` and ``TIEGCMDATA``.
For example:

.. code-block:: bash

    export TIEGCMHOME=/path/to/your/tiegcm
    export TIEGCMDATA=/path/to/your/tiegcm/data

.. note::

    The ``TIEGCMHOME`` and ``TIEGCMDATA`` environment variables are required
    for running the GTR model. They should point to the TIEGCM source code
    directory and the TIEGCM data directory, respectively.

    The TIEGCMDATA directory is located in the following locations:
        - On ``derecho``: ``/glade/campaign/hao/itmodel/tiegcm3.0/new_data``
        - On ``pleiades``: ``/nobackup/nrao3/tiegcm/tiegcm3.0/data``
        - The required data files can be downloaded from the NCAR Globus endpoint using the following link: `TIEGCM Data Files <https://app.globus.org/file-manager?origin_id=b2502c58-c3eb-470f-86d4-cbdcd0aeb6c8&origin_path=%2F>`_


Resolution guide for TIEGCM
************************************************
Two TIEGCM executables are required for running the GTR model:

    - TIEGCM Standalone
        This is the TIEGCM code that runs independently and is used for initialization of the model.
    - TIEGCM Coupled
        This is the TIEGCM code that runs in a coupled mode with the GR model, providing 
        real-time updates to the thermosphere and ionosphere conditions during the simulation.

Depending on the Gamera resolution you will need to compile different TIEGCM resolution executables:
    - For a ``D`` run
        - TIEGCM Standalone: horires = 2.5, vertres = 0.25(1/4), mres = 2
        - TIEGCM Coupled: horires = 2.5, vertres = 0.25(1/4), mres = 2
    - For a ``Q`` run
        - TIEGCM Standalone: horires = 2.5, vertres = 0.25(1/4), mres = 2
        - TIEGCM Coupled: horires = 1.25, vertres = 0.125(1/8), mres = 1
    - For a ``O`` run
        - TIEGCM Standalone: horires = 1.25, vertres = 0.125(1/8), mres = 1
        - TIEGCM Coupled: horires = 0.625, vertres = 0.0625(1/16), mres = 0.5


The TIEGCM code is built using the ``tiegcmrun`` script, which is provided in
the ``tiegcm`` code repository. The script is provided in the
``tiegcm/tiegcmrun`` directory. More information on ``tiegcmrun.py`` can be found
in the `TIEGCM Quick Start Guide <https://tiegcm-docs.readthedocs.io/en/latest/tiegcm/quickstart.html>`_.

.. important::
    Make sure to load the modules lised in the ``kaiju`` build instructions
    before running the ``tiegcmrun`` script.

Build guide for TIEGCM code for a ``Q`` run:
#########################################################################################

We will use ``tiegcmrun`` script to build the code which requires the minimum amount
of input from the user. At each prompt, you can either type in a value, or hit 
the :kbd:`Return` key to accept the default value (shown in square brackets at 
the end of the prompt).

1. First we will create a directory for the "Q" TIEGCM build in the TIEGCMHOME directory.

.. code-block:: bash

    cd $TIEGCMHOME
    mkdir tiegcm_build_Q
    cd tiegcm_build_Q

2. Next, we will build the standalone TIEGCM executable by running the ``tiegcmrun.py`` script with the ``-oc``.

.. note::

    The ``-oc`` option stands for "only compile", which means that the script will only compile the code and not run it.
    Since the Gamera resolution is ``Q``, we will set the following:

        - horizontal resolution for the standalone TIEGCM to 2.5 degrees
        - vertical resolution to 0.25(1/4) scale height
        - magnetic grid resolution to 2 degrees

.. code-block:: bash

    $TIEGCMHOME/tiegcmrun/tiegcmrun.py -oc
    Instructions:
    -> Default Selected input parameter is given in GREEN
    -> Warnings and Information are given in YELLOW
    -> Errors are given in RED
    -> Valid values (if any) are given in brackets eg. (value1 | value2 | value3) 
    -> Enter '?' for any input parameter to get a detailed description


    Run Options:
    User Mode = BASIC
    Compile = True
    Execute = False
    Coupling = False


    Name of HPC system (derecho|pleiades|linux) [pleiades]: 
    Standalone Executable [/glade/derecho/scratch/nikhilr/tiegcm_build_Q/exec/tiegcm.exe]: 
    Horizontal Resolution (Deg) (5.0|2.5|1.25|0.625) [2.5]: 
    Vertical Resolution (Scale Height) (1/2|1/4|1/8|1/16) [1/4]: 
    Magnetic grid resolution (Degree) (2|1|0.5) [2]: 

After these inputs, the script will compile the TIEGCM code and create the standalone executable and should output something like this:

.. code-block:: bash

    ..
    .. 
    gmake[1]: Leaving directory '/nobackup/nrao3/tiegcm/tiegcm_build_Q/exec'    
    Executable copied from /nobackup/nrao3/tiegcm/tiegcm_build_Q/exec/tiegcm.exe to /nobackup/nrao3/tiegcm/tiegcm_build_Q/stdout

3. Next, we will build the coupled TIEGCM executable by running the ``tiegcmrun.py`` script with the ``-oc`` and ``-co`` options.

.. note::

    The ``-co`` option stands for "coupled", which means that the script will compile the code for the coupled TIEGCM executable.
    Since the Gamera resolution is ``Q``, we will set the following:

        - horizontal resolution for the coupled TIEGCM to 1.25 degrees
        - vertical resolution to 0.125(1/8) scale height
        - magnetic grid resolution to 1 degrees

.. code-block:: bash

    $TIEGCMHOME/tiegcmrun/tiegcmrun.py -oc -co
    Instructions:
    -> Default Selected input parameter is given in GREEN
    -> Warnings and Information are given in YELLOW
    -> Errors are given in RED
    -> Valid values (if any) are given in brackets eg. (value1 | value2 | value3) 
    -> Enter '?' for any input parameter to get a detailed description

    Run Options:
    User Mode = BASIC
    Compile = True
    Execute = False
    Coupling = True

    Name of HPC system (derecho|pleiades|linux) [pleiades]: 
    Coupled Executable [/glade/derecho/scratch/nikhilr/tiegcm_build/exec/tiegcm.x]: 
    Horizontal Resolution (Deg) (5.0|2.5|1.25|0.625) [2.5]: 1.25
    Vertical Resolution (Scale Height) (1/2|1/4|1/8|1/16) [1/8]: 
    Magnetic grid resolution (Degree) (2|1|0.5) [1]: 

After these inputs, the script will compile the TIEGCM code and create the coupled executable and should output something like this:

.. code-block:: bash

    ..
    ..
    gmake[1]: Leaving directory '/nobackup/nrao3/tiegcm/tiegcm_build_Q/exec'    
    Executable copied from /nobackup/nrao3/tiegcm/tiegcm_build_Q/exec/tiegcm.x to /nobackup/nrao3/tiegcm/tiegcm_build_Q/stdout

4. You should now see the following files in your run directory:

.. code-block:: bash

    ls
    exec  hist  stdout

The executables are located in the ``stdout`` directory, and the stdout files

.. code-block:: bash

    ls stdout
    defs.h  tiegcm.exe  tiegcm.x

Build guide for TIEGCM code for a ``D`` run on ``derecho``:
#########################################################################################

1. First we will create a directory for the "D" TIEGCM build in the TIEGCMHOME directory.

.. code-block:: bash

    cd $TIEGCMHOME
    mkdir tiegcm_build_D
    cd tiegcm_build_D

2. Next, we will build the standalone TIEGCM executable by running the ``tiegcmrun.py`` script with the ``-oc``.

.. note::
    
    Since the Gamera resolution is ``D``, we will set the following:

        - horizontal resolution for the standalone TIEGCM to 2.5 degrees
        - vertical resolution to 0.25(1/4) scale height
        - magnetic grid resolution to 2 degrees

.. code-block:: bash

    $TIEGCMHOME/tiegcmrun/tiegcmrun.py -oc

3. Next, we will build the coupled TIEGCM executable by running the ``tiegcmrun.py`` script with the ``-oc`` and ``-co`` options.

.. note::

    Since the Gamera resolution is ``D``, we will set the following:

        - horizontal resolution for the coupled TIEGCM to 2.5 degrees
        - vertical resolution to 0.25(1/4) scale height
        - magnetic grid resolution to 2 degrees

.. code-block:: bash

    $TIEGCMHOME/tiegcmrun/tiegcmrun.py -oc -co

4. The executables are located in the ``stdout`` directory, and the stdout files

.. code-block:: bash

    ls stdout
    defs.h  tiegcm.exe  tiegcm.x


Build guide for TIEGCM code for a ``O`` run on ``derecho``:
#########################################################################################

1. First we will create a directory for the "D" TIEGCM build in the TIEGCMHOME directory.
.. code-block:: bash

    cd $TIEGCMHOME
    mkdir tiegcm_build_O
    cd tiegcm_build_O

2. Next, we will build the standalone TIEGCM executable by running the ``tiegcmrun.py`` script with the ``-oc``.

.. note::

    Since the Gamera resolution is ``O``, we will set the following:

        - horizontal resolution for the standalone TIEGCM to 1.25 degrees
        - vertical resolution to 0.125(1/8) scale height
        - magnetic grid resolution to 1 degrees

.. code-block:: bash

    $TIEGCMHOME/tiegcmrun/tiegcmrun.py -oc

3. Next, we will build the coupled TIEGCM executable by running the ``tiegcmrun.py`` script with the ``-oc`` and ``-co`` options.

.. note::
    Since the Gamera resolution is ``O``, we will set the following:

        - horizontal resolution for the coupled TIEGCM to 0.625 degrees
        - vertical resolution to  0.0625(1/16) scale height
        - magnetic grid resolution to 0.5 degrees

.. code-block:: bash

    $TIEGCMHOME/tiegcmrun/tiegcmrun.py -oc -co

4. The executables are located in the ``stdout`` directory, and the stdout files

.. code-block:: bash

    ls stdout
    defs.h  tiegcm.exe  tiegcm.x

.. note:: Documentation on the analysis tools is found
    :doc:`here </tools/index>`.


