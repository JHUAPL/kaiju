Building the HDF5 library
=========================

Introduction
------------

The HDF5 library is required by the ``kaiju`` for data storage. This page describes how to build and install HDF5 for use with the kaiju software. If your system already provides the HDF5 library, feel free to use that version instead of building the code from source.

MacOS
-----

These instructions describe how to build and install the HDF5 library (version 1.12.1) for use with the ``kaiju`` software on MacOS systems. These instructions were developed on a 2019 MacBook Pro with an i9 processor and 64 GB of RAM. The machine was running MacOS 11.6.5 (Big Sur), but these instructions should work with minimal changes on other versions of MacOS, or using other versions of the HDF5 source code. Several of the commands below display execution times - the times you observe may be more or less than the times shown below, depending on your hardware and software versions.

**Note**: The HDF5 software may be available through one of several MacOS package managers (HomeBrew, MacPorts, Fink, ...). To save time, you should try to install the HDF5 software using a package manager before trying to build the HDF5 software from source code. However, HDF5 built with gfortran (i.e. from HomeBrew) will not be compatible if you want to build kaiju with Intel Fortran. You must build HDF5 from source with ifort in order to build kaiju with ifort. 

**Note**: This procedure uses the Intel C and Fotran compilers. The code is built under ``$HOME/local/src``, and is installed under ``$HOME/local/hdf/1.12.1``, creating a user-only installation. Substitute your desired build and installation locations for these paths as appropriate.

**NOTE**: These instructions should work under any command shell.

.. code-block:: shell

   # Make a build tree.
   cd $HOME
   mkdir -p local/src
   cd local/src

   # Download the source tarball from (requires free account login):
   # https://www.hdfgroup.org/downloads/hdf5/

   # Unpack the source code (your version number may vary).
   tar xzvf hdf5-1.12.1.tar.gz

   # Use the Intel compilers.
   # For sh/bash/zsh-compatible shells:
   export PATH=/opt/intel/oneapi/compiler/latest/mac/bin/intel64:$PATH
   # For csh/tcsh-compatible shells:
   # setenv PATH /opt/intel/oneapi/compiler/latest/mac/bin/intel64:$PATH

   # Configure the HDF5 software with Fortran support.
   cd hdf5-1.12.1
   date; time \
     ./configure \
     FC=ifort \
     --prefix=$HOME/local/hdf5/1.12.1 \
     --enable-fortran \
     >& configure.out
   # Took 2m13.280s.
   # Examine configure.out for errors.

   # Build the HDF5 software.
   date; time make >& make.out
   # Took 13m22.799s.
   # Examine make.out for errors.

   # Test the HDF5 software.
   date; time make check >& make_check.out
   # Took 19m9.099s.
   # Check make_check.out for errors.

   # Install the HDF5 software.
   date; time make install >& make_install.out
   # Took 0m8.130s.
   # Check make_install.out for errors.


Linux
-----

These instructions describe how to build and install the HDF5 library (version 1.12.1) for use with the ``kaiju`` software on Linux systems. These instructions were developed on a Ubuntu virtual machine, but these instructions should work with minimal changes on other versions of Linux, or using other versions of the HDF5 source code. Several of the commands below display execution times - the times you observe may be more or less than the times shown below, depending on your hardware and software versions.

**Note**: The HDF5 software may be available through one of several Linux package managers (``yum``, ``dnf``, ``rpm``, ...). To save time, you should try to install the HDF5 software using a package manager before trying to build the HDF5 software from source code. You may require assistance from a system administrator to install the software using a package manager.

**Note**: This procedure uses the Intel C and Fotran compilers. The code is built under ``$HOME/local/src``, and is installed under ``$HOME/local/hdf/1.12.1``, creating a user-only installation. Substitute your desired build and installation locations for these paths as appropriate.

**NOTE**: These instructions should work under any command shell.

.. code-block:: shell

   # Make a build tree.
   cd $HOME
   mkdir -p local/src
   cd local/src

   # Download the source tarball from (requires free account login):
   # https://www.hdfgroup.org/downloads/hdf5/

   # Unpack the source code (your version number may vary).
   tar xzvf hdf5-1.12.1.tar.gz

   # Use the Intel compilers.
   # For sh/bash-compatible shells:
   export PATH=/opt/intel/oneapi/compiler/latest/mac/bin/intel64:$PATH
   # For csh/tcsh-compatible shells:
   setenv PATH /opt/intel/oneapi/compiler/latest/mac/bin/intel64:$PATH

   # Configure the HDF5 software with Fortran support.
   cd hdf5-1.12.1
   date; time \
     ./configure \
     FC=ifort \
     --prefix=$HOME/local/hdf5/1.12.1 \
     --enable-fortran \
     >& configure.out
   # Took XXX.
   # Examine configure.out for errors.

   # Build the HDF5 software.
   date; time make >& make.out
   # Took XXX.
   # Examine make.out for errors.

   # Test the HDF5 software.
   date; time make check >& make_check.out
   # Took XXX.
   # Check make_check.out for errors.

   # Install the HDF5 software.
   date; time make install >& make_install.out
   # Took XXX.
   # Check make_install.out for errors.


HPC systems
-----------

HPC systems typically provide software via a user-selectable system, such as the ``module`` command. The instructions below were developed on ``pleiades``, and thus presume use of the ``module`` command. The ``module`` commands for ``cheyenne`` and related systems are slightly different, and will be provided here as they are documented.

``pleiades``
^^^^^^^^^^^^^^^^

The ``pleiades`` system provides multiple versions of the HDF5 software. The desired version can be loaded with a ``module`` command. The primary concern is loading the correct version - the HDF5 library you use must be built with the Intel compilers, and be built for the proper version of the ``kaiju`` software (serial or MPI - more information on the differences if provided on the `build <build>`_ page). If desired, you may also build the HDF5 software from source code as described in the `Linux <LINK>`_ section.

To load the Intel-compiled HDF5 library and its prerequisites for use with the serial ``kaiju`` software on ``pleiades`` (this works under any shell):

.. code-block:: shell

   module load comp-intel/2020.4.304
   module load hdf5/1.8.18_serial


To load the Intel-compiled HDF5 library and its prerequisites for use with the MPI ``kaiju`` software on ``pleiades`` (this works under any shell):

.. code-block:: shell

   module load comp-intel/2020.4.304
   module load mpi-hpe/mpt.2.23
   module load hdf5/1.8.18_mpt
