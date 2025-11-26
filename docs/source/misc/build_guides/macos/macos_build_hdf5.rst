Building and installing the HDF5 library on MacOS
=================================================

Introduction
------------

This file describes how to build and install HDF5 on an Intel-based Mac running MacOS Ventura for use with the kaiju software.

These instructions assume that the Intel compiler suite has already been installed. If HDF5 is built with gcc, then it will not link with the kaiju code, which is built with the Intel compiler.

Building and installing the HDF5 library
----------------------------------------

.. code-block:: shell

   # Use the Intel compiler.
   export INTEL_HOME=/opt/intel
   export PATH=$INTEL_HOME/oneapi/compiler/latest/mac/bin/intel64:$PATH

   # Specify the name for this machine.
   export HOST_SYSTEM=ventura

   # Specify and create the root of the build tree.
   export BUILD_ROOT=$HOME/$HOST_SYSTEM/local/hdf5/1.14.1-2
   mkdir -p $BUILD_ROOT/src
   cd $BUILD_ROOT/src

   # Download the source tarball.
   wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.1/src/hdf5-1.14.1-2.tar.gz

   # Unpack the source code.
   tar xzvf hdf5-1.14.1-2.tar.gz

   # Move into the code directory.
   cd hdf5-1.14.1-2

   # Configure the library with Fortran support.
   date; time \
     ./configure \
     FC=ifort \
     --prefix=$BUILD_ROOT \
     --enable-fortran \
     >& configure.out; date
   # Took 2m9.603s

   # Compile the library.
   date; time make >& make.out; date
   # Took 10m33.074s

   # Test the library.
   date; time make check >& make_check.out; date
   # Took 15m17.099s

   # Install the library.
   date; time make install >& make_install.out
   # Took 0m8.438s

   # Clean the build tree.
   date; time make clean >& make_clean.out
   # Took 0m7.811s


Using the HDF5 library
----------------------

To use the HDF library, you must set environment variables:

.. code-block:: shell

   export HDF5_DIR=$HOME/ventura/local/hdf5/1.12.1
   export PATH=$HDF5_DIR/bin:$PATH
   export HDF5_INCLUDE_DIRS=$HDF5_DIR/include
   export HDF5_LIBRARIES=$HDF5_DIR/lib
   export CPATH=$HDF5_INCLUDE_DIRS
