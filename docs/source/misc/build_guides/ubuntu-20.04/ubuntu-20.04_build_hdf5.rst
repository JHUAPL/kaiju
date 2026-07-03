Building and installing the HDF5 library on Ubuntu 20.04
========================================================

Introduction
------------

This file describes how to build and install HDF5 on an Intel-based system running Ubuntu 20.04 for use with the kaiju software.

These instructions assume that the Intel compiler suite has already been installed. If HDF5 is built with gcc, then it will not link with the kaiju code, which is built with the Intel compiler.

Building and installing the HDF5 library
----------------------------------------

.. code-block:: shell

   # Specify the name for this machine.
   export HOST_SYSTEM=ubuntu-20.04

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

   # Configure the Intel compilers.
   source /opt/intel/oneapi/setvars.sh

   # Configure the library with Fortran support.
   date; time \
     ./configure \
     FC=ifort \
     --prefix=$BUILD_ROOT \
     --enable-fortran \
     >& configure.out; date
   # Took 0m26.314s

   # Compile the library.
   date; time make >& make.out; date
   # Took 8m43.088s

   # Test the library.
   date; time make check >& make_check.out; date
   # Took 9m9.848s

   # Install the library.
   date; time make install >& make_install.out
   # Took 0m4.232s

   # Clean the build tree.
   date; time make clean >& make_clean.out
   # Took 0m1.874s


Using the HDF5 library
----------------------

To use the HDF library, you must set environment variables:

.. code-block:: shell

   export HDF5_DIR=$HOME/ubuntu-20.04/local/hdf5/1.14.1-2
   export PATH=$HDF5_DIR/bin:$PATH
   export HDF5_INCLUDE_DIRS=$HDF5_DIR/include
   export HDF5_LIBRARIES=$HDF5_DIR/lib
   export CPATH=$HDF5_INCLUDE_DIRS
   export INCLUDE="-I$HDF5_INCLUDE_DIRS"
