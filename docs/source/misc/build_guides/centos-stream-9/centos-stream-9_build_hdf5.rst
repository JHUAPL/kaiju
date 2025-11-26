Building and installing the HDF5 library on CentOS-Stream 9
===========================================================

Introduction
------------

This file describes how to build and install HDF5 on an Intel-based system running CentOS-Stream 9 for use with the kaiju software.

These instructions assume that the Intel compiler suite has already been installed. If HDF5 is built with gcc, then it will not link with the kaiju code, which is built with the Intel compiler.

Building and installing the HDF5 library
----------------------------------------

.. code-block:: shell

   # Specify the name for this machine.
   host_system=centos-stream-9

   # Specify and create the root of the build tree.
   build_root=$HOME/$host_system/local/hdf5/1.14.1-2
   mkdir -p $build_root/src
   cd $build_root/src

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
     --prefix=$build_root \
     --enable-fortran \
     >& configure.out; date
   # Took 0m39.463s

   # Compile the library.
   date; time make >& make.out; date
   # Took 9m29.899s

   # Test the library.
   date; time make check >& make_check.out; date
   # Took 10m8.337s

   # Install the library.
   date; time make install >& make_install.out
   # Took 0m5.598s

   # Clean the build tree.
   date; time make clean >& make_clean.out
   # Took 0m3.342s


Using the HDF5 library
----------------------

To use the HDF library, you must set environment variables:

.. code-block:: shell

   export HDF5_DIR=$HOME/centos-stream-9/local/hdf5/1.14.1-2
   export PATH=$HDF5_DIR/bin:$PATH
   export HDF5_INCLUDE_DIRS=$HDF5_DIR/include
   export HDF5_LIBRARIES=$HDF5_DIR/lib
   export CPATH=$HDF5_INCLUDE_DIRS
   export INCLUDE="-I$HDF5_INCLUDE_DIRS"
