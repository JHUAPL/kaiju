Building the CDF library on Ubuntu 20.04
========================================

Introduction
------------

This file describes how to build and install the CDF library on Ubuntu 20.04 for use with the ``kaiju`` software.

**NOTE**: A Ubuntu 20.04-specific subdirectory of the user home directory is used for this build.

As of this date (2023-07-11), the latest version of the CDF library is:

``cdf39_0-dist-all.tar.gz``

Building and installing the CDF library
---------------------------------------

.. code-block:: shell

   # Specify the name for this system.
   export HOST_SYSTEM=ubuntu-20.04

   # Specify and create the root of the build tree.
   export BUILD_ROOT=$HOME/$HOST_SYSTEM/local/cdf/3.9.0
   mkdir -p $BUILD_ROOT/src
   cd $BUILD_ROOT/src

   # Download the source tarball.
   wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_0/linux/cdf39_0-dist-all.tar.gz

   # Unpack the source code.
   tar xzvf cdf39_0-dist-all.tar.gz

   # Move into the code directory.
   cd cdf39_0-dist

   # Build the library using the default system GNU compiler.
   date; time make OS=linux ENV=gnu CURSES=no all >& make.out; date
   # Took 0m30.425s

   # Test the library.
   date; time make test >& make_test.out; date
   # Took 0m0.153s

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$BUILD_ROOT install >& make_install.out; date
   # Took 0m0.859s

   # Clean the build tree.
   date; time make clean >& make_clean.out; date
   # Took 0m0.142s


Using the CDF library
---------------------

To use this software, you must run the setup script:

.. code-block:: shell

   source $BUILD_ROOT/bin/definitions.B
