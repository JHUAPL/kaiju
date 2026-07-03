Building and installing the NASA CDF library on MacOS
=====================================================

Introduction
------------

This page describes how to build and install CDF 3.9.0 on an Intel-based Mac running MacOS Ventura for use with the kaiju software. Version 3.9.0 is the most recent version of the `NASA CDF library <https://cdf.gsfc.nasa.gov/>`_ as of 2023-07-05. The build is performed using the default ``gcc`` compiler provided by Apple.

Building and installing the CDF library
---------------------------------------

.. code-block:: shell

   # Specify the name for this machine.
   export HOST_SYSTEM=ventura

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

   # Patch the Makefile (required for a MacOS build).
   cp Makefile Makefile.orig
   cat Makefile.orig | sed 's/LDLIB=-L\$(XCODEDir)\/Platforms\/MacOSX.platform\/Developer\/SDKs\/MacOSX.sdk\/usr\/lib/LDLIB=-L\$(XCODEDir)\/SDKs\/MacOSX.sdk\/usr\/lib/' > Makefile

   # Build the library using the Apple compiler for Intel hardware.
   date; time make OS=macosx ENV=x86_64 all >& make.out; date
   # Took 0m33.062s

   # Test the library.
   date; time make test >& make_test.out; date
   # Took 0m2.396s

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$BUILD_ROOT install >& make_install.out; date
   # Took 0m5.681s

   # Clean the build tree.
   date; time make clean >& make_clean.out; date
   # Took 0m0.837s


Using the CDF library
---------------------

To use this library, you must run the setup script:

.. code-block:: shell

   source $HOME/ventura/local/cdf/3.9.0/bin/definitions.B
