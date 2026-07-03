Building the CDF library on CentOS-Stream 9
===========================================

Introduction
------------

This file describes how to build and install the CDF library on CentOS-Stream 9 for use with the ``kaiju`` software.

**NOTE**: A CentOS-Stream 9-specific subdirectory of the user home directory is used for this build.

Building and installing the CDF library
---------------------------------------

.. code-block:: shell

   # Specify the name for this system.
   host_system=centos-stream-9

   # Specify and create the root of the build tree.
   build_root=$HOME/$host_system/local/cdf/3.9.0
   mkdir -p $build_root/src
   cd $build_root/src

   # Download the source tarball.
   wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_0/linux/cdf39_0-dist-all.tar.gz

   # Unpack the source code.
   tar xzvf cdf39_0-dist-all.tar.gz

   # Move into the code directory.
   cd cdf39_0-dist

   # Build the library using the default system GNU compiler.
   date; time make OS=linux ENV=gnu CURSES=no all >& make.out; date
   # Took 0m32.670s

   # Test the library.
   date; time make test >& make_test.out; date
   # Took 0m0.392s

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$build_root install >& make_install.out; date
   # Took 0m1.691s

   # Clean the build tree.
   date; time make clean >& make_clean.out; date
   # Took 0m0.378s


Using the CDF library
---------------------

To use this software, you must run the setup script:

.. code-block:: shell

   source $build_root/bin/definitions.B
