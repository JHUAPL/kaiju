Building the NASA CDF (Common Data Format) library
==================================================

The CDF library is required by SpacePy for full functionality. SpacePy is installed as part of the process of setting up your Python environment, which is described `here <./install_python.md>`_. If your system already provides the CDF library, feel free to use that version instead of building the code from source.

MacOS
-----

These instructions presume the use of MacOS Ventura. The instructions should also work for other versions of MacOS. These instructions illustrate steps required to build and insta;; version 3.9.0 of the CDF library. Modify these instructions as appropriate if you use a different version.

This procedure will use the default system-provided C compiler (``/usr/bin/gcc``).

This procedure builds and installs the code under ``$HOME/local``, creating a user-only installation. Substitute your desired build and installation locations for these paths as appropriate.

.. code-block:: shell

   # Make a build tree.
   cd $HOME
   mkdir -p local/cdf/3.9.0/src
   cd local/cdf/3.9.0/src

   # Download the source tarball.
   # NOTE: Source tarball is in a linux directory, not Mac.
   wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest/linux/cdf39_0-dist-all.tar.gz

   # Unpack the source code.
   tar xzvf cdf39_0-dist-all.tar.gz

   # Build the library using the Apple compiler for Intel hardware.
   cd cdf39_0-dist
   date; time make OS=macosx ENV=x86_64 all >& make.out
   # Took about 22 s on a MacBook Pro i9/2019.
   # Examine make.out for errors.

   # Test the library.
   date; time make test >& make_test.out
   # Took about 5 s on a MacBook Pro i9/2019.
   # Examine make_test.out for errors.

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$HOME/local/cdf/3.8.1 install >&
     make_install.out
   # Took about 8 s on a MacBook Pro i9/2019.
   # Examine make_install.out for errors.


In order to use the newly-compiled library, you must set up your environment (path variables and aliases) by "sourcing" the setup script:

.. code-block:: shell

   # For sh/bash/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.B

   # For csh/tcsh/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.C

   # For ksh/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.K


Linux and HPC systems
---------------------

These instructions should work for any Linux distribution. These instructions were tested on ``pleiades``, as well as an Ubuntu-based virtual machine.

This procedure will use the default system-provided C compiler (``/usr/bin/gcc``).

This procedure builds the code under ``$HOME/local/src``, and installs the compiled files under ``$HOME/local/cdf/3.8.1``, creating a user-only installation. Substitute your desired build and installation locations for these paths as appropriate.

**NOTE**: These instructions should work under any command shell.

.. code-block:: shell

   # Make a build tree.
   cd $HOME
   mkdir -p local/src
   cd local/src

   # Download the source tarball.
   wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/linux/cdf38_1-dist-all.tar.gz
   # OR:
   # curl https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/linux/cdf38_1-dist-all.tar.gz -o cdf38_1-dist-all.tar.gz

   # Unpack the source code.
   tar xzvf cdf38_1-dist-all.tar.gz

   # Build the library using the default system compiler (presumed to be gcc).
   cd cdf38_1-dist
   date; time make OS=linux ENV=gnu all >& make.out
   # Took about 33 s on pleiades.
   # Examine make.out for errors.

   # Test the library.
   date; time make test >& make_test.out
   # Took about < 1 s on pleiades.
   # Examine make_test.out for errors.

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$HOME/local/cdf/3.8.1 install >& make_install.out
   # Took about < 1 s on pleiades.
   # Examine make_install.out for errors.


In order to use the newly-compiled library, you must set up your environment (path variables and aliases) by "sourcing" the setup script:

.. code-block:: shell

   # For sh/bash/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.B

   # For csh/tcsh/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.C

   # For ksh/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.K
