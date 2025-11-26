Building the NASA CDF (Common Data Format) library
==================================================

The CDF library is required by SpacePy for full functionality. SpacePy is
installed as part of the process of setting up your Python environment for
``kaipy``. If your system already provides the CDF library, feel free to use
that version instead of building the code from source.


HPC systems (``derecho`` and ``aitken``)
------------------------------------------

This procedure builds the code under ``$HOME/local/src``, and installs the
compiled files under ``$HOME/local/cdf/3.8.1``, creating a user-only
installation. Substitute your desired build and installation locations for
these paths as appropriate.

**NOTE**: These instructions should work under any command shell.

.. code-block:: bash

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

In order to use the CDF library, you must set up your environment by
*sourcing* (not *running*) the setup script:

.. code-block:: bash

   # For sh/bash/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.B

   # For csh/tcsh/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.C

   # For ksh/compatible shells:
   source $HOME/local/cdf/3.8.1/bin/definitions.K
