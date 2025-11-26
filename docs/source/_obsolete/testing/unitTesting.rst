Unit Testing
============

Overview
--------

Kaiju supports unit testing through the pFUnit library. In order to compile
and run the built in tests, the user must add pFUnit to their kaiju
repository.

Getting pFUnit
--------------

On Cheyenne, there are pre-built binaries for pFUnit available at
``/glade/p/hao/msphere/gamshare/``. There are a few different versions of
pFUnit built with different versions of different compilers. For example, you
can copy the build of pFUnit-4.1.5 made with Intel Fortan 18 from the folder
``/glade/p/hao/msphere/gamshare/pFUnit-4.1.5/ifort-18/``.

Within a specific build folder you will find four subfolders named
``FARGPARSE-X.X``, ``GFTL-X.X``, ``GFTL_SHARED-X.X``, and ``PFUNIT-X.X``.
These four folders should either be copied or linked into your
/kaiju/external/ folder. No additional setup is required, cmake will
automatically find pFUnit once these four folders are properly located in
your external folder (note that these four folders cannot be in a subfolder
beneath the external folder, they must be directly within the external folder
itself).

If you are not on Cheyenne, or cannot use these binares for whatever reason,
there are instructions for how to :doc:`Build pFUnit <buildpFUnit>`.

Building the tests
------------------

Once pFUnit has been installed into the external folder of the repository, the
testing executables can be built with the cmake command "make allTests".

** Under Intel Fortran 21, the FC environment variable had to be manually set
to mpiifort in order for compilation to succeed.

Running the tests
-----------------

Once the test executables have been built, there are two submission scripts in
the kaiju/tests/ folder that can be used to run them on Cheyenne. The script
``runNonCaseTests.pbs`` runs a faster subset of tests which validate specific
pieces of the code. The script ``runCaseTests.pbs`` runs a slower set of tests
which run full simulations and then verify the final states. These submission
scripts should be submitted from the folder containing the test executables,
usually the build/bin folder.

Verifying test results
----------------------

When the two jobs mentioned above are run, they create output files labeled
with the job ID. They will look something like 1nonCaseTests.o12345671 and
``caseTests.o1234567`` respectively, but with different numbers. These files
contain summaries of all of the test executables that ran, and should contain
a series of lines that look like this:

.. code-block:: bash

   Gamera Tests Complete
    OK
    (10 tests)

If any of the tests don't look like this, it's an indication that some of the
tests have failed, and the more detailed output file for that executable
should be examined.
