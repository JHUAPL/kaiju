Additional Automated Tests
==========================

As new code is written and added to Kaiju, as many automated tests as possible
should be written for it. There are not yet any strict rules about what parts
of the code require tests, but if there are any parts of the code which are
known to be error-prone, those are good targets for tests. Or at an absolute
bare minimum create a few tests that run the code end-to-end, verifying the
final output (a type of
`Smoke testing <http://softwaretestingfundamentals.com/smoke-testing/>`_)

Writing tests for pFUnit
------------------------

The `pFUnit website <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>`_
has documentation about how to write tests for pFUnit. In particular, I find
the
`examples <https://github.com/Goddard-Fortran-Ecosystem/pFUnit/tree/master/Examples/MPI_Halo/tests>`_
very helpful.

Adding to Existing Test Executables
-----------------------------------

The easiest way to add new tests is to add more tests to existing test
executables. Cmake looks in the kaiju/tests/ folders for file with the
extension .pf, and all such files are built into test executables. For
example, all .pf files in the kaiju/tests/gamera/ folder get compiled into the
"gamTests" binary. A new .pf file could be added to the gamera folder with new
tests for Gamera, and they would automatically be included and run.

You can use these templates for serial and MPI .pf files as a starting point:

:doc:`Serial Test Template <serialTestTemplate>`

:doc:`MPI Test Template <mpiTestTemplate>`

Creating a New Test Executable
------------------------------

If it isn't appropriate to add new tests to an existing test executable, then
you will need to create another set of tests. This should be done by creating
a new subfolder in the kaiju/tests/ folder, which will contain all of the
pFUnit .pf test files, as well as any required supporting files.

The cmake file located at kaiju/tests/CMakeLists.txt will need to be updated
with information about this new test executable. Each executable needs three
lines of information in the cmake file. The first line contains the name of
the new subfolder, and the search string to find the new .pf files. The second
line lists the dependencies (libraries) that this test executable will need in
order to compile and run. The third line provided the desired name of the
output binary, combines that with the source files and dependencies provided
on the previous two lines, and tells pFUnit what the maximum number of MPI
ranks required is to run these tests. For non-mpi tests this number should
simply be 1.

By way of example, this is the entry for the "gamMpiTests" executable,
containing tests for MPI enabled gamera.

.. code-block:: bash

   file(GLOB gamMpiTestFiles "${CMAKE_CURRENT_SOURCE_DIR}/gamera_mpi/*.pf")
   set(gamMpiTestLibs baselib gamlib basempilib gammpilib tgiclib)
   add_pfunit_ctest (gamMpiTests TEST_SOURCES ${gamMpiTestFiles} LINK_LIBRARIES ${gamMpiTestLibs} MAX_PES 64)
