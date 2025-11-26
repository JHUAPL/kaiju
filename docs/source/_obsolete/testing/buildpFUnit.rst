Build pFUnit
============

pFUnit can be downloaded from
`https://github.com/Goddard-Fortran-Ecosystem/pFUnit <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>`_.
It is recommended that the user get the latest version of the master branch.

The pFUnit repository contains instructions for building pFUnit. pFUnit should
be built with cmake, and with support for both MPI and OPENMP enabled.

Once pFUnit is built, the compiled files (along with any external dependencies
pFUnit needs) must be placed into the "external" directory in the main kaiju
repository. As of pFUnit version 4.1, this meant that 4 folders were copied
from the pFUnit "installed" folder: "PFUNIT-4.1", "GFTL-1.1",
"GFTL_SHARED-1.0", and "FARGPARSE-0.9". These 4 folders are all placed in the
"external" folder.

NOTE
====

Due to an incompatibility, the file ``add_pfunit_ctest.cmake`` in
PFUNIT-4.1/include/ must be modified:

Line 50 must be changed from ``set (test_sources_f90)`` to
``set (test_sources_f90 "")``

Line 56 containing ``list (APPEND test_sources_f90 ${f90_file})`` must be
commented out or deleted

These line numbers may change in other versions of pFUnit.
