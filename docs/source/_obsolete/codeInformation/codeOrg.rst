Code Organization
=================

The Kaiju repository holds source code, pre- and post-processing
scripts, and run configurations for various the Kaiju elements (Gamera
MHD, CHIMP test-particle and interpolation, ReMIX ionospheric solver,
etc).

Directories from $KAIJUHOME


* src: Main Fortran source files

  * base: Common support libraries (IO, math, string parsing,
    timing, etc)

    * base/kdefs.F90: Highest-level parameter definitions
      (precision, constants, etc)
    * defs: Parameter definitions for various Kaiju elements
      (Gamera, CHIMP, etc)
    * types: Type definitions for main data structures for Kaiju
      elements

* gamera/chimp/remix: Source files for various Kaiju elements

  * ICs/: Problem files for specifying initial conditions, can be
    changed using CMake variables

* voltron: Files for running coupled magnetosphere-ionosphere runs
* drivers: Fortran program files, contain main heartbeat loops for
  executables
* examples: Run parameters for various simple example runs
* places: Run scripts and configuration informatin for running on
  specific environments (e.g. Cheyenne)
* kaipy: Python modules to support pre- and post-processing scripts
* scripts: Command line python scripts for pre- and post-processing
* cmake: CMake defaults

  * user.cmake: Optional file (not included by default) to specify
    user over-rides to CMake variables

* tests: pFUnit scripts and cases for unit testing
