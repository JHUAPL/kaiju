Development Roadmap
===================

Major Efforts
-------------

* Improve Multifluid support
* Improve timing\&profiling (timing info for all ranks, perhaps use
  Allinea)
* Incorporate Doxygen and link to wiki
* CHIMP unit tests - Particle pushing in analytic/numerical dipole,
  invariant conservation, etc
* Reorg/rewrite kaipy (better organization, and better speed of
  pipeline)
* RCM re-write to be object-oriented (Model/Grid/State), add unit
  tests, and support thread and MPI parallelism
* Grid optimization (use grid generation tools to improve default LFM
  grid?)
* GPU Implementation
* K-Decomposition of MHD solver
* Handle rotation axis by adding corotation potential in remix, and
  letting it tell gamera/rcm

Minor Efforts
-------------

* Hall MHD (Pull code from old omega repo)
* Create "default" configs to be used with XMLGenerator to create XMLs
  just-in-time
* Add simple reproducibility test/100 yard dash test
* Add auto-generation of allinea performance report data/auto-link to
  wiki page
* Use Allinea tools to do a extensive study of code
  performance/bottlenecks, investigate limitations on vectorization
  and OMP/MPI-scaling.
* Reorg wiki for application instead of serial/MPI, and assign
  responsibilities
* Auto-generation of spacecraft trajectory/MHD comparisons, pyspedas +
  sctrack.x
* Auto-generation of spacecraft trajectory/RCM comparisons
* Auto-generation of ampere/remix comparisons
* Auto-generate comparisons of sim pressure vs. empirical (TS/SST)
  pressure
* Handle remix/conductance mirror ratio based on Rin

Completed Efforts
-----------------

* MPI Implementation for "vanilla" MHD
* MPI Implementation for coupled Gamera/ReMIX (Voltron)
* Coupling implementation for Voltron "deep coupling", field-line
  tracing and heating for empirical/RCM pressure ingestion
* Make decomposed(MPI) face-centered data have a single master value
  at boundaries
* Asymmetric MPI transfers for Voltron shallow and deep coupling
* Reduce MPI transfers of unused regions of necessary variables
* Rename types module to gamtypes for consistency
* Handle units/msphutils globals better
* Added support for Intel Fortran 19&21
* Upgraded to Fortran 2008 MPI Interface
* CMake support added for multiple compilers and MPI libraries
* Git-LFS to store binary data in repo
* Add unit testing framework for python (similar to pFUnit for
  Fortran)
* Add intel mem/thread-checker autotest
