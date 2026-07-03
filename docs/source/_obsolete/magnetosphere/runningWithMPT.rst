Running with MPT
================

Running MPI jobs with MPT requires some modifications to the PBS submission
scripts on any cluster (Cheyenne, Pleiades, likely others). Those
modifications are outlined below.

Note: when compiling with MKL in addition to MPT and MPI, MKL should be
enabled in the initial cmake command when initializing the build directory to
allow the code to compile

.. code-block:: shell

   cmake -DENABLE_MPI=ON -DENABLE_MKL=ON ..

Mpirun Command
--------------

Inside a typical PBS submission script using Intel MPI, when the mpi
application is executed the command would like something like this:

.. code-block:: shell

   mpirun ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out

mpirun is the mpi launcher used with Intel MPI (and several other MPI
libraries), and it runs the executable defined by the ${EXE} environment
variable. However, MPT has three differences. First, it uses the launcher
mpiexec_mpt instead of mpirun. Second, it requires explicit information about
how to pin the applications threads to processor cores (something that Intel
MPI can do automatically). And finally, MPT requires an additional environment
variable allowing us to create complex custom data types.

The same line above must be replaced by these lines in a submission script
using MPT:

.. code-block:: shell

   export MPI_TYPE_DEPTH=32
   export MPI_IB_CONGESTED=0
   export OMP_NUM_THREADS=36
   export NODEFILE=$TMPDIR/nodefile.$PBS_JOBID
   cp $PBS_NODEFILE $NODEFILE
   mpiexec_mpt ${KAIJUHOME}/scripts/preproc/correctOMPenvironment.sh ${NODEFILE} omplace ${EXE} ${RUNID}.xml ${JNUM} >> ${RUNID}.${JNUM}.out

Note that this does require the user to have set the TMPDIR environment
variable at the top of the script, which is recommended on Cheyenne. If the
script does not use TMPDIR, then NODEFILE should be set to simply
"nodefile.$PBS_JOBID".

This also requires the user to have set the KAIJUHOME environment variable to
the root directory of their kaiju repository (as mentioned in several other
wiki pages). This is so that the PBS script can find the
"correctOMPenvironment.sh" script. This is a custom replacement for the
standard "omplace" script. Because our submission scripts typically have
different numbers of MPI processes on different ranks (2 for gamera processes
and 1 for voltron processes), the standard omplace script can't properly
identify this and pin the processes threads.
