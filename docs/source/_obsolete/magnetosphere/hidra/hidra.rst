Running HIDRA on NASA HEC Systems
=================================

What follows below is to compile and use Jupyterlab on Pleiades. Replace
rmalbarr with your specific Pleiades username. Note that you will have a
certain /nobackup directory number (for me it is /nobackupp12/rmalbarr).
Change this according to your number. Similarly, I have /home7/rmalbarr as my
home directory. Change this home directory number for you, accordingly.

Get set up on Pleiades
----------------------


#. Setup RSA tokens: copy from bitbucket to Pleiades Font End (pfe) and
`follow wiki <https://www.nas.nasa.gov/hecc/support/kb/enabling-your-rsa-securid-hard-token-fob_59.html>`_

#. Setup of ssh pass through,
`follow wiki <https://www.nas.nasa.gov/hecc/support/kb/setting-up-ssh-passthrough_232.html>`_

This should enable you to login to pfe with where the passcode is given by
SecureID mobile app (dual-factor authentication):

.. code-block:: bash

   ssh pfe


#. Setup the sup client
`from wiki <https://www.nas.nasa.gov/hecc/support/kb/using-the-secure-unattended-proxy-sup_145.html>`_

This will enable you to send large files between remote and local servers with
``shiftc``

From local machine, for example, run:

.. code-block:: bash

   sup shiftc <files> rmalbar@pfe:/nobackupp12/rmalbarr


#. Clone repo to your nobackup or home directory on pfe. Here, for example,
check out the hidra branch of kaiju. From pfe prompt, run:

.. code-block:: bash

   install home-brew
   git lsf install

   git branch
   git checkout hidra

From
`kaiju wiki <https://bitbucket.org/aplkaiju/kaiju/wiki/quickStart/prerequisites>`_\ , run:

.. code-block:: bash

   git clone https://YOUR_BITBUCKET_USERNAME@bitbucket.org/aplkaiju/kaiju.git
   export KAIJUHOME=$HOME/kaiju

where $HOME, for me, is set to /nobackupp12/rmalbarr. Change this according
to your username and desired directory. Note: Use Atlassian App password for
the clone command above. This is found in your bitbucket profile.

Running (e.g. HIDRA) on Pleiades
--------------------------------

#. Make a build directory within kaiju. For example, for me, it is
/nobackupp12/rmalbarr/kaiju/build. From this build directory run (for hidra
example):

.. code-block:: bash

   module load pkgsrc/2021Q2
   module load comp-intel/2020.4.304
   module load mpi-hpe/mpt.2.25
   module load szip/2.1.1
   module load hdf5/1.8.18_mpt
   FC=ifort FFLAGS=“-mkl” cmake -DENABLE_MPI=ON $HOME
   make hidra.x

Send this executable file to your /nobackup for queue submission. Note:
nobackup has much more storage space than home directory. Run qsub to submit a
pbs script as usual. A sample pbs script that I use on pfe is below:

.. code-block:: bash

   #!/bin/bash
   #Example Gamera PBS script
   #PBS -A rmalbarr
   #PBS -N HIDRAN-001_ML1N_ssOpO2_nHem
   #PBS -j oe
   #PBS -q normal
   #PBS -l walltime=8:00:00
   #PBS -l select=8:ncpus=28:mpiprocs=4:ompthreads=7:model=bro

   export EXE="./hidra_ML1N_ssOpO2.x"
   export RUNID=hidraN-001_ML1N_ssOpO2_nHem

   # source ~/.bashrc

   export KAIJUHOME=/nobackupp12/rmalbarr/kaiju/

   if [[ -z "/nobackupp12/rmalbarr/kaiju/" ]]; then
       # $KAIJUHOME environment variable is not set
       echo "The KAIJUHOME environment variable is not set"
       echo "You must either pass your environment with the -V option or"
       echo "  execute the kaiju/scripts/setupEnvironment script in your ~/.bashrc
   file"
       exit
   fi

   if [[ ! -z "$MODULE_LIST" ]]; then
       # user passed a list of modules to load as the environment variable MODULE_LLIST
       # call this with the flag " -v MODULE_LIST="<modules>" " to use this option
       # where <modules> is a space-separated list of modules in quotes
       # Example:
       #  qsub -v MODULE_LIST="intel/2021.2 ncarenv/1.3 ncarcompilers/0.5.0 mpt/2.222" RunMpi.pbs
       module purge
       module load $MODULE_LIST
   elif [[ ! -z "$MODULE_SET" ]]; then
       # user passed a module set name to load as the environment variable MODULE_SSET
       # call this with the flag "-v MODULE_SET=<set name>" to use this option
       # where <set_name> is a saved set of modules, as printed by "module savelistt"
       # Example:
       # qsub -v MODULE_SET=kaiju21 RunMpi.pbs
       module purge
       module restore $MODULE_SET
   else
       # user did not pass a module set, load a default set
       module purge
       module load pkgsrc/2021Q2
       module load comp-intel/2020.4.304
       module load mpi-hpe/mpt.2.25
       module load szip/2.1.1
       module load hdf5/1.8.18_mpt
   fi

   module list
   hostname
   date

   #module load arm-forge/21.1.1
   module load arm-forge/20.2

   export KMP_STACKSIZE=128M
   export MPI_TYPE_DEPTH=32
   export MPI_IB_CONGESTED=0
   #export OMP_NUM_THREADS=9
   export NODEFILE=$TMPDIR/nodefile.$PBS_JOBID
   cp $PBS_NODEFILE $NODEFILE

   #ddt --connect mpiexec_mpt ${KAIJUHOME}/scripts/preproc/correctOMPenvironment.shh
    ${NODEFILE} omplace ${EXE} hidra-001.xml > ${RUNID}.out
   mpiexec_mpt /nobackupp12/rmalbarr/kaiju/scripts/preproc/correctOMPenvironment.shh
    ${NODEFILE} omplace ${EXE} hidra-001.xml > ${RUNID}.out

   date

   qsub runHIDRA-002_ML1N_ssOpO2_nHem.pbs

Note: there are different systems to run on. Above, I run on Broadwell nodes
(model=bro). For Pleiades Broadwell nodes (28 cores per node), I consider 2
options:

.. code-block:: bash

   #PBS -l select=16:ncpus=28:mpiprocs=2:ompthreads=14:model=bro
   #PBS -l select=8:ncpus=28:mpiprocs=4:ompthreads=7:model=bro

Change these values and wall-time, etc, as needed.

#. Once run is complete, send all large output data files to Lou mass storage
(lfe). For example, for me to send ``<files>`` to my lfe:

.. code-block:: bash

   shiftc <files> rmalbarr@lfe:/u/rmalbarr/

Check to see where they go. From pfe prompt:

.. code-block:: bash

   ssh lfe

To exit back to pfe prompt, run: ``exit``

Note: I routinely reach disk quota on /nobackup so its best to send all
large output files to lfe. For Jupyter analysis, copy individual files back to
/nobackup on pfe.
