
GAMERA-REMIX
============

A GAMERA-REMIX (GR) run without RCM simulates the magnetosphere without ring
current but coupled with the ionosphere represented by the REMIX module. Note
we still need to use a Voltron executable instead of Gamera executable for a
GR run. To turn off RCM, set DtDeep = -1 and Gamera/dosrc to F.

Below are an example xml file of model parameters and an example pbs file of
running parameters.

Example xml for a double-resolution MPI run
-------------------------------------------

.. code-block:: xml

   <KAIJU>
   <?xml version="1.0"?>
   <VOLTRON>
     <time tFin="43210.0"/>
     <spinup doSpin="T" tSpin="3600.0" tIO="0.0"/>
     <output dtOut="300.0" tsOut="100">
     <coupling dt="5.0"/>
     <restart dtRes="1800.0"/>
     <threading NumTh="36"/>
   </VOLTRON>
   <Gamera>
     <sim runid="msphere" doH5g="T" H5Grid="lfmD.h5" icType="user" pdmb="0.75" rmeth="7UP"/>
     <floors dFloor="1.0e-6" pFloor="1.0e-6"/>
     <restart doRes="F" resID="msphere" nRes="-1"/>
     <physics doMHD="T" doBoris="T" Ca="10.0"/>
     <ring gid="lfm" doRing="T"/>
     <wind tsfile="bcwind.h5"/>
     <source doSource="F"/>
     <iPdir N="2" bcPeriodic="F"/>
     <jPdir N="2" bcPeriodic="F"/>
     <kPdir N="1" bcPeriodic="T"/>
     <threading NumTh="18"/>
   </Gamera>
   <REMIX>
     <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
     <conductance doStarlight="T" />
     <precipitation aurora_model_type="FEDDER" alpha="0.34" beta="9.4362323" R="0.042567956" doAuroralSmooth="F"/>
   </REMIX>
   </KAIJU>

Example pbs file for a double-resolution MPI run
------------------------------------------------

.. code-block:: bash

   #!/bin/bash
   #PBS -A P54048000
   #PBS -N 04Fe-1
   #PBS -j oe
   #PBS -q economy
   #PBS -l walltime=06:00:00
   #PBS -l select=2:ncpus=36:mpiprocs=2:ompthreads=18+1:ncpus=36:mpiprocs=1:ompthreads=36

   #Example usage

   export TMPDIR=/glade/scratch/$USER/temp
   mkdir -p $TMPDIR

   export EXE="./voltron_mpi.x"
   export RUNID="cmriD1"

   #Optional stuff to load an environment
   source ~/.bashrc

   if [[ -z "$KAIJUHOME" ]]; then
       # $KAIJUHOME environment variable is not set
       echo "The KAIJUHOME environment variable is not set"
       echo "You must either pass your environment with the -V option or"
       echo "  execute the kaiju/scripts/setupEnvironment script in your ~/.bashrc file"
       exit
   fi

   if [[ ! -z "$MODULE_LIST" ]]; then
       # user passed a list of modules to load as the environment variable MODULE_LIST
       # call this with the flag '-v MODULE_LIST="<modules>"' to use this option
       # where <modules> is a space-separated list of modules in quotes
       # Example:
       #  qsub -v MODULE_LIST="intel/2021.2 ncarenv/1.3 ncarcompilers/0.5.0 mpt/2.22" RunMpi.pbs
       module purge
       module load $MODULE_LIST
   elif [[ ! -z "$MODULE_SET" ]]; then
       # user passed a module set name to load as the environment variable MODULE_SET
       # call this with the flag '-v MODULE_SET=<set name>' to use this option
       # where <set_name> is a saved set of modules, as printed by 'module savelist'
       # Example:
       # qsub -v MODULE_SET=kaiju21 RunMpi.pbs
       module purge
       module restore $MODULE_SET
   else
       # user did not pass a module set, load a default set
       module purge
       module restore mpikaiju
   fi

   if [[ ! -z "$MPT_VERSION" ]]; then
       echo "USING MPIEXEC_MPT"
       export MPI_TYPE_DEPTH=32
       export OMP_NUM_THREADS=36
       export MPI_IB_CONGESTED=0
       export NODEFILE=$TMPDIR/nodefile.$PBS_JOBID
       cp $PBS_NODEFILE $NODEFILE
       export MPICOMMAND="mpiexec_mpt $KAIJUHOME/scripts/preproc/correctOMPenvironment.sh $NODEFILE omplace"
   else
       echo "USING MPIRUN"
       export MPICOMMAND="mpirun"
       export OMP_NUM_THREADS=18
       export I_MPI_PIN_DOMAIN="omp"
   fi


   module list
   hostname
   date
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   ${MPICOMMAND} ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date
