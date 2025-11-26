Steady-State Run
================

.. These instructions are to build the MPI gamhelio on Cheyenne supercomputer.

.. First, you must create a directory in which to compile the ``kaiju`` code for gamhelio case. Typically, this is ``$KAIJUHOME/build_helio``\ , where ``KAIJUHOME`` is the path to the directory created when you cloned the kaiju repository. Now you can build the ``kaiju`` software:

.. .. code-block:: shell

..    cd $KAIJUHOME
..    mkdir build_helio
..    cd build_helio
..    cmake  -DENABLE_MPI=ON  ..
..    make -j4 gamhelio_mpi

.. When the build is complete, you will find the compiled executable in the bin subdirectory of your build directory:

.. .. code-block:: shell

..    bash-4.2$ ls bin
..    gamhelio_mpi.x

.. Creating Grid and Boundary Conditions
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. The default grid for the inner heliosphere simulations is spherical evenly spaced in r, theta, phi directions. For now, the south and north poles are cut out, so theta is within the range [0.1pi,0.9pi]. Grid parameters such as limits in r and theta directions and the number of cells Ni, Nj, Nk are set in ``KAIJUHOME/kaipy/gamhelio/ConfigScripts/startup.config`` under ``[Grid]``. The default values for these parameters are below.

.. Other options for the inner helio spherical grid exist in ``KAIJUHOME/kaipy/gamera/gamGrids.py``. GenKSph creates a default uniform spherical grid. Other options are GenKSphNonU which makes a non-uniform grid in r-direction changing smoothly from finer grid near the inner boundary to coarser grid near the outer boundary; and GenKSphNonUG which creates a custom grid for a CME project with a fine uniform grid in the region 0.1-0.3 AU and a non-uniform coarser grid further out to 1 AU. If needed modify the preprocessing script wsa2gamera.py to use any of these options or create your own grid function in ``KAIJUHOME/kaipy/gamera/gamGrids.py``.

.. To generate h5 files with grid and boundary conditions you need to have ready a config file ``KAIJUHOME/kaipy/gamhelio/ConfigScripts/startup.config``. If needed, modify paths to the output helio grid file, output file innerbc.h5, and input WSA fits file. 

.. .. code-block:: shell

.. ..   [Gamera]
..    gameraGridFile = heliogrid.h5
..    GridDir = ./
..    gameraIbcFile = innerbc.h5
..    IbcDir = ./

.. ..   [Grid]
..    tMin = 0.1
..    tMax = 0.9
..    Rin  = 21.5
..    Rout = 220.
..    Ni   = 128
..    Nj   = 64
..    Nk   = 128

.. ..   [WSA]
..    ;wsafile is the path to the WSA fits file relative to $KAIJUHOME
..    ;Helio test uses WSA file for Carrington Rotation 2193, by default
..    wsafile = examples/helio/vel_201708132000R002_ahmi.fits
..    density_temperature_infile = no
..    gauss_smooth_width         = 0 ; 8
..    normalized                 = no

.. ..   [Constants]
..    gamma  = 1.5
..    Nghost = 4
..    Tsolar = 25.38
..    nCS = 1100. ; in [cm-3]
..    TCS = 1.e6  ; in [K]

.. ..   [Normalization]
..    B0 = 1.e-3  ; in [Gs] equals to 100 [nT]
..    n0 = 200.   ; in [cm-3]

.. Now run   

.. .. code-block:: shell

..    wsa2gamera.py ~/kaiju/kaipy/gamhelio/ConfigScripts/startup.config

.. Check that you successfully produced files heliogrid.h5 and innerbc.h5.

.. .. code-block:: shell

..    user@cheyenne5:/rundir/> h5dump -n innerbc.h5 
..    HDF5 "innerbc.h5" {
..    FILE_CONTENTS {
..     group      /
..     dataset    /br
..     dataset    /br_kface
..     dataset    /rho
..     dataset    /temp
..     dataset    /vr
..     dataset    /vr_kface
..     }
..    }

..    user@cheyenne5:/rundir> h5dump -n heliogrid.h5
..    HDF5 "heliogrid.h5" {
..    FILE_CONTENTS {
..     group      /
..     dataset    /X
..     dataset    /Y
..     dataset    /Z
..     }
..    }

.. XML input file
.. ^^^^^^^^^^^^^^

.. An example xml input file for mpi gamera helio run is as follows:

.. .. code-block:: shell

..    <?xml version="1.0"?>
..    <Kaiju>
..    <Gamera>
..        <sim runid="wsa" doH5g="T" H5Grid="heliogrid.h5" icType="user" pdmb="1.0" rmeth="7UP"/>
..        <time tFin="850."/>
..        <output dtOut="12." tsOut="50" timer="F"/>
..        <physics doMHD="T" gamma="1.5"/>
..        <prob Tsolar = "25.38"/>
..        <restart resFile = "wsa.Res.00008.h5" dtRes="1000." doRes="F"/>
..        <iPdir N="4" bcPeriodic="F"/>
..        <jPdir N="2" bcPeriodic="F"/>
..        <kPdir N="4" bcPeriodic="T"/>
..    </Gamera>
..    </Kaiju>

.. For high-resolution run 1024x512x1024 use this decomposition

.. .. code-block:: shell

..    <iPdir N="8" bcPeriodic="F"/>
..    <jPdir N="4" bcPeriodic="F"/>
..    <kPdir N="8" bcPeriodic="T"/>

.. PBS script
.. ^^^^^^^^^^

.. Here is an example pbs script to run mpi gamera (for low resolution helio run 256x128x256)

.. .. code-block:: shell

..    #!/bin/bash
..    #PBS -A ujhb0015
..    #PBS -N heliompi
..    #PBS -j oe
..    #PBS -q regular
..    #PBS -l walltime=00:20:00
..    #PBS -l select=16:ncpus=36:mpiprocs=2:ompthreads=18

..    #Example usage

..    export TMPDIR=/glade/scratch/$USER/temp
..    mkdir -p $TMPDIR

..    export EXE="./gamhelio_mpi.x"
..    export RUNID="wsa"

..    #Optional stuff to load an environment
..    source ~/.bash_profile

..    if [[ -z "$KAIJUHOME" ]]; then
..        # $KAIJUHOME environment variable is not set
..        echo "The KAIJUHOME environment variable is not set"
..        echo "You must either pass your environment with the -V option or"
..        echo "  execute the kaiju/scripts/setupEnvironment script in your ~/.bashrc file"
..        exit
..    fi
..    if [[ ! -z "$MODULE_LIST" ]]; then
..        # user passed a list of modules to load as the environment variable MODULE_LIST
..        # call this with the flag '-v MODULE_LIST="<modules>"' to use this option
..        # where <modules> is a space-separated list of modules in quotes
..        # Example:
..        #  qsub -v MODULE_LIST="intel/2021.2 ncarenv/1.3 ncarcompilers/0.5.0 mpt/2.22" RunMpi.pbs
..        module purge
..        module load $MODULE_LIST
..    elif [[ ! -z "$MODULE_SET" ]]; then
..        # user passed a module set name to load as the environment variable MODULE_SET
..        # call this with the flag '-v MODULE_SET=<set name>' to use this option
..        # where <set_name> is a saved set of modules, as printed by 'module savelist'
..        # Example:
..        # qsub -v MODULE_SET=kaiju21 RunMpi.pbs
..        module purge
..        module restore $MODULE_SET
..    else
..        # user did not pass a module set, load a default set
..        module purge
..        module restore mpikaiju
..    fi

..    if [[ ! -z "$MPT_VERSION" ]]; then
..        echo "USING MPIEXEC_MPT"
..        export MPI_TYPE_DEPTH=32
..        export OMP_NUM_THREADS=36
..        export MPI_IB_CONGESTED=0
..        export NODEFILE=$TMPDIR/nodefile.$PBS_JOBID
..        cp $PBS_NODEFILE $NODEFILE
..        export MPICOMMAND="mpiexec_mpt $KAIJUHOME/scripts/preproc/correctOMPenvironment.sh $NODEFILE omplace"
..    else
..        echo "USING MPIRUN"
..        export MPICOMMAND="mpirun"
..        export OMP_NUM_THREADS=18
..        export I_MPI_PIN_DOMAIN="omp"
..    fi


..    module list
..    hostname
..    date
..    export KMP_STACKSIZE=128M
..    export JNUM=${PBS_ARRAY_INDEX:-0}
..    echo "Running $EXE"
..    ${MPICOMMAND} ${EXE} ${RUNID}.xml ${JNUM} >> ${RUNID}.${JNUM}.out
..    date

.. The example above uses 16 computer nodes (2 MPI ranks per node) creating 32 processes for 32 MPI ranks (4x2x4 = 32 in decomposition for low resolution run above). 

.. For high resolution run 1024x512x1024 we have 8x4x8 = 256 MPI ranks so we select 128 nodes (with 2 MPI ranks/node).

.. .. code-block:: shell

..    #PBS -l walltime=11:59:00
..    #PBS -l select=128:ncpus=36:mpiprocs=2:ompthreads=36:mem=109GB

.. See PBS job basics `here <https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs>`_ on Cheyenne.

.. Using MPT
.. ^^^^^^^^^

.. If you are running with the MPT mpi library, the submission script will require some additional modifications, described in a dedicated wiki page [[Running with MPT]].

.. Submitting a run
.. ^^^^^^^^^^^^^^^^

.. Create a run directory. Link files with the grid heliogrid.h5, inner boundary conditions innerbc.h5 and gamera executible gamera_mpi.x, for example

.. .. code-block:: shell

..    ln -s ~/kaiju/build/bin/gamera_mpi.x gamera_mpi.x

.. Place into run directory pbs script gamera.pbs and input xml file wsa.xml. So here is a content of  your run directory

.. .. code-block:: shell

..    user@cheyenne5:/glade/work/user/helioRun> ls
..    gamera_mpi.x  gamera.pbs  heliogrid.h5  innerbc.h5  wsa.xml

.. Run the job 

.. .. code-block:: shell

..    qsub gamera.pbs

.. Check a status of your job in a queue

.. .. code-block:: shell

..    qstat -u username

.. Normalization in Gamera-Helio
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. The three main normalization parameters are


.. #. Length L = 1R_S = 6.955e10 cm
.. #. Magnetic field magnitude B0 = 100 nT = 1.e-3 Gs
.. #. Number density n0 = 200 cm-3

.. Velocity is normalized to the Alfven velocity V0 = B0/sqrt(4 pi rho0) ~ 154 km/s.
.. Time is normalized to t = L/V0 = 4637 s ~ 1.25 h
.. Pressure is normalized to the magnetic pressure B0^2/(4*pi).

.. Output data
.. ^^^^^^^^^^^
