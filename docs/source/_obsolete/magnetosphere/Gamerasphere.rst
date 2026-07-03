Magnetosphere Quick Start
=========================

This is a quick set of instructions to run a coupled Gamera-ReMIX
(magnetosphere-ionosphere) run on Cheyenne.  It uses the executable
"voltron.x".

OMP Magnetosphere
-----------------

Basic magnetosphere runs require five things


* Grid file: Collection of grid corners in HDF5 format
* Solar wind: Solar wind time series in HDF5 format
* Input deck: Run configuration parameters in XML format
* Executable: Duh (voltron.x)
* Run script: Job script to setup run for a single node on Cheyenne

Note, this is assuming the main repo has been downloaded to a base directory
kaiju

Grid file
---------

The simplest way of generating is to use kaiju/scripts/genLFM.py to create a
standard LFM-style grid of resolution D/Q/O/H for double/quad/oct/hex.  Note,
you'll likely have to do this on Casper, the post-processing node.

.. code-block::

   [skareem@casper26:~/kaiju/scripts]$ genLFM.py -gid D
   Generating Double LFM-style grid ...

   Output: lfmD.h5
   Size: (48,48,64)
   Inner Radius: 2.000000
   Sunward Outer Radius: 28.000106
   Tail Outer Radius: 301.011944
   Low-lat BC: 45.000000
   Ring params:
   <ring gid="lfm" doRing="T" Nr="4" Nc1="8" Nc2="16" Nc3="32" Nc4="32"/>

   Writing to lfmD.h5

The grid generation spits parameters for ring average which will be added to
the XML input deck.

Solar wind
----------

Gamera reads solar wind data using HDF5. The simplest way to generate an
appropriate HDF5 file is to use the omni2wind script, which downloads data
from the CDAS database between the specified timeframe and converts it into
the correct format to be read in by Gamera. omni2wind.py utilizes CDASWS as an
interface to the CDAS datasets as well as geopack for coordinate
transformations. The modules must be installed before the script can be run.
ai.cdas can be installed using the command:

.. code-block::

   pip install cdasws

The library has a dependency on Spacepy which also needs to be installed
(https://spacepy.github.io/install.html). The instructions to install geopack
are located in the README file in $KAIJUDIR/external/geopack-2008. You will
need to specify the compiler, and must configure to use 8-byte reals. Here is
an example installation command on an intel compiler:

.. code-block::

   python setup-geopack.py install config_fc --fcompiler=intelem --f77flags=-r8

To install the modules on Cheyenne, you will need to create your own python
library environment by cloning the NPL library and then activating it with the
commands:

.. code-block::

   ncar_pylib -c [version of library to be cloned] [location to clone library into]
   source $CLONEDLIBRARY/bin/activate

The cdasws and geopack modules can then be installed into the new library and
will be included in your path when the cloned library is activated.
omni2wind.py generates three files: bcwind.h5, the HDF5 solar wind input file
for Gamera, as well as two images, the fit to Bx used and the solar wind input
variables. Example command line output can be seen below. The time-period
selected must be larger that an hour and a half, otherwise it will not produce
a bcwind.h5 file and end in an error.

.. code-block::

   (NPL) bash-4.2$ omni2wind.py -t0 2010-02-25T12:00:00 -t1 2010-02-25T14:00:00
   Retrieving f10.7 data from CDAWeb
   100% [................................................................................] 4668 / 4668
   Average f10.7:  81.1
   Retrieving solar wind data from CDAWeb
   100% [..............................................................................] 27878 / 27878

     RECALC_08: RADIAL SOLAR WIND --&gt; GSW SYSTEM IDENTICAL HERE
     TO STANDARD GSM (I.E., XGSW AXIS COINCIDES WITH EARTH-SUN LINE)

   Bx Fit Coefficients are  [-3.3785629   0.60784241  0.11600477]
   Saving "OMNI_HRO_1MIN.txt_bxFit.png"
   /glade/u/home/adamm/dav_pylib/lib/python3.6/site-packages/matplotlib/cbook/deprecation.py:107: MatplotlibDeprecationWarning: Adding an axes using the same arguments as a previous axes currently reuses the earlier instance.  In a future version, a new instance will always be created and returned.  Meanwhile, this warning can be suppressed, and the future behavior ensured, by passing a unique label to each axes instance.
     warnings.warn(message, mplDeprecation, stacklevel=1)
   Saving "OMNI_HRO_1MIN.txt.png"
   Converting to Gamera solar wind file
       Found 15 variables and 120 lines
       Offsetting from LFM start ( 0.00 min) to Gamera start ( 0.00 min)
   Writing Gamera solar wind to bcwind.h5

Notes
-----

An HDF5 file can also be created from an existing LFM-style solar wind input
using the reWind2.py script. For an example see
kaiju/examples/earthcmi/OMNI_HRO_1MIN_16772.txt_SW-SM-DAT.

Note the GAMERA input of solar wind file requires more information than the
density, velocity and IMF in LFM. The extra variables include F10.7,
diple tilt, MJD (Julian date of each record), fitting coefficients of IMF
(Bx0, ByC, BzC). Sound speed is not taken but the temperature is needed.

The MJD is especially helpful when a different time step of solar wind input
is needed.

The default cadence of solar wind input is 1 min. In order to use higher rate
solar wind input, we used to have to modify the init-fortran code in LFM. Here
in the GAMERA world, the solar wind file MJD is the only variable to be
modified for a flexible input cadence.

.. code-block::

   [skareem@casper26:~]$ reWind2.py OMNI_HRO_1MIN_16772.txt_SW-SM-DAT
   Reading LFM solar wind from OMNI_HRO_1MIN_16772.txt_SW-SM-DAT
       Found 11 variables and 1440 lines
       Offsetting from LFM start ( 0.00 min) to Gamera start ( 0.00 min)
   Writing Gamera solar wind to bcwind.h5

Input deck
----------

The input deck is an XML file that specifies various parameters to be read at
runtime.  An example for a double resolution run is shown below.  Note that
times used to be input in code units and the conversion is 63.8 seconds. Now
we have changed the time variables in xml to be in seconds. This means

.. code-block::

   #!xml
   <?xml version="1.0"?>
   <!-- Magnetosphere params, Voltron times in seconds -->
   <!-- MJD0 is modified Julian date of T=0 in solar wind file -->
   <VOLTRON>
     <time tFin="36000.0"/>
     <output dtOut="60.0" tsOut="100" doTimer="F"/>
     <restart dtRes="1800.0"/>
     <coupling dt="5.0"/>
   </VOLTRON>
   <Gamera>
     <sim runid="msphere" doH5g="T" H5Grid="lfmD.h5" icType="user" pdmb="1.0" pFloor="1.0e-8" dFloor="1.0e-6" rmeth="7UP"/>
     <restart doRes="F" resId="msphere" nRes="0"/>
     <physics doMHD="T" doBoris="T" Ca="10.0"/>
     <prob Rho0="0.2" P0="0.001"/>
     <ring gid="lfm" doRing="T" Nr="4" Nc1="8" Nc2="16" Nc3="32" Nc4="32"/>
     <wind tsfile="bcwind.h5"/>
   </Gamera>
   <!-- Remix params -->
   <REMIX>
     <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
     <conductance pedmin="2.0" hallmin="1.0" sigma_ratio="3.0" const_sigma="True" ped0="5.0"/>
   </REMIX>

With these parameters the run will go for 10 hours (36000 s), outputting every
1 minute (dtOut=60.0 s) and writing restarts every 30 minutes
(dtRes=1800.0 s). The ReMIX-Gamera coupling is done every 5 seconds
(coupling dt=5.0). Note the coupling dt is a required input for the
ionospheric coupling.

The runid gives the name of the output mhd/mix file. H5Grid is the name of the
grid file that has to be present in the run directory. If the job to be
submitted with this xml file is not restarting from a previously saved state,
set doRes as False. Otherwise, to restart, set doRes as True and make sure
the resFile is linked to the right restarting file. As of 23 April 2020
restarts are specified with ID/# instead of filename. Instead of
restart/resFile, specify restart/resID and restart/nRes. The restart file
msphere.Res.00005.h5 would be:

``<restart resId="msphere" nRes="5"/>``

Specifying nRes="-1" will read the XXXXX symbolic link.

``<physics>`` domain:

``<prob>`` domain:

``<ring>`` domain: gid tells the model which type of grid is being
used. It is supposed to be consistent with the input grid file. Ring average
technique is implemented if doRing is True. The number of parameters for ring
average is set in Nr, and each parameter is listed as Nc1, Nc2, ... Usually
there are four parameters for double resolution grid and 8 parameters for quad
resolution grid.

``<wind>`` domain: tsfile takes the name of the solar wind file to
be used.

Under the REMIX block, the Np and Nt gives the number of grid cells in Remix
along longitude and latitude. The low latitude boundary is set at 45 deg
latitude.

This setup uses constant conductance, an example of using the conductance
model would replace the conductance block with this,

.. code-block::

   #!xml
   <conductance F107="100.0" pedmin="2.0" hallmin="1.0" sigma_ratio="3.0" const_sigma="False" ped0="10.0"/>

Note, that starting the initial state (a pure dipole) with the conductance
model turned on can sometimes be erratic.  It's better to spin up for a while
using constant conductance then do a restart where the conductance model is
turned on.  But whatever, you do you.

Executable and run script
-------------------------

Compile the target "voltron.x" and then just run it with a PBS script.  You
should probably request a full node and run with 72 threads.  Below is the PBS
script that I use, but you might need to make changes to it for your setup.
This is assuming an XML file called "cmiD.xml", that you have a preset module
savelist, and that voltron.x is in your path.  Probably some other stuff too.
I ran it with

.. code-block::

   qsub -v CHIMPEXE="voltron.x" -N cmiD RunCHIMP.pbs

and it worked.

.. code-block::

   #!bash
   #!/bin/bash
   #PBS -A P28100045
   #PBS -N KAIJU
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72

   #Example usage
   #qsub -v CHIMPEXE="slice.x" -N cmiD RunCHIMP.pbs
   #Module savelist "kaiju"
   #Currently Loaded Modules:
   #  1) git/2.9.5    (H)   4) impi/2018.4.274       7) python/2.7.16
   #  2) intel/18.0.5       5) ncarenv/1.3           8) cmake/3.14.4
   #  3) hdf5/1.10.5        6) ncarcompilers/0.5.0

   export EXE=${CHIMPEXE:-"slice.x"}
   export RUNID=${PBS_JOBNAME}

   source ~/.bashrc
   module restore kaiju

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date

Generating solar wind file manually
-----------------------------------

Solar wind data can be defined using an HDF-5 file.  This is specified in the
XML input deck under "Gamera/wind/tsfile" as seen in the example XML file
above.

Solar wind HDF5 file input units

* Time [s]
* Density [#/cc]
* Velocity [m/s]
* Pressure [nPa]
* Magnetic field [nT]

Gamera will read the HDF file and convert to its internal code units (see
normalization above).  Below is an example of a Python script to generate a
solar wind file.

.. code-block::

   #!python

   import h5py
   import numpy as np
   import sys

   TW = 1.0e+5     #Default temperature, K
   nW = 5          #Default density, #/cc
   VxW = 400.0     #Default wind, km/s
   f107val = 100.0 #Default f10.7 flux
   tilt = 0.0      #Default dipole tilt, radians
   mjd0 = 58767.0  #Default MJD, set for 2019-10-11 00:00:00

   Bx0 = 0.0 #Default Bx offset for planar front, keep at zero
   ByC = 0.0 #Default By coefficient used to calculate Bx, include if want tilted field
   BzC = 0.0 #Default Bz coefficient used to calculate Bx, include if want tilted field

   fOut = "bcwind.h5"

   #Time bounds [hours]
   tMin = 0.0
   tMax = 6.0
   dt = 60.0 #Cadence [s]

   SimT = (tMax-tMin)*60.0*60.0
   NumT = np.int( np.ceil(SimT/dt)+1 )

   print("Generating %d slices, T=[%5.2f,%5.2f]"%(NumT,tMin,tMax))

   T = np.linspace(tMin,tMax,NumT)
   D = np.zeros(NumT)
   Temp = np.zeros(NumT)
   Vx = np.zeros(NumT)
   Vy = np.zeros(NumT)
   Vz = np.zeros(NumT)
   Bx = np.zeros(NumT)
   By = np.zeros(NumT)
   Bz = np.zeros(NumT)
   f107 = np.zeros(NumT)
   ThT = np.zeros(NumT)
   mjd = np.zeros(NumT)
   symh = np.zeros(NumT)

   tWin = 1.0 #Window times [hr]
   for i in range(NumT):
       t = T[i] #Time in hours
       if (t <= tWin):
           D[i] = nW
           Vx[i] = -VxW
           Temp[i] = TW
           f107[i] = f107val
           ThT[i] = tilt
           mjd[i] = mjd0 + T[i]/24.0
       elif (t <= 3*tWin):
           D[i] = nW
           Vx[i] = -VxW
           Temp[i] = TW
           Bz[i] = -5.0
           f107[i] = f107val
           ThT[i] = tilt
           mjd[i] = mjd0 + T[i]/24.0
       elif (t <= 6.0*tWin):
           D[i] = nW
           Vx[i] = -VxW
           Temp[i] = TW
           Bz[i] = +5.0
           f107[i] = f107val
           ThT[i] = tilt
           mjd[i] = mjd0 + T[i]/24.0
       else:
           D[i] = nW
           Vx[i] = -VxW
           Temp[i] = TW
           Bz[i] = -5.0
           f107[i] = f107val
           ThT[i] = tilt
           mjd[i] = mjd0 + T[i]/24.0

   #Write solar wind
   #t,D,V,Temp,B = [s],[#/cm3],[m/s],[K],[nT]

   oTScl = (60*60.0) #hr->s
   oDScl = 1.0
   oVScl = 1.0e+3 #km/s->m/s
   oTempScl = 1.0
   oBScl = 1.0


   with h5py.File(fOut,'w') as hf:
       hf.create_dataset("T" ,data=oTScl*T)
       hf.create_dataset("symh" ,data=symh)
       hf.create_dataset("D" ,data=oDScl*D)
       hf.create_dataset("Temp" ,data=oTempScl*Temp)
       hf.create_dataset("Vx",data=oVScl*Vx)
       hf.create_dataset("Vy",data=oVScl*Vy)
       hf.create_dataset("Vz",data=oVScl*Vz)
       hf.create_dataset("Bx",data=oBScl*Bx)
       hf.create_dataset("By",data=oBScl*By)
       hf.create_dataset("Bz",data=oBScl*Bz)
       hf.create_dataset("tilt",data=ThT)
       hf.create_dataset("f10.7",data=f107)
       hf.create_dataset("MJD",data=mjd)
       hf.create_dataset("Bx0",data=Bx0)
       hf.create_dataset("ByC",data=ByC)
       hf.create_dataset("BzC",data=BzC)

Visualization
-------------

Data can be read into VisIt or Paraview using the "kaiju/scripts/genXDMF.py"
script which will create an XDMF file which either of those viewers can
natively read.

Alternatively, you can try "kaiju/scripts/msphpic.py" script which uses the
Python post-processing routines to generate a multi-panel plot of a
magnetosphere run.  It uses Python libraries that you may or may not have,
I don't know your life.

Creating gamera videos in parallel via slurm job submission and gamsphVid.py.

.. code-block::

   #!/bin/bash -l
   #SBATCH -J VidMHD
   #SBATCH --output=%x_%A_%a.out
   #SBATCH -t 24:00:00
   #SBATCH -A UJHB0010
   #SBATCH -p dav
   #SBATCH --array=1-36
   #Defaults
   export BASE="/glade/u/home/skareem/Work/gamercm/Data/"
   export RUNID="vapD9"

   export TS="60"
   export TE="1480"
   export DT="60"

   export ARG=""
   export ODIR="VidMHD"
   export JNUM=${SLURM_ARRAY_TASK_COUNT:-"1"}
   export JID=${SLURM_ARRAY_TASK_ID:-"1"}

   echo "My JobID: " $SLURM_ARRAY_TASK_ID " of " $SLURM_ARRAY_TASK_COUNT
   export IDIR="${BASE}${RUNID}"
   source ~/.bashrc
   setdav
   gamsphVid.py -d $IDIR -ts $TS -te $TE -dt $DT -Nblk $JNUM -nID $JID -o $ODIR $ARG

Modify and launch the above script using "sbatch ABOVEFILE.s"

.. code-block::

   #sbatch array#

Specifies how many parallel jobs you want (36 here)

.. code-block::

   export TS="60"
   export TE="1480"
   export DT="160"

are the time domain in which you want to create movies.  TS is the start time
(in minutes) of your run starting from T=0.  TE is time end (in minutes). DT
is the time (in seconds) between slices.  DT will attempt to find the cloest
timestep, so if outputs are 60s but you input DT=15s, then you get duplicated
frames.

.. code-block::

   source ~/.bashrc
   setdav

Set up your environment variables.  These are Kareem specific, so modify to
match whatever your Casper environment is.

MPI Magnetosphere
-----------------

That's still like a whole thing.
