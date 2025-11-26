CHIMP Quick Start Guide
=======================

Introduction
------------

This page provides instructions to help the user get up and running quickly
with `CHIMP <https://cgs.jhuapl.edu/Models/chimp.php>`_. We provide
instructions for using the CHIMP tool ``push.x``, which computes the
trajectories of test particles within the domain of the MAGE model run.

These instructions assume you are running on the ``derecho`` supercomputer,
and that you have already made a successful ``kaiju`` model run in your
current working directory. The model run uses the run identifier ``geospace``.

Running ``push.x``
------------------

Create the XML file
^^^^^^^^^^^^^^^^^^^

Create an XML file describing the input parameters for ``push.x``. It should
look something like this on ``derecho``:

.. code-block:: xml

   <?xml version="1.0" ?>
   <Kaiju>
       <Chimp>
         <sim runid="geospace"/>
         <time T0="36000.0" dt="0.5" tFin="43200.0"/>
         <fields doMHD="F" ebfile="geospace" grType="LFM" isMPI="T"/>
         <parallel Ri="4" Rj="4" Rk="1"/>
         <pusher epsht="0.05" imeth="FO"/>
         <tps Np="2500" species="Op"/>
         <units uid="LFM"/>
         <output doEQProj="T" dtOut="5.0" tsOut="100"/>
         <radius max="2.25" min="2.05"/>
         <phi max="360.0" min="0.0"/>
         <alpha max="180.0" min="90.0"/>
         <height max="90.0" min="15.0" doOutflow="T"/>
         <energy doEBInit="F" doLog="F" max="10.0" min="0.001"/>
         <domain rmax="25.0" rmin="2.005"/>
         <tracer epsds="0.05"/>
         <stream doStream="F" max="39600.0" min="36000.0"/>
      </Chimp>
   </Kaiju>

``T0`` is the start time for the test particle trajectories, in seconds from
the start time of the simulation being used for the calculation.

``tFin`` is the stop time the test particle trajectories, in seconds from
the start time of the simulation being used for the calculation.

For other parameters see the :doc:`Chimp XML <chimpXML>` page.

Create the PBS file
^^^^^^^^^^^^^^^^^^^

The ``push.x`` tool is usually run in parallel (using a PBS "job array") to
improve performance. The PBS script ``RunCHIMP.pbs`` might look something like
this:

.. code-block:: bash

   #PBS -A P28100045
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=128:ompthreads=128
   #PBS -m ae

   export RUNID=${PBS_JOBNAME}

   module list
   hostname
   date
   export OMP_NUM_THREADS=128
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   EXE=./push.x
   echo "Running $EXE"
   ${EXE} ${RUNID}.xml ${JNUM} >& ${RUNID}.${JNUM}.out
   date

Run ``push.x``
^^^^^^^^^^^^^^

Copy the executable for ``push.x`` into your working directory. For example:

.. code-block:: bash

   cp $KAIJUDIR/build/bin/push.x .

Here is an example command to submit two parallel jobs:

.. code-block:: bash

   qsub -v -J 1-2 -N geospace RunCHIMP.pbs

``-J 1-2`` means the same ``push.x`` run is split into 2 jobs that run in
parallel on separate nodes.

``-N geospace`` gives the job array the same name as the run ID
(``geospace``). It will also force ``push.x`` to use the XML file called
``geospace.xml`` for its input.

The job index in the job array is used to seed the generation of test
particles within each run. This allows the code to be reproducible. If you
want to run additional test particles with the same xml file, you need to make
sure the job numbers are not the same i.e. if your first set of runs was
submitted with the command above, you can run two additional, unique jobs with
the command:

.. code-block:: bash

   qsub -c -J 3-4 -N geospace RunCHIMP.pbs

CHIMP output format
^^^^^^^^^^^^^^^^^^^

The output file from a CHIMP run will be look something like
``runid.000000.h5part`` and contain the variables:

``x``, ``y``, ``z``
   Particle position coordinates

``K``
   particle energy [keV]

``isIn``
   ``0`` or ``1`` depending on whether the particle is alive or not (left the
   simulation domain) and whether it's been "born" yet or not.

``id``
   Particle ID number

``xeq``, ``yeq``, ``Teq``
   Position where the particle last crossed the Z=0 plane, and at what time or
   instantaneous projection depending on ``doTrc`` in the XML file.

``Keq``/``ebKeq``
   Energy at the equator in lab and ExB frame

``Aeq``
   Equatorial pitch angle

``OCb``
   Set to ``0``/``1``/``2`` based on the topology of the field line at the
   particle position: ``0`` = IMF, ``1`` = open, ``2`` = closed. This will be
   all zeroes if you didn't tell it to do field-line tracing (computationally
   expensive) in the XML file.

``T0p``
   The "birthday" of the particle, when it was/will be "born".
