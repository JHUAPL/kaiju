Debugging
=========

Compile in debug mode
---------------------

These instructions are specific for debugging on Cheyenne.

The first step of debugging is to compile the code in the debug mode.
For GAMERA or Voltron, add the debug flag to cmake when compiling, e.g.:

.. code-block:: bash

   mkdir build
   cd build
   cmake -DENABLE_MPI=ON -DENABLE_MKL=OFF -DCMAKE_BUILD_TYPE=DEBUG ..
   make voltron_mpi.x


For TIEGCM, turn the debug flag to TRUE in the job file, e.g.:

.. code-block:: bash

   set debug=TRUE


Detailed explanation:

Ensure that gamera (or whichever application you're debugging) has been
built with debugging information included. For gamera, this can be done
by setting the CMAKE option CMAKE_BUILD_TYPE to either RELWITHDEBINFO
or DEBUG. RELWITHDEBINFO will build the standard optimized version of
gamera, but include some debugging information. Not all debugging
information can be included in the RELWITHDEBINFO version, so if you
find that you can't get all of the information that you want out of
Allinea, build with DEBUG instead. This is an unoptimized version of the
application with full debugging information.

Launch Alinea
-------------

Load the arm-forge module. As of 2022, the default version is
arm-forge/20.2.

.. code-block:: python

   module load arm-forge


Load "arm-forge/20.2" when compiling gamera with intel 21 compilers and
"arm-forge/19.1" for older versions of intel. Also remember to load the
same modules when launching the debugging job.

Launch the Alinea debugger with the command "ddt". If the test case is
multi-threaded or uses multiple MPI processes, it is recommended to run
in an interactive job. Select the "Run" option, and specify appropriate
settings for the application, working directory, MPI, OpenMP, etc....


.. image:: https://bitbucket.org/repo/kMoBzBp/images/1029093052-allineaDDT.PNG
   :target: https://bitbucket.org/repo/kMoBzBp/images/1029093052-allineaDDT.PNG
   :alt: Allinea DDT

To use OpenMP, simply check the "OpenMP" box and specify how many OpenMP
threads you want to use


.. image:: https://bitbucket.org/repo/kMoBzBp/images/576706594-ddtopenmp.PNG
   :target: https://bitbucket.org/repo/kMoBzBp/images/576706594-ddtopenmp.PNG
   :alt: ddtopenmp


To use MPI, check the "MPI" box and specify how many MPI processes you
want to use. If this is gamera built with the recommended configuration,
it should be set to use "Intel MPI(MPMD)", and there should be no
additional arguments to mpiexec.hydra


.. image:: https://bitbucket.org/repo/kMoBzBp/images/1369899180-allineampi.PNG
   :target: https://bitbucket.org/repo/kMoBzBp/images/1369899180-allineampi.PNG
   :alt: allineampi


Once it is configured, click "Run" to start debugging. For more
information about using DDT, such as setting breakpoints and checking
the values of variables and arrays, there are DDT guides available
online `DDT Guide <https://developer.arm.com/docs/101136/latest/ddt>`_.

Attaching DDT to batch jobs
---------------------------

For particularly complicated or processor intensive cases, it is
possible to submit a job normally through the qsub submission system,
and have DDT attach to it once the job has started to run. This is
usually as simple as changing the application command inside the
submission script to have "ddt --connect" before it.

For example, if you are submitting a job which has the command:

.. code-block:: bash

   mpirun gamera.x


Change it to be the command:

.. code-block:: bash

   ddt --connect mpirun gamera.x

Before you submit this job (or at least before the job begins to run
after being in the queue), open DDT and let it sit at the main screen.
Once the job begins, a connection request dialog will appear in the DDT
GUI, and accepting it will begin the debugging session.

Running DDT With a Remote Client
--------------------------------

It is strongly recommended that you run DDT with the GUI local to you,
no matter where the software being debugged happens to be. This is
called using a Remote Client. The performance of the DDT GUI over x11 is
extremely poor, making the entire process slow and frustrating, and
rendering some options (such as plotting data) impossible to use.

Setting up DDT to use a Remote Client is straightforward, and an example
will be given here to run the GUI locally while debugging a job on
cheyenne.

First, you must locally download the EXACT same version of the Remote
Client DDT software as will be used remotely for debugging. So you must
download either the 20.2 or the 19.1 version of the Remote Client. In
the examples above, we are using arm-forge 19.1. The downloads page
`here <https://developer.arm.com/tools-and-software/server-and-hpc/downloads/arm-forge/older-versions-of-remote-client-for-arm-forge>`_
has all previous versions of the Remote Client for multiple operating
systems. You can download either version for your machine from that
page.

Once you have downloaded the Remote Client, it requires configuration to
be used with Cheyenne. The official documentation for this process can
be found
`here <https://developer.arm.com/documentation/101136/2010/Arm-Forge/Connecting-to-a-remote-system>`_.

In order to connect specifically to cheyenne, click on the "Remote
Launch" dropdown box, and then on "Configure...". Click on the "Add"
button to create a new configuration, and then set it up to look like
this, inserting your own username in the second box.


.. image:: https://bitbucket.org/repo/kMoBzBp/images/20146052-CheyenneRemote.PNG
   :target: https://bitbucket.org/repo/kMoBzBp/images/20146052-CheyenneRemote.PNG
   :alt: cheyenneremote


Once you've completed this, you can select "Cheyenne" at any time from
the "Remote Launch" dropdown on the DDT Remote Client main page, and it
will connect to Cheyenne, asking you to authenticate. Once it has
connected you launch your software using the "ddt --connect" option as
described above, and the connection request dialog will automatically
appear in your local DDT Remote Client, where you can perform the entire
debugging process.

The rest of the debugging process is the exact same as if the client
were not on your local machine. You can perform all functions, including
adding breakpoints and watchpoints, examining variables, restarting
sessions, etc....
