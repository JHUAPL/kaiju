Running Gamera Helio with an embedded CME
=========================================

.. Compile ``gamhelio.x`` or ``gamhelio_mpi.x`` as described in 
.. :doc:`Helio Quickstart </quickStart/helioQuickStart>` using the ``giblow``
.. branch. 

.. You will then need to add the ``<CME />`` element to the input xml fle and fill out the appropriate CME model parameters. Currently, only the Gibson-Low flux rope model is implemented and you may supply related model parameters in the ``<prob />`` element.  

.. Some reasonable defaults are below representing a slightly elongated un-tethered spheromak launched centered on the equatorial plane at longitude = 2.95 [rad].

.. .. code-block:: xml

..    <?xml version="1.0"?>
..    <Kaiju>
..    <Gamera>
..        <sim runid="wsa256cmeLon295" doH5g="T" H5Grid="heliogrid.h5" icType="user" pdmb="1.0" rmeth="7UP"/>
..        <time tFin="72.0"/>
..        <spinup doSpin="T" tSpin="160.0" tIO="0.0"/>
..        <output dtOut="1.0" tsOut="100" doTimer="T"/>
..        <physics doMHD="T" gamma="1.5"/>
..        <prob doCME="T" rotateCME="T" isSpSymSW="F" DeltaT="0.0" model="monopole"/>
..        <helio Tsolar ="25.38" vrin="385.0" vrkin="385.0" tin="4.d5" rhoin="800.0" brin="0.000" brkfin="0.000" />
..        <restart resId = "wsa" nRes = "00008" dtRes="6." doRes="F"/>
..        <!--- Next 3 lines only necessary for MPI runs, ignored for serial runs ---> 
..        <iPdir N="4" bcPeriodic="F"/>
..        <jPdir N="4" bcPeriodic="F"/>
..        <kPdir N="4" bcPeriodic="T"/>
..    </Gamera>
..    <CME>
..        <sim isLoud="T" isDebug="F" isTopomorph="F" isStandalone="F" model="GL" scaleBmax="T"/>
..        <time Tstart_transient="0.0" />
..        <output doTimer="T" />
..        <prob Den_CME = "800.0" T_CME = "1000000.0" orientation="1.57079" lat="0.0" lon="2.95" alpha="0.0" frontheight="21.5" legsang="60.0" apar="0.05" Bmax="0.003" vel_fh="900.0" />
..    </CME>
..    </Kaiju>
