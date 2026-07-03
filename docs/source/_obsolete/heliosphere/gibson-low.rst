Running the Gibson & Low Model
==============================

.. Compiling Standalone
.. --------------------

.. Load modules for serial runs.
.. In your build directory

.. .. code-block:: bash

..    make gl.x

.. Model Configuration
.. -------------------

.. Specify model parameters in ``<prob ...>`` XML config, defaults parameters for spheromak are set to Schmitt et. al 2011 in Table 2. 

.. .. code-block:: xml

..    <?xml version="1.0"?>
..    <Kaiju>
..    <CME>
..        <sim runid="sphere64" isLoud="F" isDebug="F" isTopomorph="F" isStandalone="T" isAtmosphere="T" model="GL"/>
..        <time Tstart_transient="0.0" />
..        <output doTimer="T" />
..        <prob StateCoord="Sphere" orientation="1.57079" cmer="-1.57079" alpha="0.0" tfin="72000.0" dt="7200.0" frontheight="1.35" bmax="1.0" legsang="45.
..    0" vel_fh="50"/>
..        <idir min="1.0" max="21.0" N="64"/>
..        <jdir min="0.1" max="0.9" N="64"/>
..        <kdir min="0.0" max="1.0" N="64"/>
..    </CME>
..    </Kaiju>

.. The last three lines in the CME block define a 3D grid. Here, we set a spherical shell 1-21 R_Sun in radius, 0.1\ *pi-0.9*\ pi in polar direction(theta angle), 0-316 degrees in azimuthal direction (phi angle).

.. Input parameters have the following meanings. ``orientation`` means the orientation of a spheromak, e.g. ``0`` means the spheromak axis of asymmetry is perpendicular to the solar equatorial plane, ``1.57`` means that the axis of asymmetry is perpendicular to the meridional slice. ``frontheight`` is the distance from the center of the Sun to the front of a spheromak. ``bmax`` is the maximum value of the magnetic field strength in a spheromak set in [Gs], the default value is 0.001. ``legsang`` is the angular width of a CME in degrees. ``vel_fh`` is the velocity at the CME front in [km/s]. ``tfin`` is the total simulation time in [s]. ``dt`` is the time interval in [s] between data output.

.. See more examples of XML files here

.. .. code-block::

..    ~/kaiju/examples/gl

.. Postprocessing
.. --------------

.. You may then use the ``genXDMF.py`` script in order to process the standalone h5 output to visualize in Paraview. The "grid" is formatted such that it has coordinates X,Y,Z on a real sphere in units of rsun.
