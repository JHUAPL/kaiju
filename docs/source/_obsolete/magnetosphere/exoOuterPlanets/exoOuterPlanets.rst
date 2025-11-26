Non-Earth Planetary Magnetospheres
==================================

Instructions for running GAMERA-REMIX for non-Earth planetary magnetosphers.

By default, the voltron executable is designed to simulate Earth's
magnetosphere, and the config file assumes that the planet is Earth. This page
outlines how to run the model for non-Earth planetary magnetospheres, and what
assumptions are made in the model when doing so.

Model choices you should know about
-----------------------------------

Coordinate systems
^^^^^^^^^^^^^^^^^^

GAMERA and REMIX assume that the magnetic dipole axis is along the +z
direction. This means **the magnetic moment should always be positive**. The
sign of the corotation potential should be with respect to the orientation of
the rotation axis with the dipole axis. i.e., the corotation potential should
be positive if the two axes are aligned (e.g. Earth) and negative if they are
anti-aligned (e.g. Jupiter and Saturn). There is currently no correction to
the solar wind parameters with regard to the planet's coordinate system, so
you must do that conversion ahead of time.

IC files
--------

The default IC file, (src/voltron/ICs/earthcmi.F90), is specific to Earth's
magnetosphere. An alternative IC file, src/voltron/ICs/planetcmi.F90, can be
used as a more generic starting point. This may be used as-is for some simple
applications where only the planet parameters are changed.

You can configure voltron to build with this IC file either by setting the
``VOLTIC`` variable in ccmake, or by setting it as an argument with cmake:

``cmake <other_args> -DVOLTIC=$KAIJUHOME/src/voltron/ICs/planetcmi.F90 ..``

More specific configurations (e.g. Jupiter & Saturn with special initial
conditions, continuous mas loading from moons, etc.) should be defined in new
IC files using planetcmi.F90 as a starting point.

Config/XML file options
-----------------------

The planet parameters must also be specified in the xml config file. Here is
an example snippet:

.. code-block:: xml

   <Kaiju>
     <Gamera>
       <prob planet="OTHER" x0="4913961" Rion="1.1" M0="0.1" Psi0="0"/>
     </Gamera>
   </Kaiju>

where

``Planet`` - Planet name. Mercury, Earth, Jupiter, and Saturn all have defined
parameters that don't need to be specified in the xml. "Other" can be used to
simulate a custom planet, where all parameters default to Earth unless
overwritten by any of the following:

``x0`` - Planet radius in meters

``Rion`` - Planetary ionosphere radius in planetary radii

``M0`` - Magnetic moment in Gauss (\ **Should always be positive.** See above
about coordinate systems)

``G0`` - Gravitational acceleration at planet surface in m/s^2 (e.g. 9.81 for
Earth)

``doGrav`` - uses ``G0`` to add gravitational term in GAMERA

``Psi0`` - Corotation potential in kV (See above about coordinate systems to
determine if it should be positive or negative)

``doCorot`` - uses ``Psi0`` to enforce corotation

It is also highly advised to tell REMIX to use a constant conductance, as the
more complex models are specific to Earth:

.. code-block:: xml

     <REMIX>
       <conductance const_sigma="T" ped0="1"/>
     </REMIX>
