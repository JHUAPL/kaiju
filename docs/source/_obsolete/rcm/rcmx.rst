Running rcm.x
=============

These are instructions to run the RCM stand-alone driver.

Follow the instructions in the quick start section, once the standard drivers
can be compiled you can use the same build system to run "make rcm.x" which
will create the driver.

In addition to the binary "rcm.x" there are two additional files that are
needed. The first is the RCM config which can be generated using the
"genRCM.py" script under "kaiju/scripts", this will create a file
"rcmconfig.h5". The second is an input deck XML, an example of which is below

.. code-block:: xml

   <RCM>
     <output debug="F" toRCM="F" toMHD="F"/>
   </RCM>

Assuming the XML deck is called "rcmtest.xml" then the executable can be run
simply with "rcm.x rcmtest.xml". The experimental new clawpack version can be
enabled by editing "kaiju/src/rcm/rcm_include.h" and changing "doClaw95" to
true.

.. code-block:: fortran

       INTEGER (iprec), PARAMETER :: &
          isize = 200, &
          jsize = 101, &
          jwrap =   3, &
          ksize  = 090,  kcsize = 090, &
          nptmax = 50000, &
          iesize =   2, &
          ncoeff =   5
       LOGICAL, PARAMETER :: asci_flag = .FALSE.
       LOGICAL, PARAMETER :: isGAMRCM = .TRUE. !Whether running coupled to Gamera
       LOGICAL, PARAMETER :: doQuietRCM = .TRUE.
       LOGICAL, PARAMETER :: doClaw95 = .FALSE.
