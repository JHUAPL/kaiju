VOLTRON Block - Coupler
=======================

VOLTRON
-------

Example XML:

.. code-block:: xml

   <VOLTRON>
     <time tFin="3601.0"/>
     <output dtOut="60" tsOut="100" doTimer="F"/>
     <spinup doSpin="T" tSpin="1800.0" tIO="0.0"/>
     <restart dtRes="1800.0"/>
     <imag doInit="T"/>
     <coupling dt="5.0" doGCM="T" dtDeep="15.0" rTrc="40.0" imType="RCM" doQkSquish="T"/>
   </VOLTRON>

VOLTRON/time
------------

tFin: [seconds] - What time (in seconds) from start of simulation to stop run
at

VOLTRON/output
--------------

* dtOut: [seconds]
* tsOut: [?]
* doTimer: [True/False] - ?

VOLTRON/restart
---------------

dtRes: [seconds] - What cadence we will output a restart file? 1800.0 =
every 30 minutes

VOLTRON/imag
------------

* doInit: [True/False] - Do we do Kareem's initialized ring current?

VOLTRON/coupling - The most important section
---------------------------------------------

* dt: [seconds] - How often do we couple with REMIX (and GCM)
* doGCM: [True/False] - Are we coupling to a GCM?
* dtDeep: [seconds] - How often do we couple to RCM?
* rTrc: [?] - ?
* imType: [RCM] - ?
* doQkSquish: [True/False] - this has to do with the very expensive operation
  of taking every (basically) cell on the gamera grid and projecting it to the
  northern hemisphere
* doSerial: [True/False] - Voltron runs concurrently with Gamera by default.
  This can be set to True to force Voltron and Gamera to run serially, taking
  turns and waiting for each other.

Gamera Block - The MHD stuff
----------------------------

Example XML:

.. code-block::

   <Gamera>
     <sim runid="msphere" doH5g="T" H5Grid="lfmQ.h5" icType="user" pdmb="1.0" pFloor="1.0e-8" dFloor="1.0e-6" rmeth="7UP"/>
     <restart doRes="F" nRes="-1" resID="msphere"/>
     <physics doMHD="T" doBoris="T" Ca="10.0"/>
     <ring gid="lfm" doRing="T"/>
     <wind tsfile="bcwind.h5"/>
     <iPdir N="4" bcPeriodic="F"/>
     <jPdir N="4" bcPeriodic="F"/>
     <kPdir N="1" bcPeriodic="T"/>
     <source doSource="T"/>
   </Gamera>

Gamera/sim
----------

REMIX Block - The thing that solves for Potential
-------------------------------------------------

Example XML:

.. code-block::

   <!-- Remix params -->
   <REMIX>
     <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
     <conductance F107="100.0" pedmin="2.0" hallmin="1.0" sigma_ratio="3.0" const_sigma="False" ped0="10.0"/>
     <precipitation aurora_model_type="RCMONO" alpha="0.2" beta="0.4"/>
   </REMIX>

REMIX/grid
----------

* Np: [integer] - How many bins in longitude direction (360/Np where Np=360
  means 1°).
* Nt: [integer] - how many bins in colatitude direction ( (90-LowLat)/Nt where
  Nt=45 and LowLat=45 means 1°).
* LowLatBoundary: [degree] - What colatitude degree is your Low Latitude
  Boundary?

REMIX/conductance
-----------------

There are a lot of conditionals within this section. So this section may
become complicated.

* const_sigma: [True/False] - Do we want to use uniform constant conductance?

If const_sigma = True

* ped0: [float] - set pedersen conductance to a uniform ped0 value

If const_sigma = False

F107: [float] - What SFU is the F10.7 for the run?

pedmin: [float] - What is the minimum pedersen conductance

hallmin: [float] - what is the minimum hall conductance

sigma_ratio: [float] - What is the maximum Ped/Hall that we allow?

ped0: [float] - Even if const_sigma = False, this is used to specify the
conductance during spin-up (-60min) period.

If doGCM = True

pedmin: [float] - minimum pedersen conductance [min(gcmped,pedmin)]

hallmin: [float] - minimum hall conductnace [min(gcmhall,hallmin)]

REMIX/precipitation
-------------------

auroral_model_type: [FEDDER,RCMONO] - Do we use Fedder precipitation or
RCM+Zhang Mono precipitation? Keep in mind RCMONO requires RCM or it will do
weird things.

alpha: [float] - Alpha parameter that specifies the Ti/Te temperature ratio.
[default RCMONO=0.2, FEDDER=1.0332467]

beta: [float] - Beta parameter that approximately translates to loss cone
filling factor. [default RCMONO = 0.4, FEDDER = 0.4362323]

R: [float] - Only used by FEDDER. [default FEDDER=0.083567956]

CHIMP Block - Crazy fast field line tracer
------------------------------------------

Example XML:

.. code-block:: xml

   <!-- EB-Tracer (CHIMP) params -->
   <CHIMP>
     <units uid="EARTHCODE"/>
     <fields grType="LFM"/>
     <domain dtype="SPH" rmin="2.0" rmax="40.0"/>
     <tracer epsds="0.05"/>
   </CHIMP>

RCM Block - The Rice Convection Model
-------------------------------------

Example XML:

.. code-block:: xml

   <RCM>
     <output debug="F" toRCM="F" toMHD="F"/>
     <clawpack doKaiClaw="F"/>
     <ellipse xSun="15.0" xTail="-15.0" yDD="15.0" isDynamic="T" dRadMHD="1.0"/>
   </RCM>

RCM/ellipse
-----------

xSun: [float] - maximum positive X value of RCM grid in RE

xTail: [float] - furthest down tail RCM grid can go.

yDD: [float] - how wide can RCM's flanks be in the Y direction

isDynamic: [True/False] - Turn on dynamic plasmasphere vs gallagher.

dRadMHD: [float] - Magic number

Full Example XML - By your powers combined, I am VOLTRON/MiniMAGE!
------------------------------------------------------------------

Example XML:

.. code-block:: xml

   <Kaiju>
   <!-- Example XML file for coupled RCM+Gamera+ReMIX+CHIMP -->
   <?xml version="1.0"?>
   <!-- Magnetosphere params, Voltron times in seconds -->
   <VOLTRON>
     <time tFin="3601.0"/>
     <output dtOut="60" tsOut="100" doTimer="F"/>
     <restart dtRes="1800.0"/>
     <imag doInit="T"/>
     <coupling dt="5.0" doGCM="F" dtDeep="15.0" rDeep="8.0" rTrc="16.0" imType="RCM"/>
   </VOLTRON>
   <Gamera>
     <sim runid="msphere" doH5g="T" H5Grid="lfmQ.h5" icType="user" pdmb="1.0" pFloor="1.0e-8" dFloor="1.0e-6" rmeth="7UP"/>
     <restart doRes="T" nRes="0" resID="msphere"/>
     <physics doMHD="T" doBoris="T" Ca="10.0"/>
     <ring gid="lfm" doRing="T"/>
     <wind tsfile="bcwind.h5"/>
     <iPdir N="4" bcPeriodic="F"/>
     <jPdir N="4" bcPeriodic="F"/>
     <kPdir N="1" bcPeriodic="T"/>
     <source doSource="T"/>
   </Gamera>
   <!-- Remix params -->
   <REMIX>
     <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
     <conductance F107="100.0" pedmin="2.0" hallmin="1.0" sigma_ratio="3.0" const_sigma="False" ped0="10.0"/>
     <precipitation aurora_model_type="RCMONO" alpha="0.2" beta="0.4"/>
   </REMIX>
   <!-- EB-Tracer (CHIMP) params -->
   <CHIMP>
     <units uid="EARTHCODE"/>
     <fields grType="LFM"/>
     <domain dtype="SPH" rmin="2.0" rmax="25.0"/>
     <tracer epsds="0.05"/>
   </CHIMP>
   <RCM>
     <output debug="F" toRCM="F" toMHD="F"/>
     <clawpack doKaiClaw="F"/>
     <ellipse xSun="15.0" xTail="-15.0" yDD="105.0" isDynamic="T" dRadMHD="1.0"/>
   </RCM>
   </Kaiju>
