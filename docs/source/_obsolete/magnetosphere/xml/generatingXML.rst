Generating XML Files
====================

Kaiju requires XML files to provide key information on how a run is to occur.
This wiki page should help ensure that users can generate their own XML files,
convert old XML files into new ones, and use the new XMLGenerator script to
use configuration (.ini) files to make the XML generating process easier.

Converting Old XML To The New Format
------------------------------------

A good portion of users likely already have XML files that they have been
using for a while. Unfortunately, those are likely broken now. *Technically*,
they were never XML files to begin with because *technically* XML must have a
single root node.

As a result, all XML files must now contain the ``<Kaiju></Kaiju>`` root node.
To convert an XML that previously worked:

Add ``<Kaiju>`` to the top of the file alone on the first line.

Add "</Kaiju>" to the bottom of the file all alone on the last line.

Indent everything in between by one level.

So an old XML file that will no longer work might look like this:

.. code-block:: xml

   <voltron>
       <sim tFin="677.115987461"/>
       <spinup doSpin="F"/>
       <coupling dt="5.0"/>
   </voltron>
   <Gamera>
       <iPdir N="1"/>
       <jPdir N="2"/>
       <kPdir N="1"/>
       <ibc bc="user" jperiod="F" kperiod="T"/>
       <jbc bc="user" iperiod="F" kperiod="T"/>
       <kbc bc="periodic" iperiod="F" jperiod="F"/>
       <sim runid="msphere" doH5g="T" H5Grid="lfmD.h5" icType="user" pdmb="4.0" rmeth="7UP" pFloor="1e-06" dFloor="0.0001"/>
       <restart dtRes="14.1065830721" doRes="F" resFile="msphere.Res.XXXXX.h5"/>
       <output dtOut="0.235109717868" tsOut="100" timer="F"/>
       <physics doMHD="T" doBoris="T" Ca="10.0"/>
       <prob Rho0="0.1" P0="0.001" rCut="16.0" lCut="8.0"/>
       <ring gid="lfm" doRing="T" Nr="3" Nc1="8" Nc2="16" Nc3="32"/>
       <wind tsfile="bcwind.h5"/>
   </Gamera>
   <REMIX>
       <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
       <conductance const_sigma="T" ped0="10.0"/>
   </REMIX>

Whereas the corrected XML file will look like this:

.. code-block:: xml

   <Kaiju>
       <voltron>
           <sim tFin="677.115987461"/>
           <spinup doSpin="F"/>
           <coupling dt="5.0"/>
       </voltron>
       <Gamera>
           <iPdir N="1"/>
           <jPdir N="2"/>
           <kPdir N="1"/>
           <ibc bc="user" jperiod="F" kperiod="T"/>
           <jbc bc="user" iperiod="F" kperiod="T"/>
           <kbc bc="periodic" iperiod="F" jperiod="F"/>
           <sim runid="msphere" doH5g="T" H5Grid="lfmD.h5" icType="user" pdmb="4.0" rmeth="7UP" pFloor="1e-06" dFloor="0.0001"/>
           <restart dtRes="14.1065830721" doRes="F" resFile="msphere.Res.XXXXX.h5"/>
           <output dtOut="0.235109717868" tsOut="100" timer="F"/>
           <physics doMHD="T" doBoris="T" Ca="10.0"/>
           <prob Rho0="0.1" P0="0.001" rCut="16.0" lCut="8.0"/>
           <ring gid="lfm" doRing="T" Nr="3" Nc1="8" Nc2="16" Nc3="32"/>
           <wind tsfile="bcwind.h5"/>
       </Gamera>
       <REMIX>
           <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
           <conductance const_sigma="T" ped0="10.0"/>
       </REMIX>
   </Kaiju>

Generating New XML Files From Scratch
-------------------------------------

Config Files
^^^^^^^^^^^^

The XMLGenerator script uses configuration (.ini) files to generate XML files.
A config file looks like this:

.. code-block:: python

   ## MODULE NAME ##
   [Section Name 1]
   optionName1 = T

   [Section Name 2]
   optionName2 = 20.0
   optionName3 = hello

   ## MODULE NAME 2 ##
   [Section Name 3]
   optionName4 = junk # Comments

"MODULE NAME" Corresponds to things like VOLTRON, Gamera, and CHIMP.

"[Section Name]" Corresponds to certain sections underneath modules. For
example, time, output, and spinup under VOLTRON.

"optionName" Corresponds to an option underneath section headers. For
example, tFin, dtOut, and doTimer. The values that can be assigned to these
can be strings, integers, floats, or even comments.

Anything that follows a "#" is treated as a comment in the vanilla
configuration file parser. However, we have added some more functionality that
will be covered later.

Let's take a look at an XML that can be used with Kaiju and then see what it
looks like as a config file. Here is is what "kaiju/tests/remix/cmiD.xml"
looks like as an XML:

.. code-block:: xml

   <Kaiju>
       <voltron>
           <sim tFin="677.115987461"/>
           <spinup doSpin="F"/>
           <coupling dt="5.0"/>
       </voltron>
       <Gamera>
           <iPdir N="1"/>
           <jPdir N="2"/>
           <kPdir N="1"/>
           <ibc bc="user" jperiod="F" kperiod="T"/>
           <jbc bc="user" iperiod="F" kperiod="T"/>
           <kbc bc="periodic" iperiod="F" jperiod="F"/>
           <sim runid="msphere" doH5g="T" H5Grid="lfmD.h5" icType="user" pdmb="4.0" rmeth="7UP" pFloor="1e-06" dFloor="0.0001"/>
           <restart dtRes="14.1065830721" doRes="F" resFile="msphere.Res.XXXXX.h5"/>
           <output dtOut="0.235109717868" tsOut="100" timer="F"/>
           <physics doMHD="T" doBoris="T" Ca="10.0"/>
           <prob Rho0="0.1" P0="0.001" rCut="16.0" lCut="8.0"/>
           <ring gid="lfm" doRing="T" Nr="3" Nc1="8" Nc2="16" Nc3="32"/>
           <wind tsfile="bcwind.h5"/>
       </Gamera>
       <REMIX>
           <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
           <conductance const_sigma="T" ped0="10.0"/>
       </REMIX>
   </Kaiju>

Here is the exact same file, but written as a config file:

.. code-block:: ini

   ## voltron ##
   [sim]
   tFin = 677.115987461

   [spinup]
   doSpin = F

   [coupling]
   dt= 5.0

   ## Gamera ##
   [iPdir]
   N = 1

   [jPdir]
   N = 2

   [kPdir]
   N = 1

   [ibc]
   bc = user
   jperiod = F
   kperiod = T

   [jbc]
   bc = user
   iperiod = F
   kperiod = T

   [kbc]
   bc = periodic
   iperiod = F
   jperiod = F

  [sim]
   runid = msphere
   doH5g = T
   H5Grid = lfmD.h5
   icType = user
   pdmb = 4.0
   rmeth = 7UP
   pFloor = 1e-06
   dFloor = 0.0001

   [restart]
   dtRes = 14.1065830721
   doRes = F
   resFile = msphere.Res.XXXXX.h5

   [output]
   dtOut = 0.235109717868
   tsOut = 100
   timer = F

   [physics]
   doMHD = T
   doBoris = T
   Ca = 10.0

   [prob]
   Rho0 = 0.1
   P0 = 0.001
   rCut = 16.0
   lCut = 8.0

   [ring]
   gid = lfm
   doRing = T
   Nr = 3
   Nc1 = 8
   Nc2 = 16
   Nc3 = 32

   [wind]
   tsfile = bcwind.h5

   ## REMIX ##
   [grid]
   Np = 360
   Nt = 45
   LowLatBoundary = 45.0

   [conductance]
   const_sigm = T
   ped0 = 10.0

Creating A Usable XML File With XMLGenerator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you already have a .ini file that is formatted properly, then the steps for
generating a usable XML file are fairly simple.

Ensure that you have your $KAIJUDIR/scripts folder on your $PATH, $KAIJUDIR on
your $PYTHONPATH, and your python environment setup according to the
[[Quick Start]] guide.

Enter the following command: "XMLGenerator.py myConfigFile.ini
myOutputXMLFile.xml"

If you get a message saying "Template creation complete!", then it succeeded
and the "myOutputXMLFile.xml" (or whatever you decided to call it) should
exist in your folder, ready for use.

Using Config Files To Modify Existing XML Files
-----------------------------------------------

Ensure that you have your $KAIJUDIR/scripts folder on your $PATH, $KAIJUDIR on
your $PYTHONPATH, and your python environment setup according to the
[[Quick Start]] guide.

\Ensure that you have a valid XML file to use as a template for the new one
you are about to generate.

Enter the following command:
"XMLGenerator.py myTemplateXMLFile.xml myConfigFile.ini myOutputXMLFile.xml"

If you get a message saying "XML generation complete!", then it succeeded and
the "myOutputXMLFile.xml" (or whatever you decided to call it) should exist in your folder, ready for use.

Other Functionality
-------------------

There are a couple of other things you can do with your .ini configuration
files.

Time Unit Conversion
^^^^^^^^^^^^^^^^^^^^

Some settings are in units of time, like tFin, and the XML requires that to
be in seconds. However, a run can be set to run for a while and that means
this number can get pretty big pretty quickly. As such, the XMLGenerator
script allows you to set units for time settings in seconds, minutes, or
hours. So instead of writing 43,200 seconds, you can just write 12 hours.

Here is an example of a time setting in an .ini file without any special
options:

.. code-block:: ini

   [sim]
   tFin = 43200

Here is the same line but with a special flag that says "This value is in
minutes":

.. code-block:: ini

   [sim]
   tFin = 720 # [min]

And finally, one written in hours:

.. code-block:: ini

   [sim]
   tFin = 12 # [hrs]

To add these flags to your .ini files follow these rules:

The flag must be placed on the same line as the value you want to convert

There must be a "#" placed after the value AND before the flag

The following options are valid: "[sec]", "[min]", "[hrs]".

If you use an incorrect unit like "dogs", the following message will appear:

.. code-block:: ini

   ERROR: Incorrect unit type for conversion: dogs

All of these flags will ensure that the final output of the XML file has your
time values in seconds; just what Kaiju needs!

Deleting Options
^^^^^^^^^^^^^^^^

Let's say there is an option in a template XML file that I don't want in my
final, output XML file. In that case, you can use the "!DEL" flag to tell the
script you don't want this value included.

There are two ways to use this flag: As the value of the option you want to
delete, or as a comment like the unit conversion flags.

Here is an example of the first method:

.. code-block:: ini

  [sim]
   tFin = !DEL

And here is an example of the second:

.. code-block:: ini

  [sim]
   tFin = 12 # !DEL

Both methods will work to make sure that tFin will not appear in your final
XML file.

Dependencies
------------

The XMLGenerator script requires the following dependencies:

`Config Parser <https://pypi.org/project/config-parser/>`_
