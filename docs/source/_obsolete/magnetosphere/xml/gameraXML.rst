GAMERA XML
==========

Example XML file

.. code-block:: xml

   <OMEGA>
     <sim tFin="677.115987461"/>
   </OMEGA>
   <Gamera>
     <sim runid="msphere" doH5g="T" H5Grid="lfmD.h5" icType="user" pdmb="1.0" pFloor="1.0e-8" dFloor="1.0e-6" rmeth="7UP"/>
     <restart dtRes="28.2131661" doRes="F" resFile="XXXXX.h5" doReset="F" tReset="0.0"/>
     <output dtOut="0.940438871" tsOut="100" timer="F"/>
     <physics doMHD="T" doBoris="T" Ca="10.0"/>
     <prob Rho0="0.2" P0="0.001"/>
     <ring gid="lfm" doRing="T" Nr="4" Nc1="8" Nc2="16" Nc3="32" Nc4="32"/>
     <wind tsfile="bcwind.h5"/>
   </Gamera>

Note, XML keys are not case sensitive but input strings (like filenames) are.

Overall run time of simulation (in code units) specified in omega/sim/tFin.

gamera/sim

runid: String to prepend to output files from this run

doH5g/H5Grid: Whether or not to read grid from HDF5 file and if so the grid
filename

icType: Whether to do a canned problem file from prob.F90 or use the
user-specified IC routine (GAMIC cmake variable)

pdmb: Value of beta for partial donor method, 1.0 is standard (see paper)

pFloor/dFloor: Global floors (in code units) for pressure/density

rmeth: Spatial reconstruction method, standard choices are 8CENT/7UP (see
paper)

gamera/restart

dtRes: Time interval (in code units) to output restart dumps

doRes/resFile: Whether to load a restart dump for this run

doReset/tReset: If restarting whether to reset the code time to tReset

gamera/output


dtOut: Time interval (in code units) to output 3D data slices

tsOut: Interval (in timesteps) to write output to console

timer: Whether to output to console timing data

gamera/physics

doMHD: If this is an MHD run

gamma: Adiabatic index

doBoris/Ca: Whether to use Boris correction and if so the Boris speed of light
(in code units)

gamera/prob: Various parameters specified by the IC file

planet: Specify parameters for "EARTH", "SATURN", "JUPITER", "MERCURY",
"NEPTUNE", or "OTHER"

If "OTHER", you may specify values for:

x0: planet radius [m]

v0: velocity scale [m/s]

G0: gravity [m/s2]

M0: magnetic moment [gauss]

Psi0: corotation potential [kV]

Rion: ionosphere radius [x0]

doGrav: True or False

If any values are not specified, they will default to Earth values.

gamera/ring

doRing: Whether to use ring-average

Nr: Number of rings to use followed by Nc1,Nc2 ... specifying the chunking for
each ring
