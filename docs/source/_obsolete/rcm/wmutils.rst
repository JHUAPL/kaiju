Diffuse Precipitation
=====================

To run MAGE with wave model input:

Make sure DWang_chorus_lifetime.h5 and tauTDS.txt is inside
kaiju/kaipy/rcm/wmutils

Activate the NLP environment:

.. code-block:: python

   ncar_pylib casper_satcomp_pylib

To generate rcmconfig.h5 that contains the electron lifetime based on the wave
model, run kaiju/scripts/preproc/genRCM.py with 'waveModel' option on.

.. code-block:: python

   python genRCM.py -waveModel True

To add the electron lifetime to an existing rcmconfig.h5 file, run
kaiju/scripts/preproc/genRCM.py with 'addWM' option on and enter input file
name.

.. code-block:: python

   python genRCM.py -addWM True -i rcmconfig.h5

The generated rcmconfig.h5 should contain


* Eki Dataset {155}
* Kpi Dataset {7}
* Li Dataset {41}
* MLTi Dataset {25}
* Tau1i Dataset {155, 41, 25, 7}
* Tau2i Dataset {155, 41, 25, 7}
* EkTDSi           Dataset {109}
* TauTDSi          Dataset {109}
* alamc Dataset {160}
* dktable Dataset {936}
* fudgec  Dataset {160}
* ikflavc Dataset {160}

In the XML file, to run MAGE with wave model input (Dedong Wang Chorus +
Orlova16 Hiss), set the loss method in the RCM section as

.. code-block:: xml

   <RCM>
        <loss eLossMethod="WM"/>
   </RCM>

Go to Step 1. if you see “Wave model is missing in rcmconfig.h5” in .out file
of the MAGE run.
