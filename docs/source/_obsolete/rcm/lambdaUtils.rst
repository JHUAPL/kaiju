Generate and Test Lambda Distributions
======================================

Starting instructions for using kaipy/rcm/lambdautils to generate and test
lambda distributions for use in rcmconfigs.

Generating an rcmconfig.h5 file using defaults:

.. code-block:: python

   import kaipy.rcm.lambdautils.genAlam as genAlam
   from kaipy.rcm.lambdautils.AlamData import AlamParams

   params = AlamParams()   # All params can be set in constructor
                           # Otherwise use defaults in AlamParams init
   genAlam.genh5("rcmconfig.h5", params)       # Writes data to file

AlamParams is used to store the parameters used to generate electron and
proton lambda channels. The parameters, along with their defaults (as of
06-16-2021) are:

distType   : 'wolf' - Distribution type. Accepts 'lin', 'log', and 'wolf'

num_e      : 50     - Number of electron channels

num_p      : 149    - Number of proton channels

aMin_e     : -1.0   - Lower electron energy bound [eV]

aMin_p     : -1.0   - Lower proton energy bound [eV]

ktMax      : 50000  - Highest energy in [eV] that should be resolved at
GAM-RCM coupling boundary

L_kt       : 10     - L shell in [R_E] where ktMax should be resolved

tiote      : 4.0    - T_i/T_e ratio

p1         : 3.0    - p1 value used in 'wolf' lambda distribution generator

p2         : 1.0    - p2 value used in 'wolf' lambda distribution generator

addPsphere : True - Add 0-energy plasmasphere channel

In ``genAlam.genh5(params)``\ , ``doShowPlot=True`` can be set to show the
resulting lambda value vs. k channel plot. ``doTests=<bool>`` may also be set.
This will run the generated lambda channels through a series checks found in
``alamTester.py``. The output looks like:

.. code-block:: bash

   Smear test at L = 2.5: Passed
     Worst smear/cellWidth: 8.04e-01
   kT min/max range within 2.0% tolerance:
     Maxwellian : 4.33e-02/8.11e+01 [keV]
       Variance: D = 5.37e-09 P = 4.24e-05
     kappa =   5: 3.05e-01/4.04e+01 [keV]
       Variance: D = 4.45e-05 P = 2.01e-05

The smear test checks if any lambda channel is so wide that lambda- and
lambda+ would drift farther apart than a single grid cell over one GAM-RCM
coupling tilmestep. Ideally, ``smear/cellWidth < 1``.

The second check will try to find the range of temperatures between which the
lambda channels accurately reproduce the input density and pressure within a
given tolerance. If everything is -1, then the first checked temperature
(1 keV) was out of the tolerance range.

Many of the scripts in kaipy/rcm/lambdautils can be executed directly if
desired.
