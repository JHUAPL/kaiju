
SuperMAGE: SuperMAG indices
===========================

The Python module ``supermage.py`` (SuperMAGE) is part of the ``kaipy``
package. SuperMAGE is a collection of Python functions to compare simulated
ground magnetic field data generated from a MAGE magnetosphere run (using
``calcdb.x``), with `SuperMAG <https://supermag.jhuapl.edu/>`_ indices
(auroral SME/U/L and SMR). SuperMAGE provides functions to create index plots
and contour plots. Also included is a crude 1D E-Field calculator.

Requirements
------------

To use the SuperMAGE module, just import it into your Python code:

.. code-block:: python

   import supermage as sm

Reading simulated ground magnetic field perturbations from a MAGE run
---------------------------------------------------------------------

This step assumes you have run ``calcdb.x`` to compute simulated ground
magnetic field perturbations from your MAGE magnetosphere simulation results.
Instructions are provided `here <./groundMag>`_.

Assume you have the file ``geospace.deltab.h5``, created using ``calcdb.x``.
The data required for the comparison with SuperMAG measurements can be read in
as follows:

.. code-block:: python

   simulated_data = sm.ReadSimData("geospace.deltab.h5")

``simulated_data`` is a Python dictionary which contains the simulated data as
NumPy arrays. The following dictionary keys are available:

``td`` (``datetime.datetime``): Timestamps for the individual simulated data
points.

``glon`` (degrees): Geographic longitude

``glat`` (degrees): Geographic latitude

``mlat`` (degrees): Magnetic latitude

``mlt`` (hours): Magnetic local time

``smlon`` (degrees): Solar magnetic longitude

``dBt`` (nT): Spherical (polar angle) magnetic field perturbation at Earth
surface

``dBp`` (nT): Spherical (azimuthal angle) magnetic field perturbation at Earth
surface

``dBr`` (nT): Spherical (radial) magnetic field perturbation at Earth surface

``dBn`` (nT): Northward magnetic field perturbation at Earth surface

Fetching SuperMag magnetic indices
----------------------------------

We can use ``FetchSMIndices()`` to retrieve the indices corresponding to the
time period of our simulation:

.. code-block:: python

   SMI  = sm.FetchSMIndices(user, start, numofdays)

where


``user`` is the user name you created on the
`SuperMAG <https://supermag.jhuapl.edu/>`_ web site.

``start`` is a Python ``datetime.datetime`` object representing the starting
point of the desired data.

``numofdays`` is the number of days of data to fetch, starting at ``start``.

This call returns a Python dictionary of all the SuperMAG indices for the
chosen time period, stored as NumPy arrays. This call should only take a few
seconds.

As an example, assume we have read in the simulated data from the file
``storm1.deltab.h5`` described above. We can request 3 days of all SuperMAG
magnetic indices for this time period with the call:

.. code-block:: python

   SMI  = sm.FetchSMIndices(user, SIM['td'][0], 3)

Fetching SuperMAG magnetic observatory data
-------------------------------------------

The function ``FetchSMData()`` retrieves all of the available SuperMAG data
(in the form of observatory locations, magnetic field components, and
measurement times) for the desired time period, as shown below:

.. code-block:: python

   import os
   savefolder = os.path.join(os.environ['HOME'], 'supermag')
   SM = sm.FetchSMData(user, start, numofdays, savefolder, badfrac=0.1)

``user`` is the user name you created on the
`SuperMAG <https://supermag.jhuapl.edu/>`_ web site.

``start`` is a Python ``datetime.datetime`` object representing the starting
point of the desired data.

``numofdays`` is the number of days of data to fetch, starting at ``start``.

``savefolder`` is a string containing the path to a local cache folder for
SuperMAG data. If requested data is already available in this folder, it will
not be re-downloaded.

``badfrac`` (optional) is a threshold value (0 <= ``badfrac`` <= 1, default
``0.1``) for detecting bad data. If more than this fraction of the returned
data is missing or bad for an observatory, ignore that observatory.

This call returns a Python dictionary for all of the SuperMag observatory
data. The dictionary contains keys which list the observatory identifiers and
locations (in geographic and magnetic coordinates), measurement times, and the
observed magnetic field components.

This call can take several minutes per day of data to complete.

**Note**: The magnetic North component for this result is found using the
dictionary key ``BNm``.

Calculating magnetic indices for MAGE simulation results
--------------------------------------------------------

We can calculate a set of magnetic indices (SMR, SMU, SML, SMR) from our
simulation data using ``InterpolateSimData()``. This function performs several
operations:

#. Reject any simulated data which does not overlap with SuperMAG data.

#. Interpolate (in time) the simulated B data to the datetimes for the
SuperMAG indices already retrieved.

#. Interpolate (in space) the simulated B-data to the positions of the
SuperMag observatories

#. Calculates the indices listed above.

.. code-block:: python

   SMinterp = sm.InterpolateSimData(SIM, SM)

This call should complete in a few seconds.

The Python dictionary returned by this call contains the times from the
simulated results, the locations of each SuperMAG observatory

The dictionary key ``dBn`` will return the simulation data for the SuperMag
locations, whereas ``dBnsmall`` returns the simulation data interpolated to.

Making Plots
^^^^^^^^^^^^

There are two functions to make comparison plots: **MakeIndicesPlot** and
**MakeContourPlot**. These can be called as follows:

.. code-block:: python

   sm.MakeIndicesPlot(SMI, SMinterp, fignumber = 1)
   sm.MakeContourPlots(SM, SMinterp, maxx = 1000, fignumber = 2)

These return plots like the following respectively (\ **MakeContourPlots()**
make take while):

.. image:: https://bitbucket.org/repo/kMoBzBp/images/2400869799-ComparisonPlots.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/2400869799-ComparisonPlots.png
   :alt: ComparisonPlots.png

"Interpolated" is simulated data interpolated at SuperMag locations, "super"
is simulated data using all simulated data and "real" is actual SuperMag data.

Calculating E-Field
^^^^^^^^^^^^^^^^^^^

Included is the **EField1DCalculation()** function which calculates horizontal
magnetic field components using a resistive 1-D ground resistivity model (the
Quebec model found in DOI:10.1046/j.1365-246x.1998.00388.x).

Note: this function cannot handle nan values in the time-series. Any nan in
the time-series will return all nans for that individual location. Also, this
function is quite slow...

.. code-block:: python

   BX, BY, TD = SM['BNm'], SM['BEm'], SM['td']
   EX, EY = sm.EField1DCalculation(BX, BY, TD)
   SM['Ex'] = EX
   SM['Ey'] = EY

Example Run
-----------

.. code-block:: python

   # NEED IMPORT STATEMENTS HERE.

   user = 'USERNAME' # username for Supermag downloads
   savefolder = '/glade/u/home/sblake/indices/FINAL/DATA/' # folder where the Supermag jsons are stored/fetched

   # Working Example
   # Read in all needed data from .h5 simulation output
   filename = '/glade/u/home/skareem/Work/gemmie/stpatex/stpatdb.deltab.h5'
   #filename = '/glade/u/home/skareem/Work/dass/Data/newdassh/dassdb.deltab.h5'
   print("Reading SIM data")
   SIM = ReadSimData(filename)
   start = SIM['td'][0] # start time of simulation data
   numofdays = 3        # going to

   print("Fetching SM indices")
   SMI  = FetchSMIndices(user, start, numofdays)

   print("Fetching SM data")
   SM = FetchSMData(user, start, numofdays, savefolder, badfrac = 0.1)

   print("Interpolating SIM data") # interpolates and calculates SM indices
   SMinterp = InterpolateSimData(SIM, SM)

   print("Making Indices Plot")
   MakeIndicesPlot(SMI, SMinterp, fignumber = 1)

   print("Making Contour Plot")
   MakeContourPlots(SM, SMinterp, maxx = 1000, fignumber = 2)

   print("Calculating E-Field for SM data")
   BX, BY, TD = SM['BNm'], SM['BEm'], SM['td']
   EX, EY = EField1DCalculation(BX, BY, TD)
   SM['Ex'] = EX
   SM['Ey'] = EY
