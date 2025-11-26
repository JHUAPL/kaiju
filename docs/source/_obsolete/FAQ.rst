Frequently Asked Questions
==========================

This page contains frequently asked questions about running the ``kaiju``
software.

Q1. I got Intel compiler warning #11021 (unresolved references in dynamic
libraries) when linking the Fortran executable ``gamera.x``. For example:

.. code-block:: bash

    $ make gamera.x
    ...
    'ipo: warning #11021: unresolved files_mp\ *checkandkill*
    Referenced in libgamlib.a(output.F90.o)
    Referenced in libgamlib.a(gioH5.F90.o)
    Referenced in libgamlib.a(init.F90.o)'

A1. Invalid directories may be included in ``$PATH`` and/or ``$PYTHONPATH``.

Q2. What's the code unit for current density from MHD output?

A2. ``gB0/gx0/mu0 = 4.58nT/6.38e6m/(4pi*1e-7)=5.7139e-10 A/m^2``

Check the units of other variables in the HDF5 file with ``h5info`` in Matlab
or similar functions in ``python``.
