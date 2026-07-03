Building the ``kaiju`` software
===============================


Introduction
------------

This section describes how to build the ``kaiju`` software on two different
supercomputers - ``derecho`` and ``aitken``. If you are trying to build the
``kaiju`` software on a different system, use these instructions as a starting
point.


Before you begin
----------------

The ``kaiju`` software is typically built with the Intel Fortran compiler
(although the software can also be built with other Fortran compilers,
e.g., GNU). Building the software also requires the ``cmake`` build tool, the
``HDF5`` library, and an ``MPI`` library. Instructions for loading these
packages are provided in the ``module`` commands for each HPC system.


Getting the source code
-----------------------

The ``kaiju`` source code can be obtained by cloning the ``kaiju`` repository
on GitHub:

.. code-block:: bash

    git clone https://github.com/JHUAPL/kaiju.git

.. important::

    The ``kaiju`` repository on GitHub uses ``git-lfs`` to support the
    use of large binary files. You *must* make sure ``git-lfs`` is
    available in your ``git`` installation to ensure a complete clone
    of the ``kaiju`` repository.

.. note:: The ``kaiju`` software can be built in serial or MPI versions. The
    serial version is best for single-processor machines such as a laptop,
    while supercomputers such as ``derecho`` and ``aitken`` typically use
    the MPI version, to take advantage of multiple compute nodes. These
    instructions describe how to build the MPI version of ``kaiju``. The build
    instructions for the single-machine serial version are essentially the
    same as for the MPI version - typically all that is required is a Fortran
    compiler and an HDF5 library.


Building Guides
-----------------

.. toctree::
    :maxdepth: 1

    buildDerecho_GR
    buildAitken_GR
    buildDerecho_GTR
    buildAitken_GTR
