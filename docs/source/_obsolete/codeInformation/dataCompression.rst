Data compression in HDF5 files
==============================

Datasets in each output time step (HDF5 Group) can be compressed using a
number of algorithms yielding varying levels of performance and precision
(lossless vs. lossy). The default algorithm is
`SZIP <https://support.hdfgroup.org/documentation/hdf5/latest/group___s_z_i_p.html>`_.
Other algorithms are available: `zlib <https://www.zlib.net/>`_,
`zstandard <https://facebook.github.io/zstd/>`_ and
`zfp <https://zfp.readthedocs.io/en/release1.0.0/>`_. In order to use
compression, the HDF5 library must have been built with SZIP or zlib, or one
may load 3rd party filter plugins.

Compression Performance
-----------------------

.. list-table::
   :header-rows: 1

   * - Algorithm
     - Settings
     - Compression Ratio
     - Built in to HDF5
   * - ZFP Lossy
     - Fixed-Rate: 4.0
     - 6.10
     - N
   * - ZFP Lossy
     - Fixed-Accuracy: 0.00001
     - 1.95
     - N
   * - ZFP Lossless
     - N/A
     - 1.94
     - N
   * - SZIP
     - ``pixels_per_block = 16``
     - 1.55
     - Y
   * - ZStandard
     - Level 22
     - 1.32
     - N
   * - zlib
     - Level 6
     - 1.21
     - Y
