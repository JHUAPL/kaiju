Computational Costs
===================

GAMERA-HELIO+GL model
---------------------

Estimates of computational costs in cpu-hours (NCAR) or SBUs (NASA) for
different model resolutions based on Helio runs

On the Derecho supercomputer:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

AMD Epyc Milan (128 cores)
~~~~~~~~~~~~~~~~~~~~~~~~~~


256x128x256 - 64 ranks on 8 nodes
    160 Hour Spinup + 100 Hrs CME propagation - ~1000 cpu-hours

512x256x512 - 256 ranks on 32 nodes
    160 Hour Spinup + 50 Hrs CME propagation - ~8-10,000 cpu-hours

1024x512x1024  - 512 ranks on 64 nodes
    160 Hour Spinup + 50 hours of CME propagation time - ~85k cpu-hours

2048x1024x2048 - 1024 ranks on 128 nodes
    160 Hour Spinup ~490k cpu-hours, 32 hours of CME propagation
    ~125k cpu-hours

On the Pleiades supercomputer:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cascade Lake Nodes (40 cores, 80 threads) Using HT and AVX-512
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

256x128x256 - 64 ranks on 8 nodes
    160 Hour Spinup + 100 Hrs CME propagation - ~200 SBUs

512x256x512 - 256 ranks on 16 nodes
    160 Hour Spinup + 100 Hrs CME propagation - ~3200 SBUs

1024x512x1024  - 512 ranks on 32 nodes
    160 Hour Spinup + 50 hours of CME propagation time - ~52k SBUs

Broadwell Nodes (28 cores, 28 threads)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

256x128x256 - 64 ranks on 16 nodes
    160 Hour Spinup + 100 Hrs CME propagation - ~300 SBUs

512x256x512 - 256 ranks on 64 nodes
    160 Hour Spinup + 100 Hrs CME propagation - ~4000 SBUs Broadwell Nodes
