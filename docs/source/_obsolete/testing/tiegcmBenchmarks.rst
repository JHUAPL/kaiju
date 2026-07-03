Instructions for running TIEGCM benchmarks on Cheyenne
======================================================

Introduction
------------

These instructions are intended to be an update and augment to the standard
`TIEGCM Benchmark instructions <https://www.hao.ucar.edu/modeling/tgcm/tiegcm2.0/userguide/html/benchmarks.html>`_
that are bit dated and relate the yellowstone computer system.  Using
repository tag TBD will get you the updates for the scripting system for
cheyenne.

In order to produce the graphics you will need to install the tgcmprocf90 code
base that utilizes NCAR graphics to produce plots.  This needs to be run on
casper with specified set of modules.

Benchmark runs
--------------

The TIEGCM benchmarks fall into three categories


#. Full year climatologies (``run_climatology``)
#. Seasons (``run_seasons``)
#. Storms (``run_storms``)

for each category there is a tcsh script that creates run directories,
complies the code, and submits the jobs at 5 and 2.5 degree resolutions. These
scripts, located in the ``benchmarks`` subdirectory of the code distribution,
use the tgcmrun program which is written in python 2 to create the batch
scripts.

Compiling and running the code on cheyenne has been tested with the following
set of modules

.. code-block:: bash

   Currently Loaded Modules:
     1) git/2.22.0     5) ncarcompilers/0.5.0   9) mpt/2.25
     2) ncarenv/1.3    6) mkl/2021.2           10) esmf_libs/8.3.0
     3) cmake/3.14.4   7) conda/latest         11) esmf-8.3.0-ncdfio-mpt-O
     4) intel/2021.2   8) netcdf/4.9.0

Visualizing the results
-----------------------

As previously mentioned to visualize the results we are currently utilizing an
older NCAR graphics package so the work will need to be conducted on casper
using a specific module and conda environment.  We are working on python
utilities and will have that available shortly.

On casper you will need to have the following set of modules loaded.

.. code-block:: bash

   Currently Loaded Modules:
     1) ncarenv/1.2    6) conda/latest         11) esmf-8.3.0-ncdfio-mpiuni-O
     2) ffmpeg/4.1.3   7) intel/19.1.1         12) ncl/6.6.2
     3) git/2.22.0     8) ncarcompilers/0.5.0  13) nco/5.1.4
     4) netcdf/4.9.1   9) openmpi/4.1.1
     5) cmake/3.14.4  10) esmf_libs/8.3.0

and your conda environment can be created from this YML file

.. code-block:: yaml

   name: tgcmproc
   channels:
     - conda-forge
     - ncar
     - r
   dependencies:
     - _libgcc_mutex=0.1=conda_forge
     - _openmp_mutex=4.5=2_gnu
     - backports=1.0=pyhd8ed1ab_3
     - backports.functools_lru_cache=1.6.1=py_0
     - backports_abc=0.5=py_1
     - bzip2=1.0.8=h7f98852_4
     - c-ares=1.19.1=hd590300_0
     - ca-certificates=2023.5.7=hbcca054_0
     - certifi=2019.11.28=py27h8c360ce_1
     - cftime=1.1.1=py27h588f082_0
     - curl=7.87.0=h6312ad2_0
     - cycler=0.10.0=py_2
     - dbus=1.13.6=hfdff14a_1
     - enum34=1.1.10=py27h8c360ce_1
     - expat=2.5.0=hcb278e6_1
     - fontconfig=2.14.2=h14ed4e7_0
     - freetype=2.12.1=hca18f0e_1
     - functools32=3.2.3.2=py_3
     - futures=3.3.0=py27h8c360ce_1
     - gettext=0.21.1=h27087fc_0
     - glib=2.66.3=h58526e2_0
     - gst-plugins-base=1.14.5=h0935bb2_2
     - gstreamer=1.14.5=h36ae1b5_2
     - hdf4=4.2.15=h9772cbc_5
     - hdf5=1.10.5=nompi_h7c3c948_1111
     - icu=64.2=he1b5a44_1
     - jpeg=9e=h0b41bf4_3
     - keyutils=1.6.1=h166bdaf_0
     - kiwisolver=1.1.0=py27h9e3301b_1
     - krb5=1.20.1=hf9c8cef_0
     - ld_impl_linux-64=2.40=h41732ed_0
     - libblas=3.9.0=8_openblas
     - libcblas=3.9.0=8_openblas
     - libclang=9.0.1=default_hb4e5071_5
     - libcurl=7.87.0=h6312ad2_0
     - libedit=3.1.20191231=he28a2e2_2
     - libev=4.33=h516909a_1
     - libexpat=2.5.0=hcb278e6_1
     - libffi=3.2.1=he1b5a44_1007
     - libgcc-ng=12.2.0=h65d4601_19
     - libgfortran-ng=7.5.0=h14aa051_20
     - libgfortran4=7.5.0=h14aa051_20
     - libglib=2.66.3=hbe7bbb4_0
     - libgomp=12.2.0=h65d4601_19
     - libiconv=1.17=h166bdaf_0
     - liblapack=3.9.0=8_openblas
     - libllvm9=9.0.1=default_hc23dcda_7
     - libnetcdf=4.7.3=nompi_h9f9fd6a_101
     - libnghttp2=1.51.0=hdcd2b5c_0
     - libopenblas=0.3.12=pthreads_hb3c22a3_1
     - libpng=1.6.39=h753d276_0
     - libsqlite=3.42.0=h2797004_0
     - libssh2=1.10.0=haa6b8db_3
     - libstdcxx-ng=12.2.0=h46fd767_19
     - libuuid=2.38.1=h0b41bf4_0
     - libxcb=1.15=h0b41bf4_0
     - libxkbcommon=0.10.0=he1b5a44_0
     - libxml2=2.9.10=hee79883_0
     - libzlib=1.2.13=h166bdaf_4
     - matplotlib=2.2.5=ha770c72_3
     - matplotlib-base=2.2.5=py27h250f245_1
     - ncurses=6.3=h27087fc_1
     - netcdf4=1.5.3=nompi_py27hd35fb8e_102
     - nspr=4.35=h27087fc_0
     - nss=3.89=he45b914_0
     - numpy=1.16.5=py27h95a1406_0
     - openssl=1.1.1t=h0b41bf4_0
     - pcre=8.45=h9c3ff4c_0
     - pip=20.1.1=pyh9f0ad1d_0
     - pthread-stubs=0.4=h36c2ea0_1001
     - pyparsing=2.4.7=pyh9f0ad1d_0
     - pyqt=5.12.3=py27hcca6a23_1
     - python=2.7.15=h5a48372_1011_cpython
     - python-dateutil=2.8.1=py_0
     - python_abi=2.7=1_cp27mu
     - pytz=2020.1=pyh9f0ad1d_0
     - qt=5.12.5=hd8c4c69_1
     - readline=8.2=h8228510_1
     - scipy=1.2.1=py27h921218d_2
     - setuptools=44.0.0=py27_0
     - singledispatch=3.6.1=pyh44b312d_0
     - six=1.16.0=pyh6c4a22f_0
     - sqlite=3.42.0=h2c6b66d_0
     - subprocess32=3.5.4=py27h516909a_0
     - tk=8.6.12=h27826a3_0
     - tornado=5.1.1=py27h14c3975_1000
     - wheel=0.37.1=pyhd8ed1ab_0
     - xorg-libxau=1.0.11=hd590300_0
     - xorg-libxdmcp=1.1.3=h7f98852_0
     - xz=5.2.6=h166bdaf_0
     - zlib=1.2.13=h166bdaf_4
     - pip:
         - pyqt5-sip==4.19.18
         - pyqtwebengine==5.12.1
