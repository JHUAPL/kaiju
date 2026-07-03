
Installing the Intel Fortran compiler
=====================================

----

Introduction
------------

These instructions are valid for all systems which use a Linux-style interface.

Install the compiler
--------------------

Download and install the Intel Fortran compiler as described `here <https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran>`_).

Add the newly-installed compiler to your command path as follows:

If you use ``sh``, ``bash``, ``zsh``, or a compatible shell:

.. code-block:: shell

   export PATH=/opt/intel/oneapi/compiler/latest/mac/bin/intel64:${PATH}


**Note**: Your installation location may differ from the default location shown above. If so, adjust the ``PATH`` setting appropriately.

**Note**: If using MacOS Big Sur or later, you may need to install ``xcode`` (in its entirety, not just the command-line tools). The ``xcode`` package is available in the MacOS App Store.
