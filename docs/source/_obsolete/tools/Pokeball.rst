Project pokeball information
============================

Motivation
----------

we want to use Containers to make the Kaiju source code distribution that
includes all dependencies, easy running the containerized kaiju, we
furthermore want to host interactive tutorials using
`Jupyter Hub <https://jupyter.org/>`_

General Setup
-------------

(1) There is a dedicated branch in bitbucket ``pokeball/clean``, any push to
which will trigger the build of a container.

(2) After the container is built, a series of tests could (!) be run.

(3) After the tests, the container image is pushed to a (private) registry
(i.e. a storage/distribution place for OCI compliant images). The image will
be tagged with two handles:

(a) the SHA of the build:

.. image:: https://bitbucket.org/repo/kMoBzBp/images/2637046703-Screenshot%202022-11-07%20at%2018.15.22.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/2637046703-Screenshot%202022-11-07%20at%2018.15.22.png
   :alt: Screenshot 2022-11-07 at 18.15.22.png

(b) "latest" if it is the latest image:

(4) After the push, the image can now be retrieved ("pulled") by anyone who
has the credentials. To make testing user-friendly, currently each push to
bitbucket will trigger the deployment of the newest kaiju-image into
Jupiter-hub, such that it can be tried out immediately afterwards, without
having to modify kubernetes**.

.. image:: https://bitbucket.org/repo/kMoBzBp/images/2354825168-Screenshot%202022-11-07%20at%2014.34.00.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/2354825168-Screenshot%202022-11-07%20at%2014.34.00.png
   :alt: Screenshot 2022-11-07 at 14.34.00.png {width=50%}

.. image:: https://bitbucket.org/repo/kMoBzBp/images/3615671085-Screenshot%202022-11-07%20at%2014.46.27.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/3615671085-Screenshot%202022-11-07%20at%2014.46.27.png
   :alt: Screenshot 2022-11-07 at 14.46.27.png {width=50%}

Jupyter documentation
---------------------

`Documentation <https://docs.jupyter.org/en/latest/>`_ of the Jupyter project
with all the components.

`Introduction to Jupyter lab <https://nocomplexity.com/documents/jupyterlab/intro.html>`_
skip installation, userful not sure if everything works exactly the same on
our notebooks.

`Python Notebook Introduction <https://realpython.com/jupyter-notebook-introduction/>`_
you can skip the installation and the the notebook is run from withing
*jupyter lab*.

'ulimit -s'
-----------

In the terminal the 'ulimit -s' will now be automatically ``unlimited`` for
the notebook this is not used. Still looking if there is a way to have it
automatically but what can be done is:

.. code-block:: python

   import resource
   resource.setrlimit(resource.RLIMIT_STACK, (-1, -1))

in the top cell and a ``!ulimit -s`` in a following cell or execution with a
 magic ``%%bash`` will have ``unlimited`` as well.

Edit a file in jupyter
----------------------

Jupyter-lab has a file editor. It also provides a multitude of syntax
highlighting (View -> Text Editor Syntax Highlighting).



File browser
~~~~~~~~~~~~

The file browser will not allow you to navigate *above* the root directory
that it starts in. In our case this is ``/home/jovyan/``.  As a result it is
not possible to navigate to the ``/app/`` folder where the source is currently
stored.

It is still possible to edit the files in jupyter with the magic commands:

Edit via magic commands
~~~~~~~~~~~~~~~~~~~~~~~

In a notebook type

.. code-block:: python

   %open /paht/to/file

and

.. code-block:: python

   %save /path/to/file

In case of the ``/app/`` folder this will not help you much as nobody has the
rights to write there.

Path
~~~~

If the kaiju *home* directory is moved, you can use

.. code-block:: bash

   export KAIJUHOME=/path/to/new/dir
   source /app/source_kaiju

to update the path with the following:

.. code-block:: bash

   export PATH="${KAIJUHOME}/scripts/datamodel:${PATH}"
   export PATH="${KAIJUHOME}/scripts/preproc:${PATH}"
   export PATH="${KAIJUHOME}/scripts/postproc:${PATH}"
   export PATH="${KAIJUHOME}/scripts/OHelio:${PATH}"
   export PATH="${KAIJUHOME}/scripts/quicklook:${PATH}"
   export PATH="${KAIJUHOME}/scripts/legacy:${PATH}"
   export PYTHONPATH="${KAIJUHOME}:${PYTHONPATH}"

Storage
~~~~~~~

You have three directories that are for different usage:


``~/myhome/<yourid>`` is a persistent storage that should be used similar to
a HPC cluster home

\``~/group_shares/jupyter_research_pokeball_student`` is a share that all of
you have. You can exchange files here. If you want to have a folder that only
you can delete set the sticky bit: ``chmod 1750 others_can_not_delete_me``

``~/work`` this is your scratch. It will exist for some time but not forever
it has 10GB atm.

Intel version
-------------

Currently the package version ``2022.2.1`` is installed for mkl, python and
the compilers. MPI and the dev tools use ``2021.7.1``. Note that the version
of the compiler will display ``2021.7.1``.

LFortran for interactive fortran
--------------------------------

LFortran is added to the image you will see a Fortran icon next to the Python
icon in *Notebooks*. Be aware, ``Lfortran`` is not like a normal fortran in
terms of the scope, see
https://stackoverflow.com/questions/70825597/lfortran-in-a-jupyter-notebook-kills-kernel

Plotting in notebook
--------------------

As there is currently something not working with the usual ``ipywidgets`` and
``matplotlib`` and the interactive part I looked up how to do interactive with
``plotly`` and it works not to bad. An example can be found in the
``Plotly_examples.ipynb``
<https://bitbucket.org/aplkaiju/kaiju/src/2ff1f3321aa3c0c58e4b5012d26bfcda309c5951/Plotly_examples.ipynb?at=pokeball%2Fclean&viewer=nbviewer>`_.
Fyi: as the plotly stuff is ``javascript`` it can not be previewed in the
bitbucket ``nbviewer``. Therefore the stuff is slightly prolonged in the
preview.

Important is that at the beginning of the first plot or before jupyter knows
that the plots should be displayed in the notebook. This is done with
``%matplotlib inline``.

ParaView and HDF5
-----------------

The analysis notebook was extended with two new features:

HDF5 viewer

ParaView Kernel

Note: The HDF5 viewer is also in the "normal" notebook

**Note**: Extensions need to be activated to properly work. Forth symbol on
the far left.

HDF5 viewer
~~~~~~~~~~~

With `Jupyterlab-h5web <https://github.com/silx-kit/jupyterlab-h5web>`_ based
on the
`H5Web: React components for data visualization and exploration <https://h5web.panosc.eu/>`_
framework it is possible to view hdf5 files in various ways.  Depending on the
size of the file this might take a while.

ParaView Kernel
~~~~~~~~~~~~~~~

The recent development of Kitware of an
`Jupyter ParaView Kernel <https://gitlab.kitware.com/paraview/iparaview-kernel>`_
that allows interactive access to a ParaView Server out of jupyter.
Unfortunately, the list of requirements is not complete. The kernel uses
ParaView itself for building and there are additional dependencies, see the
Dockerfile for the information.

In the background a ParaView server is started (currently >5.11.0 with QT) and
one can interact with it. The
`ParaView Python <https://kitware.github.io/paraview-docs/latest/python/>`_
can also be used to interact with ParaView. A sample notebook is in the group
share
(``~/group_shares/jupyter_research_pokeball_student/ParaviewTest.ipynb``).

When opening a notebook from storage it might not have the correct kernel
attached. Make sure that the kernel in the right top corner says
``IParaView Kernel``:

.. image:: https://bitbucket.org/repo/kMoBzBp/images/3536368658-Bildschirmfoto%20vom%202023-02-14%2011-43-58.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/3536368658-Bildschirmfoto%20vom%202023-02-14%2011-43-58.png
   :alt: Bildschirmfoto vom 2023-02-14 11-43-58.png


Depending on the setup there might be an error with:

.. code-block:: bash

   QStandardPaths: XDG_RUNTIME_DIR not set, defaulting to '/tmp/runtime-'

This does not prevent work as a default is set but it is not nice to see.
Nevertheless, in the ``Dockerfile_analysis_jupyter`` file this is now set to

.. code-block:: Docker

   ENV XDG_RUNTIME_DIR="/tmp/runtime-${NB_USER}"

where

.. code-block:: Docker

   ARG NB_USER="jovyan"
