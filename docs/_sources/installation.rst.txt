.. _installation:

Installation of QUIP, quippy and GAP
************************************

These instructions provide more details on the compilation and
installation of ``QUIP`` (Fortran library and main programs) and
and ``GAP`` (Fortran add-on to QUIP), and the ``quippy`` (Python interface).

Precompiled Containers
----------------------

If you have access to `Docker <https://hub.docker.com>`_ or
`Singularity <http://singularity.lbl.gov>`_, you can try one of the
`precompiled images <https://github.com/libAtoms/quip-docker>`_
to get up and running quickly.

Compilation Instructions
------------------------

First try the quickstart below, which should work with most Linux systems.

Quick start
^^^^^^^^^^^

Install [#]_ the prerequisites: GCC, gfortran, Python, and the linear algebra
libraries.  For example, on Ubuntu, do (in a terminal):

::

    $ sudo apt-get install gcc gfortran python python-pip libblas-dev liblapack-dev

For other systems, replace the ``apt-get`` part with your system package manager.
Beware that the packages might also have slightly different names; these can
usually be found with a quick search.

Don't forget the ``quippy`` prerequisites:

::

    $ pip install numpy ase f90wrap

Now you can get the code and compile:

::

    $ git clone --recursive https://github.com/libAtoms/QUIP.git
    $ export QUIP_ARCH=linux_x86_64_gfortran
    $ export QUIPPY_INSTALL_OPTS=--user  # omit for a system-wide installation
    $ make config

Answer all the questions with their defaults (by pressing enter) for now, just
to get things working.

::

    $ make
    $ make install-quippy

And now open a Python terminal and see if it works:

::

    $ python
    >>> import quippy
    >>>

If the import completes successfully (i.e. with no output) then the
installation was successful.  You may want to continue with `Installing the
Jupyter notebook`_ to run the interactive tutorials.

.. [#] If this isn't your machine and you don't have root access, these
   packages might already be installed by the system administrator.  If not,
   ask them.


Step by step
^^^^^^^^^^^^

If that didn't work, try these step-by-step instructions
instructions excerpted from the top-level `README
<https://github.com/libAtoms/QUIP/blob/public/README.md>`_.  The ``README`` file
is the most up-to-date source of installation information.

  .. include:: ../../../README.md
    :start-after: Compilation Instructions
    :end-before: ### Mac OS

If that still doesn't work or you're using a nonstandard architecture, try
looking at `Custom settings`_ and `Common Problems`_.  As a last resort you can
consult the `issue tracker on Github`_.


Installing the Jupyter notebook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Jupyter`_ is an environment for interactive computing that makes using Python
much easier and more intuitive.  Especially useful is its notebook environment,
which provides a handy way to experiment with code, see the results, and have a
record of your progress.  The interactive getting-started tutorial is a Jupyter
notebook that you can run and modify yourself.

To get Jupyter up and running, the following should suffice [#]_:

::

    $ pip install jupyter
    $ jupyter notebook

This will open a new window in your browser that you can use to navigate
through your filesystem.  To access the interactive tutorials, you can run the
``jupyter notebook`` command from your ``QUIP/src/GAP/doc/Examples`` directory (or any
enclosing directory) then navigate to the notebooks and open
``Introduction.ipynb`` to get started.

.. [#] This assumes you've already run ``sudo apt-get install python-pip; pip
   install numpy; pip install ase`` as in the `Quick start`_.


Custom settings
---------------

:makevar:`MATHS_LINKOPTS`
   Library options needed to link to BLAS and LAPACK libraries. Any working
   BLAS/LAPACK installation is fine. If you are using Linux, ATLAS is
   a good option, and you should use something like the following::

     -L/usr/local/atlas -llapack -lf77blas -lcblas -latlas

   On Mac OS X, there are build in LAPACK libraries in the Accelerate
   framework, which you can use by entering

     -framework Accelerate

:makevar:`FOX_LIBDIR`, :makevar:`FOX_INCDIR` and :makevar:`FOX_LIBS`
  Directories containing FoX libraries and header files, and required link options.
  Should be read automatically from QUIP Makefiles.

:makevar:`QUIPPY_INSTALL_OPTS`
   Installation options, e.g. specify ``--user`` to install for the current
   user ``--prefix=${PREFIX}`` to install in a non-default location.

:makevar:`HAVE_NETCDF4`
  Should be set to 1 to enable NetCDF4 support. Should be read automatically from QUIP.

:makevar:`NETCDF4_LIBS`, :makevar:`NETCDF4_FLAGS`
  Linker flags for compiling with NetCDF4 support, and flags for finding
  header files. Should be read automatically from QUIP.


.. _install_faq:

Common Problems
---------------

Permission errors when installing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are installing as root, you may need to make sure the value of
the :envvar:`QUIP_ARCH` gets through to the install script, e.g. ::

   sudo QUIP_ARCH=darwin_x86_64_gfortran make install-quippy

ImportError when importing
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get an :exc:`ImportError` with a message about unresolved
dependancies then something went wrong with the linking process -
check that all the libraries you're linking against are correct. You
can used `ldd` on Linux of `otool -L` on Mac OS X to check which
libraries the :file:`_quippy.so` Python extension is linked against.

Possible problems installing atomeye module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get an :exc:`ImportError` with a message ::
   >>> import atomeye
   ImportError: dlopen(/Users/silvia/lib/python/_atomeye.so, 2): Symbol not found: _Config_load_libatoms
   Referenced from: /Users/silvia/lib/python/_atomeye.so
   Expected in: flat namespace
   in /Users/silvia/lib/python/_atomeye.so

be sure that you have set :envvar:`QUIP_ROOT` variable before starting the compilation.
If not make clean and recompile again

If you get an :exc:`ImportError` with a message ::
   >>> import atomeye
   ImportError: dlopen(/Users/silvia/lib/python/_atomeye.so, 2): Symbol not found: __gfortran_adjustl
   Referenced from: /Users/silvia/lib/python/_atomeye.so
   Expected in: flat namespace
   in /Users/silvia/lib/python/_atomeye.so

be sure that the gfortran libraries are properly set in :makevar:`ATOMEYE_LIBS` in Makefile.atomeye

Error compiling IPModel_GAP
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you get the following error during compilation::

   /src/Potentials/IPModel_GAP.f95:51.22:

   use descriptors_module
                         1
   Fatal Error: Can't open module file 'descriptors_module.mod' for reading at (1): No such file or directory

The `GAP` module is not publicly available, so the
:file:`Makefile.inc` must contain :makevar:`HAVE_GAP` = 1.

.. _`issue tracker on Github`: https://github.com/libAtoms/QUIP/issues
.. _`Jupyter`: http://jupyter.org/
