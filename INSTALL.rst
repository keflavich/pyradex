Installation
------------

Installation *should* be relatively easy, but it may break down depending on
your compiler.

Installation requires astroquery_.  It can be pip installed:

.. code-block:: bash

   pip install astroquery

Compiling RADEX
~~~~~~~~~~~~~~~

Compiling the fortran code requires ``gfortran`` and, on numpy>=1.26
(the default on any currently-supported Python), the ``meson`` and
``ninja`` build tools that numpy's f2py uses to compile extensions:

.. code-block:: bash

   pip install meson ninja

This command should install the python-wrapped-fortran version of RADEX:

.. code-block:: bash

   $ python setup.py install_radex install_myradex install

If it doesn't work, there are a number of ways things could have gone wrong.

If the ``install_radex`` step fails with ``FileNotFoundError: ... 'meson'``,
you're missing the build tools mentioned above; install them and re-run.

If the ``install_myradex`` step fails to compile ``opkda1.f`` with an error
like ``Type mismatch in argument 'iwk' at (1); passed REAL(8) to
INTEGER(4)``, your ``gfortran`` is version 10 or newer, which rejects
argument-type mismatches that this legacy Fortran 77 code relies on.  This
is fixed in recent versions of the myRadex_ submodule (which adds
``-fallow-argument-mismatch``); run ``git submodule update --remote
myRadex`` to pick up the fix if your checkout predates it.

If compilation fails with ``ld: -rpath can only be used when targeting Mac
OS X 10.5 or later`` (seen with conda-forge's ``gfortran`` on very recent
macOS versions that the compiler doesn't recognize, e.g. macOS 26), set an
explicit deployment target before building:

.. code-block:: bash

   export FFLAGS='-mmacosx-version-min=11.0'
   export LDFLAGS="$LDFLAGS -mmacosx-version-min=11.0"

If you try to `import pyradex` and see the error::

    ImportError: No module named radex 
   
that means that the "shared object" file `radex.so` (which is what python
actually imports) never got built.  In this case, you should have seen an
error message at the `python setup.py install_radex` step.  At this point,
look in the `radex_build.log` file - it may contain useful information.

I don't know how general this advice is, but for installation on 2 macs and 1
linux machine, the trick was using the right `gfortran` versions with the right
flags.

If you're using mac-native compilers (e.g., gcc-4.2 from Xcode), you can use
the `-arch` flags.  I used this really long command to make sure everything got
set right:

.. code-block:: bash

    FFLAGS='-arch i686 -arch x86_64 -fPIC' CFLAGS='-fno-strict-aliasing -fno-common -dynamic -arch i386 -arch x86_64 -g -O2' LDFLAGS='-arch i686 -arch x86_64 -undefined dynamic_lookup -bundle' python setup.py install_radex

If you're using a different compiler, say from hpc.sourceforge.net, you need a different
set of flags, and you need to make sure that `gcc --version` and `gfortran --version` match.

.. code-block:: bash

    FFLAGS='-m64 -fPIC' CFLAGS='-fno-strict-aliasing -fno-common -dynamic -m64 -g -O2' LDFLAGS='-m64 -undefined dynamic_lookup -bundle' python setup.py install_radex install

Note that if you're using a 32 bit machine or 32 bit python, you should use
`-m32` instead of `-m64` in those flags.

You may also need to set your C_INCLUDE_PATH, e.g.:

.. code-block:: bash

    C_INCLUDE_PATH /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/include

(I didn't have to do this, but someone else did - the error that pointed us in this direction was that `stdio.h` couldn't be found)

For linux, the build failed with gcc-4.0.1 but succeeded with gcc-4.3.6, which
suggests that gcc>=4.2 might be required (since gcc-4.2 worked on a mac).


Environment Setup
~~~~~~~~~~~~~~~~~

The shell environmental variable ``RADEX_DATAPATH`` should be set to the
directory in which you store collision rate data files.  You can also use the
``datapath`` keyword argument to specify the path during initialization of the
RADEX object.


Installation Troubles
~~~~~~~~~~~~~~~~~~~~~

If you still cannot install pyradex given the above information, please post on
the issues_ page, and include the following information:

 * python version
 * numpy version
 * astropy version
 * gfortran version

e.g.:

.. code-block:: bash

    $ python -c "import sys, astropy, numpy; print(sys.version); print(numpy.__version__,astropy.__version__)"
    2.7.2 (v2.7.2:8527427914a2, Jun 11 2011, 15:22:34)
    [GCC 4.2.1 (Apple Inc. build 5666) (dot 3)]
    ('1.6.1', '0.3.dev6331')

    $ gfortran --version
    GNU Fortran (GCC) 4.2.3
    Copyright (C) 2007 Free Software Foundation, Inc.
   

.. _issues: https://github.com/keflavich/pyradex/issues

.. _astroquery: astroquery.readthedocs.org
.. _myRadex: https://github.com/keflavich/myRadex
