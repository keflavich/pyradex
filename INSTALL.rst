Installation
------------

Installation *should* be relatively easy, but it may break down depending on
your compiler.

Installation requires astroquery_ and specutils_.  These can be pip installed:

.. code-block:: bash

   pip install astroquery
   pip install https://github.com/astropy/specutils/archive/master.zip

Compiling RADEX
~~~~~~~~~~~~~~~

This command should install the python-wrapped-fortran version of RADEX:

.. code-block:: bash

   $ python setup.py install_radex install_myradex install

If it doesn't work, there are a number of ways things could have gone wrong.

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
.. _specutils: https://github.com/astropy/specutils
