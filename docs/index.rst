PyRADEX
=======


``pyradex`` is a python-wrapped version of RADEX_, using f2py to make grid
building *nearly* as fast as the Fortran code, but much more convenient.

``pyradex`` also wraps Mark Krumholz's DESPOTIC_ code with an interface
identical to pyradex so that the codes can be directly compared.  However, I
have achieved no success yet in matching the code results!  I have some
:doc:`notes` about this process, but they are messy and incomplete.

Using pyradex
-------------

pyradex is centered around the `~pyradex.core.Radex` class.  To initialize it, a
valid ``.dat`` file must be available in the appropriate `directory <directory_setup>`_.

Then, initialize the `~pyradex.core.Radex` class by giving it initial
conditions and a chemical species:

    >>> import pyradex
    >>> rdx = pyradex.Radex(species='hco+', collider_densities={'H2':1e3},
    ...                     column=1e13)

To get values from the object, you need to run RADEX:

    >>> rdx.run_radex()
    >>> print(rdx.tex)
    >>> tbl = rdx.get_table()


.. _directory_setup:

Directory Configuration
-----------------------


.. _RADEX: http://home.strw.leidenuniv.nl/~moldata/radex.html
.. _DESPOTIC: https://sites.google.com/a/ucsc.edu/krumholz/codes/despotic
