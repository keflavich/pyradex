Python RADEX interface
======================

A wrapper for RADEX (www.sron.rug.nl/~vdtak/radex/) in python.


Recommended installation procedure:

1. `make` radex as normal, but create two executables: `radex_sphere`, `radex_lvg`, and `radex_slab` by
   building with one of these three lines commented out each time::

    c      parameter (method = 1)  ! uniform sphere
          parameter (method = 2)  ! expanding sphere (LVG)
    c      parameter (method = 3)  ! plane parallel slab (shock)

2. Copy these to your system path
3. `python setup.py install` to install pyradex


.. TODO:: Write a wrapper script to build the different RADEX versions
