Python RADEX interface
======================

A wrapper for RADEX (www.sron.rug.nl/~vdtak/radex/) in python.

As of v0.2, created October 26, 2013, this package includes both a python
wrapper of the command-line program and a direct wrapper of the fortran code
created with f2py.

Installation procedure for the f2py-wrapped version
---------------------------------------------------

You need to have `gfortran` and `f2py` on your path.  If you've successfully
built numpy from source, you should have both.

All you need to do is:

.. code-block:: bash

   $ python setup.py install

This will call a procedure `install_radex` that downloads the latest version of
RADEX from the radex homepage, patches the source, and builds a file `radex.so`,
which is a python shared object that can be imported.  

Using the f2py-wrapped version
------------------------------

The direct wrapper of the fortran code uses a class `Radex` as its underlying
structure.  This class is useful for direct manipulation of RADEX inputs and
direct access to its outputs.

Example:
.. code-block:: python

    import pyradex
    import numpy as np
    Rlvg = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=1e16, species='co',method='lvg')
    Rlvg.run_radex()
    print Rlvg.level_population[:3],Rlvg.tex[:3],Rlvg.tau[:3],Rlvg.tex[:3]*(1-np.exp(-Rlvg.tau[:3]))
    Rslab = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=1e16, species='co',method='slab')
    Rslab.run_radex()
    print Rslab.level_population[:3],Rslab.tex[:3],Rslab.tau[:3],Rslab.tex[:3]*(1-np.exp(-Rslab.tau[:3]))
    Rsphere = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=1e16, species='co',method='sphere')
    Rsphere.run_radex()
    print Rsphere.level_population[:3],Rsphere.tex[:3],Rsphere.tau[:3],Rsphere.tex[:3]*(1-np.exp(-Rsphere.tau[:3]))

Result::
    
    [ 0.21720238  0.45362191  0.27314034] [ 15.27471017  10.86732113   8.30670325] [ 0.93769234  2.74275176  2.01021824] [  9.29419812  10.16754271   7.19394197]
    [ 0.17957399  0.39486278  0.31297916] [ 17.80769375  14.88651187  11.44840706] [ 0.68134195  1.96024231  2.03949857] [  8.79811204  12.79012934   9.95903882]
    [ 0.23532836  0.4805592   0.24340073] [ 14.38256087   9.28920338   7.50189024] [ 1.06765592  3.16666395  1.84556901] [ 9.43764227  8.89771958  6.31707599]
    
Note that because of how RADEX was written, i.e. with common blocks, the values
stored in each of these objects is identical!  You cannot have two independent
copies of the RADEX class *ever*.

Recommended installation procedure for the command-line version
---------------------------------------------------------------

1. `make` radex as normal, but create two executables: `radex_sphere`, `radex_lvg`, and `radex_slab` by
   building with one of these three lines commented out each time::

    c      parameter (method = 1)  ! uniform sphere
          parameter (method = 2)  ! expanding sphere (LVG)
    c      parameter (method = 3)  ! plane parallel slab (shock)

2. Copy these to your system path
3. `python setup.py install` to install pyradex


Simple example
--------------
Using some trivial defaults::

    In [1]: import pyradex

    In [2]: T = pyradex.radex(collider_densities={'H2':1000})
    WARNING: Assumed thermal o/p ratio since only H2 was given but collider file has o- and p- H2 [pyradex.core]

    In [3]: T.pprint(show_units=True)
    J_up J_low E_UP   FREQ      WAVE    T_EX    TAU      T_R   POP_UP POP_LOW FLUX_Kkms    FLUX_Inu
                K     GHz        um      K                K                    K km / s erg / (cm2 s)
    ---- ----- ---- -------- --------- ----- --------- ------- ------ ------- --------- -------------
       1     0  5.5 115.2712 2600.7576 5.044 0.0004447 0.00086 0.4709    0.47 0.0009155     1.806e-11

    In [4]: T.meta
    Out[4]:
    {'Column density [cm-2]': '1.000E+12',
     'Density of H2  [cm-3]': '1.000E+03',
     'Density of oH2 [cm-3]': '3.509E-04',
     'Density of pH2 [cm-3]': '1.000E+03',
     'Geometry': 'Uniform sphere',
     'Line width     [km/s]': '1.000',
     'Molecular data file': '/Users/adam/repos/Radex/data/co.dat',
     'Radex version': '20nov08',
     'T(background)     [K]': '2.730',
     'T(kin)            [K]': '10.000'}




Timing information
------------------
i.e., how fast is it?::

    %timeit T = pyradex.radex(collider_densities={'H2':1000})
    1 loops, best of 3: 149 ms per loop


    for n in 10**np.arange(6):
       %timeit T = pyradex.radex(collider_densities={'H2':n})

    10 loops, best of 3: 149 ms per loop
    10 loops, best of 3: 150 ms per loop
    10 loops, best of 3: 149 ms per loop
    10 loops, best of 3: 151 ms per loop
    10 loops, best of 3: 150 ms per loop
    10 loops, best of 3: 149 ms per loop

    for n in 10**np.arange(12,18):
       ....:     %timeit T = pyradex.radex(collider_densities={'H2':1000}, column_density=n)

    10 loops, best of 3: 149 ms per loop
    10 loops, best of 3: 149 ms per loop
    10 loops, best of 3: 149 ms per loop
    10 loops, best of 3: 150 ms per loop
    10 loops, best of 3: 152 ms per loop
    10 loops, best of 3: 157 ms per loop
    
These results indicate that, even in highly optically thick cases where more
iterations are required, the execution time is dominated by the python
overheads.

If you redo these tests comparing the fortran wrapper to the "naive" version,
the difference is enormous.  The following tests can be seen in `timing.py
<examples/timing.py>`__:

::

    Python:  0.892609834671
    Fortran:  0.0151958465576
    py/fortran:  58.7403822016
    Python:  0.902825832367
    Fortran:  0.0102920532227
    py/fortran:  87.7206727205
    Python:  0.876524925232
    Fortran:  0.0730140209198
    py/fortran:  12.0048850096
    Python:  0.836034059525
    Fortran:  0.0925290584564
    py/fortran:  9.03536762906
    Python:  0.880390882492
    Fortran:  0.0725519657135
    py/fortran:  12.1346248008
    Python:  0.96048283577
    Fortran:  0.0753719806671
    py/fortran:  12.7432346512
    
    

Making Grids
------------
Is more efficient with the other script, but you can still do it...  ::

    for n in 10**np.arange(12,18):
        T = pyradex.radex(collider_densities={'H2':1000}, column_density=n)
        T.pprint()
    
    Row# Line# E_UP   FREQ      WAVE    T_EX    TAU      T_R   POP_UP POP_LOW FLUX_Kkms  FLUX_Inu
    ---- ----- ---- -------- --------- ----- --------- ------- ------ ------- --------- ---------
       1     0  5.5 115.2712 2600.7576 5.044 0.0004447 0.00086 0.4709    0.47 0.0009155 1.806e-11
    Row# Line# E_UP   FREQ      WAVE    T_EX   TAU      T_R    POP_UP POP_LOW FLUX_Kkms  FLUX_Inu
    ---- ----- ---- -------- --------- ----- -------- -------- ------ ------- --------- ---------
       1     0  5.5 115.2712 2600.7576 5.047 0.004444 0.008589  0.471  0.4698  0.009143 1.803e-10
    Row# Line# E_UP   FREQ      WAVE    T_EX   TAU     T_R   POP_UP POP_LOW FLUX_Kkms  FLUX_Inu
    ---- ----- ---- -------- --------- ----- ------- ------- ------ ------- --------- ---------
       1     0  5.5 115.2712 2600.7576 5.075 0.04415 0.08473 0.4721  0.4681    0.0902 1.779e-09
    Row# Line# E_UP   FREQ      WAVE    T_EX  TAU    T_R   POP_UP POP_LOW FLUX_Kkms  FLUX_Inu
    ---- ----- ---- -------- --------- ----- ------ ------ ------ ------- --------- ---------
       1     0  5.5 115.2712 2600.7576 5.336 0.4152 0.7475 0.4817  0.4527    0.7957 1.569e-08
    Row# Line# E_UP   FREQ      WAVE    T_EX  TAU  T_R  POP_UP POP_LOW FLUX_Kkms  FLUX_Inu
    ---- ----- ---- -------- --------- ----- ----- ---- ------ ------- --------- ---------
       1     0  5.5 115.2712 2600.7576 6.929 2.927 3.49 0.5057  0.3745     3.715 7.327e-08
    Row# Line# E_UP   FREQ      WAVE    T_EX  TAU  T_R  POP_UP POP_LOW FLUX_Kkms  FLUX_Inu
    ---- ----- ---- -------- --------- ----- ----- ---- ------ ------- --------- ---------
       1     0  5.5 115.2712 2600.7576 9.294 18.09 5.96 0.4696  0.2839     6.345 1.252e-07

If you want to create a grid with the directly wrapped version, do loops with
constant temperature: every time you load a new temperature, RADEX must read in
the molecular data file and interpolate across the collision rate values, which
may be a substantial overhead.

.. image:: https://d2weczhvl823v0.cloudfront.net/keflavich/pyradex/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

