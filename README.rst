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

   $ python setup.py install_radex install

This will call a procedure `install_radex` that downloads the latest version of
RADEX from the radex homepage, patches the source, and builds a file `radex.so`,
which is a python shared object that can be imported.  

See the install_ page for more details.

If you want pyradex to look in a specific directory for the molecular data
files, you can specify an environmental variable `RADEX_DATAPATH` prior to
starting python.  It can also be specified interactively with the `datapath`
keyword.


Using the f2py-wrapped version
------------------------------

The direct wrapper of the fortran code uses a class `Radex` as its underlying
structure.  This class is useful for direct manipulation of RADEX inputs and
direct access to its outputs.

Example:

.. code-block:: python

    import pyradex
    import numpy as np
    R = pyradex.Radex()
    Tlvg = R(collider_densities={'oH2':900,'pH2':100}, column=1e16, species='co',method='lvg')
    Tslab = R(collider_densities={'oH2':900,'pH2':100}, column=1e16, species='co',method='slab')
    Tsphere = R(collider_densities={'oH2':900,'pH2':100}, column=1e16, species='co',method='sphere')
    Tlvg[:3].pprint()
    Tslab[:3].pprint()
    Tsphere[:3].pprint()

Result::
    
         Tex           tau        frequency  upperstateenergy upperlevel lowerlevel  upperlevelpop    lowerlevelpop         flux
    ------------- -------------- ----------- ---------------- ---------- ---------- ---------------- --------------- -----------------
    15.2747101724 0.937692338925 115.2712018             5.53          2          1   0.273140336953  0.453621905471 2.93964536078e-14
    10.8673211326  2.74275175782     230.538             16.6          3          2  0.0518618367484  0.273140336953 9.26125039465e-14
    8.30670325364  2.01021823976 345.7959899            33.19          4          3 0.00379591658449 0.0518618367484 8.16324298598e-14
         Tex           tau        frequency  upperstateenergy upperlevel lowerlevel  upperlevelpop   lowerlevelpop         flux
    ------------- -------------- ----------- ---------------- ---------- ---------- ---------------- -------------- -----------------
    17.8076937528 0.681341951256 115.2712018             5.53          2          1   0.312979158313 0.394862780876 2.89304678735e-14
    14.8865118666  1.96024230849     230.538             16.6          3          2   0.102821702575 0.312979158313 1.38012283784e-13
     11.448407058  2.03949857132 345.7959899            33.19          4          3 0.00920322307626 0.102821702575  1.6139902821e-13
         Tex           tau       frequency  upperstateenergy upperlevel lowerlevel  upperlevelpop   lowerlevelpop         flux
    ------------- ------------- ----------- ---------------- ---------- ---------- ---------------- -------------- -----------------
      14.38256087 1.06765591906 115.2712018             5.53          2          1   0.243400727834 0.480559204909 2.93394133644e-14
    9.28920337666  3.1666639484     230.538             16.6          3          2   0.037299201561 0.243400727834 7.24810556601e-14
    7.50189023571 1.84556901411 345.7959899            33.19          4          3 0.00307839203073 0.037299201561 6.19215196139e-14

    
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
Is more efficient with scripts, but you can still do it...  ::

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

If you want to build a grid, *do not* make an astropy table each time!  That
appears to dominate the overhead at each iteration.

A note on self-consistency in LVG calculations
----------------------------------------------

LVG computations have weird units.  The opacity of a line only depends on the
velocity-coherent column along the line of sight, i.e. the column per km/s.

The key assumption in the LVG Sobolev approximation is that each "cell" can be
treated independently such that there are no nonlocal radiative effects.

This independence implies that there is a separation between the local volume
density and the total line-of-sight column density.

However, the quantities reported by typical codes - RADEX, DESPOTIC - are
integrated line-of-sight values.  The column density, abundance, and local
volume density are not independent, then.

In order to have a self-consistent cloud (or line of sight), you must assume
some length scale.  Usually, one specifies a velocity gradient per length scale
rather than an absolute length scale, but the length scale is important.

If a total column density of hydrogen `N(H)` is specified along with a density
`n(H)`, the length scale is trivial: `N(H)/n(H) = L`.  If you increase the
density, this length scale decreases - so far all is fine.

Within RADEX, the standard free variable is the column of the molecule of
interest.  
If you change the column of the molecule, which is possible to do explicitly,
and hold everything else fixed in RADEX (`n(H)`, `dV`), the change can be
interpreted as a change in the size scale or the column.

One could consider the alternative possibility of treating the length scale as
a free parameter, but this approach contains a danger of changing the
interpretation of the processes involved: if the length scale is decreased for
a fixed delta-V, the velocity gradient `dv/dl` must be larger.  This
interpretation should be avoided as it bears the risk of breaking the LVG
assumption.  The velocity gradient is also often an imposed constraint via the
observed linewidth, while the length scale is only weakly constrained in most
situations.

In DESPOTIC, the free variables are the total column density, the density,
the abundance, and the velocity gradient.  Length is therefore left as the
dependent variable, consistent with the above.

The Classes (`Despotic` & `Radex`) are constructed such that length is a
dependent variables and all the others can be changed.  Since abundance is not
an explicit input into RADEX, this is done with some property machinery behind
the scenes.
    

.. image:: https://d2weczhvl823v0.cloudfront.net/keflavich/pyradex/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

.. _install: install.rst
