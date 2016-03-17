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

You need to clone this repository first with `--recursive` enabled so that `myRadex <https://github.com/fjdu/myRadex>`_ is downloaded:

.. code-block:: bash

   git clone --recursive https://github.com/keflavich/pyradex.git

Then `cd` to the source directory and run:

.. code-block:: bash

   $ python setup.py install_radex install_myradex install

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
    R = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=1e16, species='co', temperature=20)
    Tlvg = R(escapeProbGeom='lvg')
    Tslab = R(escapeProbGeom='slab')
    Tsphere = R(escapeProbGeom='sphere')
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

Examples
--------
There is a rich examples gallery.  We have a few notebooks:

    http://nbviewer.ipython.org/github/keflavich/pyradex/blob/master/examples/pH2CO_interactive.ipynb
    http://nbviewer.ipython.org/github/keflavich/pyradex/blob/master/examples/FittingTheGrid.ipynb
    http://nbviewer.ipython.org/github/keflavich/pyradex/blob/master/examples/Interactive.ipynb
    http://nbviewer.ipython.org/github/keflavich/pyradex/blob/master/examples/oH2CO-interactive.ipynb
    http://nbviewer.ipython.org/github/keflavich/pyradex/blob/master/examples/pH2CO_interactive.ipynb
    http://nbviewer.ipython.org/github/keflavich/pyradex/blob/master/examples/ph2co_interactive_mm.ipynb

and a series of more involved examples:

 * examples/ch3cn_110_synthspec.py
 * examples/h2co_grids.py
 * examples/h2cs_thermometer.py
 * examples/interactive_setup_mm.py
 * examples/oh2co_density_grid.py
 * examples/oh2co_distributions.py
 * examples/oh2co_grids_2.py
 * examples/ph2co_grid_computation.py
 * examples/ph2co_grid_computation_mm.py
 * examples/ph2co_grids.py
 * examples/ph2co_grids_2.py
 * examples/ph2co_required_sn.py
 * examples/simple_co.py
 * examples/simple_co_column.py
 * examples/synthspec_ch3cn.py
 * examples/timing.py

Most of these were written to make sensitivity estimates for observing proposals.

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

    %timeit T = pyradex.pyradex(collider_densities={'H2':1000})
    10 loops, best of 3: 31.8 ms per loop

    for n in 10**np.arange(6):
       %timeit T = pyradex.pyradex(collider_densities={'H2':n})

    10 loops, best of 3: 32.1 ms per loop
    10 loops, best of 3: 32.5 ms per loop
    10 loops, best of 3: 32 ms per loop
    10 loops, best of 3: 32.1 ms per loop
    10 loops, best of 3: 32.4 ms per loop
    10 loops, best of 3: 31.9 ms per loop

    for n in 10**np.arange(12,18):
        %timeit T = pyradex.pyradex(collider_densities={'H2':1000}, column=n)

    10 loops, best of 3: 31.8 ms per loop
    10 loops, best of 3: 32.2 ms per loop
    10 loops, best of 3: 32.5 ms per loop
    10 loops, best of 3: 32.2 ms per loop
    10 loops, best of 3: 32.7 ms per loop
    10 loops, best of 3: 33.1 ms per loop
    

If you redo these tests comparing the fortran wrapper to the "naive" version,
the difference can be enormous.  The following tests can be seen in `timing.py
<examples/timing.py>`__:

::


    Python external call:               0.0323288917542
    Fortran-wrapped:                    0.0183672904968
    Fortran-wrapped, no reload:         0.000818204879761
    Fortran-wrapped, no reload, reuse:  0.000756096839905
    Fortran (call method):  0.0270668029785
    py/fortran:                    1.76013395986
    py/fortran, __call__ method:   1.1944111678
    py/fortran, no reload:         39.5119762224
    py/fortran, no reload, reuse:  42.7576072904
    Python external call:               0.0332223176956
    Fortran-wrapped:                    0.0169018030167
    Fortran-wrapped, no reload:         0.000811815261841
    Fortran-wrapped, no reload, reuse:  0.000753211975098
    Fortran (call method):  0.0275466918945
    py/fortran:                    1.96560790957
    py/fortran, __call__ method:   1.20603656594
    py/fortran, no reload:         40.9234948605
    py/fortran, no reload, reuse:  44.1075272221
    Python external call:               0.0312483787537
    Fortran-wrapped:                    0.0216565847397
    Fortran-wrapped, no reload:         0.00535380840302
    Fortran-wrapped, no reload, reuse:  0.000751805305481
    Fortran (call method):  0.031253194809
    py/fortran:                    1.44290427735
    py/fortran, __call__ method:   0.999845901985
    py/fortran, no reload:         5.83666362361
    py/fortran, no reload, reuse:  41.5644562839
    Python external call:               0.0316061973572
    Fortran-wrapped:                    0.0228497028351
    Fortran-wrapped, no reload:         0.00549430847168
    Fortran-wrapped, no reload, reuse:  0.000753903388977
    Fortran (call method):  0.031331205368
    py/fortran:                    1.38322137427
    py/fortran, __call__ method:   1.00877693615
    py/fortran, no reload:         5.75253419427
    py/fortran, no reload, reuse:  41.9234053319
    Python external call:               0.0318208932877
    Fortran-wrapped:                    0.0216773033142
    Fortran-wrapped, no reload:         0.00544350147247
    Fortran-wrapped, no reload, reuse:  0.000751280784607
    Fortran (call method):  0.0315539121628
    py/fortran:                    1.46793597093
    py/fortran, __call__ method:   1.0084611101
    py/fortran, no reload:         5.84566633234
    py/fortran, no reload, reuse:  42.3555266415
    Python external call:               0.0322543859482
    Fortran-wrapped:                    0.0225975990295
    Fortran-wrapped, no reload:         0.00569999217987
    Fortran-wrapped, no reload, reuse:  0.00075900554657
    Fortran (call method):  0.0314954996109
    py/fortran:                    1.42733685583
    py/fortran, __call__ method:   1.02409507221
    py/fortran, no reload:         5.65867196486
    py/fortran, no reload, reuse:  42.4955866185
    [ 0.006951  0.006911  0.006956]
    [ 0.006951  0.006911  0.006956]
    [ 0.006951  0.006911  0.006956]
    pyradex.pyradex timing for a 3^4 grid:  [2.6063590049743652, 2.598068952560425, 2.592205047607422]
    [ 0.00694859  0.00690934  0.00695345]
    [ 0.00694859  0.00690934  0.00695345]
    [ 0.00694859  0.00690934  0.00695345]
    pyradex.Radex() timing for a 3^4 grid:  [3.8620870113372803, 3.838628053665161, 3.805685043334961]
    [ 0.00694859  0.00690934  0.00695345]
    [ 0.00694859  0.00690934  0.00695345]
    [ 0.00694859  0.00690934  0.00695345]
    pyradex.Radex() class-based timing for a 3^4 grid:  [3.1014058589935303, 3.2805678844451904, 3.160888195037842]
    [ 0.00694859  0.00690934  0.00695345]
    [ 0.00694859  0.00690934  0.00695345]
    [ 0.00694859  0.00690934  0.00695345]
    pyradex.Radex() class-based timing for a 3^4 grid, using optimal parameter-setting order:  [0.9963750839233398, 1.0024840831756592, 0.9699358940124512]
    

Making Grids
------------
Is more efficient with scripts, but you can still do it...  ::

    R = pyradex.Radex(species='co', collider_densities={'H2':1000}, column=1e15)
    for n in 10**np.arange(12,18):
        T = R(collider_densities={'H2':1000}, column=n)
        T[:1].pprint()
    
             Tex             tau         frequency  upperstateenergy upperlevel lowerlevel upperlevelpop  lowerlevelpop       brightness           T_B
          K                             GHz            K                                                             erg / (cm2 Hz s sr)        K
    ------------- ----------------- ----------- ---------------- ---------- ---------- -------------- -------------- ------------------- ----------------
    11.0274813968 0.000166783361591 115.2712018             5.53          1          0 0.540537331305 0.297561763825   5.20877418593e-18 0.00127591598469
         Tex            tau         frequency  upperstateenergy upperlevel lowerlevel upperlevelpop  lowerlevelpop       brightness           T_B
          K                            GHz            K                                                             erg / (cm2 Hz s sr)        K
    ------------- ---------------- ----------- ---------------- ---------- ---------- -------------- -------------- ------------------- ---------------
    11.0274813968 0.00166783361591 115.2712018             5.53          1          0 0.540537331305 0.297561763825    5.2048669339e-17 0.0127495888324
         Tex            tau        frequency  upperstateenergy upperlevel lowerlevel upperlevelpop  lowerlevelpop       brightness          T_B
          K                           GHz            K                                                             erg / (cm2 Hz s sr)       K
    ------------- --------------- ----------- ---------------- ---------- ---------- -------------- -------------- ------------------- --------------
    10.9980972475 0.0166790919823 115.2712018             5.53          1          0 0.538730147174 0.296964688622   5.14681095066e-16 0.126073777202
         Tex           tau        frequency  upperstateenergy upperlevel lowerlevel upperlevelpop  lowerlevelpop       brightness          T_B
          K                          GHz            K                                                             erg / (cm2 Hz s sr)       K
    ------------- -------------- ----------- ---------------- ---------- ---------- -------------- -------------- ------------------- -------------
    11.7797140751 0.150601068675 115.2712018             5.53          1          0 0.530489509066 0.282823341198   4.78772386104e-15 1.17277754545
         Tex           tau        frequency  upperstateenergy upperlevel lowerlevel upperlevelpop  lowerlevelpop       brightness          T_B
          K                          GHz            K                                                             erg / (cm2 Hz s sr)       K
    ------------- -------------- ----------- ---------------- ---------- ---------- -------------- -------------- ------------------- -------------
    15.0692631019 0.955344506002 115.2712018             5.53          1          0 0.454752879863 0.218821739485   2.92170292028e-14 7.15686133711
         Tex           tau       frequency  upperstateenergy upperlevel lowerlevel upperlevelpop  lowerlevelpop       brightness          T_B
          K                         GHz            K                                                             erg / (cm2 Hz s sr)       K
    ------------- ------------- ----------- ---------------- ---------- ---------- -------------- -------------- ------------------- -------------
    22.6356250741 4.17742617995 115.2712018             5.53          1          0 0.318586709967 0.135596426565   7.69430015071e-14 18.8475833332

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
dependent variable and all the others can be changed.  Since abundance is not
an explicit input into RADEX, this is done with some property machinery behind
the scenes.  In v0.3, the length in Radex has been fixed to 1 pc.
    

.. image:: https://d2weczhvl823v0.cloudfront.net/keflavich/pyradex/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

.. _install: install.rst

