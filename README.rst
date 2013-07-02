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
    

Making Grids
------------
Is more efficient with the other script, but you can still do it...

::
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

