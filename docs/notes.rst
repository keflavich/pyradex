A comparison of DESPOTIC and RADEX
----------------------------------

This is a collection of notes on the internal implementation of the escape
probability approximation within RADEX and DESPOTIC.  It's mostly a collection
of 'notes to self' to help me understand how the individual parameters are
defined.

Math versions of despotic and radex:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Eqn 39 of Krumholz 2013 vs Eqn 21 of van der Tak 2007:

.. math::

   \tau_D = \tau_{ij,LVG} = \frac{g_i}{g_j} \frac{A_{ij} \lambda_{ij}^3}{8 \pi |dv_r/dr|} x_s N_H f_j \left(1-\frac{f_i g_j}{f_j g_i}\right)

   \tau_R = \tau = \frac{c^3}{8\pi \nu_{ul}^3} \frac{A_{ul} N_{mol}}{1.064 \Delta V} \left[ x_l \frac{g_u}{g_l} - x_u \right]

Footnote 6 in RADEX, combined with footnote 11 in DESPOTIC:

.. math::

    \tau_R = 1.5 \tau_D   (spherical)

    \tau_R = \tau_D (LVG, slab?)

Assuming the two are equal, we get the constants dropping out, and are left with

.. math::

   \frac{x_s N_H}{|dv_r/dr|} = \frac{N_{mol}}{1.064 \Delta V}

which approximately defines what is required to get equivalent results out of
DESPOTIC and RADEX.


RADEX:
""""""

The input parameters are :math:`cdmol`, :math:`deltav`, :math:`dens`.  The
density only affects the collision rates; it only affects the optical depth
through :math:`crate` and the implicit coupling to column density via the
abundance.

.. math::
   xt = xnu^3 = (f/c)^3

   cdmol \equiv  X_{mol} n(H_2) 

   X_{mol} = cdmol / n(H_2)

   cddv = cdmol / deltav = \frac{N_{mol}}{dV}

   dVdR_{RADEX} = deltav = FWHM

    \tau_{RADEX} = \frac{cddv}{(1.0645 (8.0 \pi))} A_{ul} (c/\nu)^{3} (n_l \frac{g_u}{g_l}  - n_u)

                 = \frac{1}{1.0645 (8.0 \pi) dV_{RADEX}} A_{ul} (c/\nu)^{3} (n_l \frac{g_u}{g_l}  - n_u) N_{mol}

    \tau_h = \tau_{RADEX}/2

    \beta = 2.0 \frac{1.0 - e^{-2.34 \tau_h}}{4.68 \tau_h}

          = \frac{1.0 - e^{-1.17 \tau_{RADEX}}}{1.17 \tau_{RADEX}}

We can factor the 2.35 into the tau equation to get something a little more self-consistent

.. math::

     \tau_{RADEX,f} = \frac{2.35}{2 \cdot 1.0645 (8.0 \pi) dV_{RADEX}} A_{ul} (c/\nu)^{3} (n_l \frac{g_u}{g_l}  - n_u) N_{mol}

     \tau_{RADEX,f} = \frac{8 \log(2)}{2 \sqrt(2 \pi) (8.0 \pi) dV_{RADEX}} A_{ul} (c/\nu)^{3} (n_l \frac{g_u}{g_l}  - n_u) N_{mol}

DESPOTIC:
"""""""""


.. math::

   dV = thisCloud.dVdr = \frac{dV}{dR} = \sigma

   n \cdot dR = N

   X_{H2} = 2 X_H

   n_H = n_{H2} / 2

    \tau_{DESPOTIC} = \frac{1}{(8 \pi dV_D/dR)} \frac{g_u}{g_l} A_{ul} (c/\nu)^3 \left(1-\frac{n_u g_l}{n_l g_u}\right) X_{H2} n_H n_l

                    = \frac{dR}{(8 \pi dV_D)}  A_{ul} (c/\nu)^3 \left(n_l\frac{g_u}{g_l}-n_u\right) \frac{n_{mol}}{2}

                    = \frac{1}{(8 \pi dV_D)}  A_{ul} (c/\nu)^3 \left(n_l\frac{g_u}{g_l}-n_u\right) \frac{N_{mol}}{2}


    \beta = \frac{1.0 - e^{-\tau_D}}{\tau_D}


RADEX and DESPOTIC are equivalent if

.. math::
   \tau(RADEX,f) = \tau(DESPOTIC) \frac{2.35}{1.0645}

but so far it looks like


If `dV` is given as a fixed quantity to both, in order to make them equivalent, we can do:

.. math::
   \tau(RADEX,f) = \tau_N \frac{2.35}{2 \cdot 1.0645}

   \tau(DESPOTIC) = \tau_N \frac{1}{2}

   \tau(RADEX,f)/\tau(DESPOTIC) = \frac{2.35}{1.0645}

To fix it, free the dV variable

.. math::

   \tau(RADEX,f)/\tau(DESPOTIC) = 1
        = \frac{2.35}{2 \cdot 1.0645 dV_{RADEX}} /  \frac{1}{2 dVdr_D}

   dVdr_D = \frac{1.0645}{2.35} dV_{RADEX}


One factor of 2 comes from defining a cloud *radius* vs a cloud *diameter*.
Because DESPOTIC was written do determine line transfer from the center of a
cloud outward, it uses the radius.  RADEX uses the diameter.  However, DESPOTIC
and RADEX also define their abundances differently

1.0645 is the ratio of the integral of a Gaussian to the FWHM, i.e.

.. math::
   \int_{-\infty}^{\infty}\frac{1}{\sqrt{8 \log(2)}} e^{-x^2/2} dx = \sqrt{\frac{2 \pi}{8 \log(2)}} = 1.0644670194312262

Code versions:
~~~~~~~~~~~~~~

DESPOTIC:
"""""""""


.. code-block:: python

        elif escapeProbGeom == "LVG":
            tau = (self.data.levWgt[u]/self.data.levWgt[l]) * \
                self.data.EinsteinA[u,l] * (c/self.data.freq[u,l])**3 / \
                (8.0*np.pi*abs(thisCloud.dVdr)) * \
                self.abundance * thisCloud.nH * self.levPop[l] * \
                (1.0 - self.levPop[u]*self.data.levWgt[l] / \
                     (self.levPop[l]*self.data.levWgt[u]))
            # Note that, in computing escape probabilities, we need to
            # handle the case of very small tau with care to avoid
            # roundoff problems and make sure we correctly limit to
            # beta -> 1 as tau -> 0
            idx = tau > 1e-6
            self.escapeProb[u[idx], l[idx]] = \
                (1.0 - np.exp(-tau[idx])) / tau[idx]
            idx = tau <= 1e-6
            self.escapeProb[u[idx], l[idx]] = \
                1.0 - tau[idx]/2.0

RADEX:
""""""


.. code-block:: fortran


            xt  = xnu(iline)**3.0
            m   = iupp(iline)
            n   = ilow(iline)

            taul(iline) = cddv*(xpop(n)*gstat(m)/gstat(n)-xpop(m))
     $           /(fgaus*xt/aeinst(iline))

            beta = escprob(taul(iline))

      taur = tau/2.0

         else if (method.eq.2) then
   c     Expanding sphere = Large Velocity Gradient (LVG) or Sobolev case.
   C     Formula from De Jong, Boland and Dalgarno (1980, A&A 91, 68)
   C     corrected by factor 2 in order to match ESCPROB(TAU=0)=1
           if (abs(taur).lt.0.01) then
             beta = 1.0
           else if(abs(taur).lt.7.0) then
             beta = 2.0*(1.0 - dexp(-2.34*taur))/(4.68*taur)
           else
             beta = 2.0/(taur*4.0*(sqrt(log(taur/sqrt(pi)))))
           endif

