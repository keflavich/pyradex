import pyradex
import numpy as np
import astropy.units as u

def test_synthspec():
    R = pyradex.Radex(column=1e15)
    R.run_radex()
    wcs = np.linspace(103.0, 103.1, 1000)*u.GHz
    S = pyradex.synthspec.SyntheticSpectrum.from_RADEX(wcs, R,
                                                       linewidth=10*u.km/u.s)

