import pyradex
import pyradex.fjdu
import pytest
import os
import distutils.spawn
import numpy as np
from astropy import units as u

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

def test_thin_co():
    density = {'oH2': 100,
               'pH2': 900,
              }
    RR = pyradex.Radex(species='co', column=1e10, density=density, temperature=20)
    FF = pyradex.fjdu.Fjdu(species='co', column=1e10, density=density, temperature=20)

    rtbl = RR()
    ftbl = FF()

    diff = rtbl['upperlevelpop'] - ftbl['upperlevelpop']
    #np.testing.assert_allclose(diff, 0)
