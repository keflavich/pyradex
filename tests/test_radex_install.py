import numpy as np

def test_radex_install():
    from pyradex.radex import radex

    assert radex.setup.radat.dtype == np.dtype('S1')
