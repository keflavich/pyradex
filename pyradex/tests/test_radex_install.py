import numpy as np

def test_radex_install():
    from pyradex.radex import radex

    # older f2py exposes this as an array of chars (S1), newer as one string (S120)
    assert radex.setup.radat.dtype.kind == 'S'
