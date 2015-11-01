# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# not sure if this is needed: it's not a normal part of package-template from .version import __version__
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from .core import pyradex,write_input,parse_outfile,call_radex,Radex

    from . import utils
    from . import despotic_interface
    from . import radex
    from . import synthspec
