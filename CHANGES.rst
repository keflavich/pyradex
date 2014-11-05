Changes
=======

0.4.0 - November 1, 2014: Add Fujun Du's "myradex" with an identical interface
        to RADEX.  Factored out common functions to base_class.py

0.3.1 - October 28, 2014: allow 'validate_colliders' to be disabled so that
        fast grids can be run.  Fast grids require validate_colliders=False,
        reload_molfile=False.
      - Major refactor of quantity/units use: factor of ~10 speedup for default
        R.__call__ approach.

0.3.0 - August 28, 2014: MAJOR internal refactor and API change.  `astropy`
        is now required, and there is no longer a 'total h2 column' concept
        associated with the RADEX code.  For compatibility with DESPOTIC, there
        will need to be an additional interface layer.
        Changed default species too, since `hco+` doesn't exist on the LAMDA
        servers.

0.2.2 - August 27, 2014: Minor bugfixes including update of specutils version
        inclusion

0.2.1 - May 29, 2014 - bugfix release.  First to include CHANGES.  Intended to
        work when pip-installing

0.2 - first python-wrapped version?
