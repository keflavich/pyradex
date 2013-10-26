import timeit
import numpy as np
import textwrap

for n in 10**np.arange(12,18):

    setup = "import pyradex"
    ptiming = timeit.Timer(stmt="pyradex.pyradex(collider_densities={'oH2':900,'pH2':100},column_density=%e)" % n,setup=setup).repeat(3,10)
    print "Python: ",np.min(ptiming)
    setup = """
    import pyradex
    R = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=%e)""" % n
    ftiming = timeit.Timer(stmt="R.run_radex(); T = R.tex",setup=textwrap.dedent(setup)).repeat(3,10)
    print "Fortran: ",np.min(ftiming)
    print "py/fortran: ",np.min(ptiming)/np.min(ftiming)
