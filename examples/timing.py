import timeit
import numpy as np
import textwrap

for n in 10**np.arange(12,18):

    setup = "import pyradex"
    ptiming = timeit.Timer(stmt="pyradex.pyradex(collider_densities={'oH2':900,'pH2':100},column=%e)" % n,setup=setup).repeat(3,10)
    print "Python: ",np.min(ptiming)
    setup = """
    import pyradex
    R = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=%e)""" % n
    ftiming = timeit.Timer(stmt="R.run_radex(); T = R.tex",setup=textwrap.dedent(setup)).repeat(3,10)
    print "Fortran: ",np.min(ftiming)
    #dominated by array creation...
    #ftiming2 = timeit.Timer(stmt="R(collider_densities={'oH2':900,'pH2':100}, column=%e)" % n, setup=textwrap.dedent(setup)).repeat(3,10)
    #print "Fortran (call method): ",np.min(ftiming2)
    print "py/fortran: ",np.min(ptiming)/np.min(ftiming)#,"py/fortran (call method): ",np.min(ptiming)/np.min(ftiming2)


gridtest = """
# build a small grid
for ii,T in enumerate([5,10,20]):
    for jj,column in enumerate([1e13,1e15,1e17]):
        for kk,density in enumerate([1e3,1e5,1e7]):
            for mm,opr in enumerate([1e-2,0.1,1]):
                fortho = opr/(1+opr)
                print density,fortho,column,T
                result = {caller}(collider_densities={{'H2':density,'oH2':density*fortho,'pH2':density*(1-fortho)}},
                                  column=column,
                                  temperature=T,
                                  species='co',
                                  minfreq=100,
                                  maxfreq=400)
                grid[ii,jj,kk,mm] = result['TAU'][0]
"""

setup = """
import pyradex
import numpy as np
grid = np.empty([3,3,3,3])
"""
ptiming = timeit.Timer(stmt=gridtest.format(caller='pyradex.pyradex'),setup=setup).repeat(3,1)
print "pyradex.pyradex timing for a 3^4 grid: ",ptiming
setup += "R = pyradex.Radex()"
ftiming = timeit.Timer(stmt=gridtest.format(caller='R'),setup=setup).repeat(3,1)
print "pyradex.Radex() timing for a 3^4 grid: ",ftiming
