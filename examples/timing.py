import timeit
import numpy as np
import textwrap

for n in 10**np.arange(12,18):

    setup = "import pyradex"
    nreps = 10
    ptiming = timeit.Timer(stmt="pyradex.pyradex(collider_densities={'oH2':900,'pH2':100},column=%e,temperature=20)" % n,setup=setup).repeat(3,nreps)
    print "Python: ",np.min(ptiming)/nreps
    setup = """
    import pyradex
    R = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=%e, temperature=20)""" % n
    ftiming = timeit.Timer(stmt="R.run_radex(); T = R.tex",setup=textwrap.dedent(setup)).repeat(3,nreps)
    print "Fortran-wrapped: ",np.min(ftiming)/nreps
    ftiming3 = timeit.Timer(stmt="R.run_radex(validate_colliders=False, reload_molfile=False); T = R.tex",setup=textwrap.dedent(setup)).repeat(3,nreps)
    print "Fortran-wrapped, no reload: ",np.min(ftiming3)/nreps
    #dominated by array creation...
    #ftiming2 = timeit.Timer(stmt="R(collider_densities={'oH2':900,'pH2':100}, column=%e)" % n, setup=textwrap.dedent(setup)).repeat(3,10)
    #print "Fortran (call method): ",np.min(ftiming2)
    print "py/fortran: ",np.min(ptiming)/np.min(ftiming)#,"py/fortran (call method): ",np.min(ptiming)/np.min(ftiming2)
    print "py/fortran, no reload: ",np.min(ptiming)/np.min(ftiming3)#,"py/fortran (call method): ",np.min(ptiming)/np.min(ftiming2)

    # Can check correctness too:
    # x = pyradex.pyradex(collider_densities={'oH2':900,'pH2':100},column=n, temperature=20)
    # print x['pop_up'][0], x['t_ex'][0]
    # R.column=n; R.temperature=20
    # R.run_radex(reload_molfile=False, validate_colliders=False)
    # print R.level_population[1], R.tex[0]
    # R.run_radex()
    # print R.level_population[1], R.tex[0]


gridtest = """
# build a small grid
for ii,T in enumerate([5,10,20]):
    for jj,column in enumerate([1e13,1e15,1e17]):
        for kk,density in enumerate([1e3,1e5,1e7]):
            for mm,opr in enumerate([1e-2,0.1,1]):
                fortho = opr/(1+opr)
                result = {caller}(collider_densities={{'oH2':density*fortho,'pH2':density*(1-fortho)}},
                                  column=column,
                                  temperature=T,
                                  species='co')
                grid[ii,jj,kk,mm] = result['tau'][0]
"""

setup = """
import pyradex
import numpy as np
grid = np.empty([3,3,3,3])
"""
ptiming = timeit.Timer(stmt=gridtest.format(caller='pyradex.pyradex'),setup=setup).repeat(3,1)
print "pyradex.pyradex timing for a 3^4 grid: ",ptiming
setup += "R = pyradex.Radex(column=1e15)"
ftiming = timeit.Timer(stmt=gridtest.format(caller='R'),setup=setup).repeat(3,1)
print "pyradex.Radex() timing for a 3^4 grid: ",ftiming

gridtest_class = """
# build a small grid
for ii,T in enumerate([5,10,20]):
    for jj,column in enumerate([1e13,1e15,1e17]):
        for kk,density in enumerate([1e3,1e5,1e7]):
            for mm,opr in enumerate([1e-2,0.1,1]):
                fortho = opr/(1+opr)
                R.density = {'oH2':density*fortho,'pH2':density*(1-fortho)}
                R.column=column
                R.temperature=T
                grid[ii,jj,kk,mm] = R.tau[0]
"""

ftiming2 = timeit.Timer(stmt=gridtest_class,setup=setup).repeat(3,1)
print "pyradex.Radex() class-based timing for a 3^4 grid: ",ftiming2
