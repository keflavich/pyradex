import timeit
import numpy as np
import textwrap

# Check correctness before doing timing tests
import pyradex
py_pop = [pyradex.pyradex(collider_densities={'oH2':900,'pH2':100},column=n, temperature=20)['pop_up'][0]
          for n in 10**np.arange(12,18)]

R = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=1e15, temperature=20)

R_noreload_pop = []
for n in 10**np.arange(12,18):
    R.column = n
    R.run_radex(reload_molfile=False, validate_colliders=False)
    R_noreload_pop.append(R.level_population[1])

R_pop = []
for n in 10**np.arange(12,18):
    R.column = n
    R.run_radex(reload_molfile=True, validate_colliders=True)
    R_pop.append(R.level_population[1])

R_reuse_pop = []
for n in 10**np.arange(12,18):
    R.column = n
    R.run_radex(reload_molfile=False, validate_colliders=False, reuse_last=True)
    R_reuse_pop.append(R.level_population[1])

for p1,p2,p3,p4 in zip(py_pop, R_noreload_pop, R_pop, R_reuse_pop):
    np.testing.assert_almost_equal(p1, p2, decimal=4)
    np.testing.assert_almost_equal(p2, p3)
    np.testing.assert_almost_equal(p3, p4)



for n in 10**np.arange(12,18):

    setup = "import pyradex"
    nreps = 10
    ptiming = timeit.Timer(stmt="pyradex.pyradex(collider_densities={'oH2':900,'pH2':100},column=%e,temperature=20)" % n,setup=setup).repeat(3,nreps)
    print "Python external call:              ",np.min(ptiming)/nreps
    setup = """
    import pyradex
    R = pyradex.Radex(collider_densities={'oH2':900,'pH2':100}, column=%e, temperature=20)""" % n
    ftiming = timeit.Timer(stmt="R.run_radex(); T = R.tex",setup=textwrap.dedent(setup)).repeat(3,nreps)
    print "Fortran-wrapped:                   ",np.min(ftiming)/nreps
    ftiming3 = timeit.Timer(stmt="R.run_radex(validate_colliders=False, reload_molfile=False); T = R.tex",setup=textwrap.dedent(setup)).repeat(3,nreps)
    print "Fortran-wrapped, no reload:        ",np.min(ftiming3)/nreps
    ftiming4 = timeit.Timer(stmt="R.run_radex(validate_colliders=False, reload_molfile=False, reuse_last=True); T = R.tex",setup=textwrap.dedent(setup)).repeat(3,nreps)
    print "Fortran-wrapped, no reload, reuse: ",np.min(ftiming4)/nreps
    # dominated by array creation...
    ftiming2 = timeit.Timer(stmt="R(collider_densities={'oH2':900,'pH2':100}, column=%e)" % n, setup=textwrap.dedent(setup)).repeat(3,nreps)
    print "Fortran (call method): ",np.min(ftiming2)/nreps
    print "py/fortran:                   ",np.min(ptiming)/np.min(ftiming)
    print "py/fortran, __call__ method:  ",np.min(ptiming)/np.min(ftiming2)
    print "py/fortran, no reload:        ",np.min(ptiming)/np.min(ftiming3)
    print "py/fortran, no reload, reuse: ",np.min(ptiming)/np.min(ftiming4)

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
print grid[0,0,0,:]
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
                R.run_radex(validate_colliders=False,
                            reload_molfile=False,
                            reuse_last=True)
                grid[ii,jj,kk,mm] = R.tau[0]
print grid[0,0,0,:]
"""

ftiming2 = timeit.Timer(stmt=gridtest_class,setup=setup).repeat(3,1)
print "pyradex.Radex() class-based timing for a 3^4 grid: ",ftiming2
