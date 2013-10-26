import pyradex
import pylab as pl
import numpy as np

R = pyradex.Radex(column=1e16,temperature=20)
R.maxiter=1000
#print R.radex.impex.molfile,R.molpath

for column in np.linspace(1e14,1e17,10):
    R.column = column

    #R.debug = 1

    R.run_radex(reuse_last=True)
    
    pl.figure(1)
    pl.plot(R.level_population, label="c=%e" % column)

    pl.figure(2)
    pl.plot(R.frequency, R.tau, label="c=%e" % column)

    pl.figure(3)
    pl.plot(R.frequency, R.tex, label="c=%e" % column)

f1 = pl.figure(1)
ax1 = f1.gca()
ax1.set_xlim(0,10)
ax1.set_xlabel("Energy Level")
ax1.set_ylabel("Population")

f2 = pl.figure(2)
ax2 = f2.gca()
ax2.set_xlabel("Frequency")
ax2.set_ylabel("Optical Depth")
ax2.set_xlim(0,1000)

f3 = pl.figure(3)
ax3 = f3.gca()
ax3.set_xlabel("Frequency")
ax3.set_ylabel("Excitation Temperature")
ax3.axis([0,1000,0,65])

print("Same as above, but not starting from last position.")
for column in np.linspace(1e14,1e17,10):
    R.column = column

    R.run_radex(reuse_last=False)
