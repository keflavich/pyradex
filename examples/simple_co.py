import pyradex
import pylab as pl

R = pyradex.Radex(column=1e16)
R.maxiter=1000
#print R.radex.impex.molfile,R.molpath

for temperature in [5,10,20,30,40,50,60]:
    R.temperature = temperature
    
    R.run_radex()
    
    pl.figure(1)
    pl.plot(R.level_population, label="T=%i" % R.temperature)

    pl.figure(2)
    pl.plot(R.frequency, R.tau, label="T=%i" % R.temperature)

    pl.figure(3)
    pl.plot(R.frequency, R.tex, label="T=%i" % R.temperature)

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
