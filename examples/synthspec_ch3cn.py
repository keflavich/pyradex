"""
Example:
    Plot the relative intensities of the 3mm CH3CN (methyl cyanide) lines
"""
import pyradex
import pylab as pl
from astropy import units as u
import numpy as np
import matplotlib as mpl

R = pyradex.Radex(species='ch3cn',abundance=1e-11,column=None)

# There are 5 lines of interest in this band
nlines = 5
fluxes = {ii:[] for ii in xrange(nlines)}

# Temperature range: 20-500 K is allowed (by CH3CN data file)
temperatures = np.linspace(20,500,8)
temperatures = [20,50,100,300,500]

# set up figure
pl.figure(1)
pl.clf()

for ii,temperature in enumerate(temperatures):

    R.temperature = temperature
    R.run_radex()
    #wcs = pyradex.synthspec.FrequencyArray(91.95*u.GHz,92*u.GHz,1000)
    wcs = np.linspace(91.95,92,1000)*u.GHz
    S = pyradex.synthspec.SyntheticSpectrum.from_RADEX(wcs,R,linewidth=5*u.km/u.s)#R.get_table(),'ch3cn')
    
    # spectral colors
    color = mpl.cm.spectral(float(ii)/(len(temperatures)))
    S.plot(label='%i K' % temperature, color=color, linewidth=2, alpha=0.9)

    for ii in xrange(nlines):
        fluxes[ii].append(S.table[ii]['T_B'])

pl.legend(loc='best')
pl.savefig("CH3CN_synthetic_spectra.pdf",bbox_inches='tight')

linenames = {ii:S.table[ii]['upperlevel']+" - "+S.table[ii]['lowerlevel'] for ii in xrange(nlines)}

pl.figure(2)
pl.clf()
pl.subplot(1,2,1)
for ii in xrange(nlines):
    pl.plot(temperatures,fluxes[ii],label=linenames[ii])

# Line #4 is the "reference line" at lowest energy
pl.subplot(1,2,2)
for ii in xrange(nlines):
    pl.plot(temperatures,np.array(fluxes[ii])/np.array(fluxes[4]),label=linenames[ii])

pl.savefig("CH3CN_flux_ratios.pdf",bbox_inches='tight')

pl.figure(3)
pl.clf()
for ii in xrange(nlines):
    pl.plot(temperatures,np.array(fluxes[4])/np.array(fluxes[ii]),label=linenames[ii])

pl.legend(loc='best')
