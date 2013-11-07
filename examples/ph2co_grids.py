"""
Create some simple grids for the low-frequency para-H2CO lines
"""
import pyradex
import pylab as pl
import numpy as np
import matplotlib

ntemp,ndens = 20,20

temperatures = np.linspace(10,50,ntemp)
densities = np.logspace(2.5,7,ndens)
abundance = 10**-8.5
opr = 0.01 # assume primarily para
fortho = opr/(1+opr)

taugrid_71M = np.empty([ndens,ntemp])
texgrid_71M = np.empty([ndens,ntemp])
fluxgrid_71M = np.empty([ndens,ntemp])
taugrid_145 = np.empty([ndens,ntemp])
texgrid_145 = np.empty([ndens,ntemp])
fluxgrid_145 = np.empty([ndens,ntemp])
taugrid_355M = np.empty([ndens,ntemp])
texgrid_355M = np.empty([ndens,ntemp])
fluxgrid_355M = np.empty([ndens,ntemp])
columngrid = np.empty([ndens,ntemp])

import os
if not os.path.exists('ph2co-h2.dat'):
    import urllib
    urllib.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/ph2co-h2.dat')

R = pyradex.Radex(species='ph2co-h2', abundance=abundance)
R.run_radex()

# get the table so we can look at the frequency grid
table = R.get_table()

# Target frequencies:
table[np.array([6,1,11])].pprint()

for ii,tt in enumerate(temperatures):
    R.temperature = tt
    for jj,dd in enumerate(densities):
        R.density = {'oH2':dd*fortho,'pH2':dd*(1-fortho)}
        R.abundance = abundance # reset column to the appropriate value
        R.run_radex(reuse_last=False, reload_molfile=True)

        TI = R.total_intensity
        taugrid_71M[jj,ii] = R.tau[6]
        texgrid_71M[jj,ii] = R.tex[6].value
        fluxgrid_71M[jj,ii] = TI[6].value
        taugrid_145[jj,ii] = R.tau[1]
        texgrid_145[jj,ii] = R.tex[1].value
        fluxgrid_145[jj,ii] = TI[1].value
        taugrid_355M[jj,ii] = R.tau[11]
        texgrid_355M[jj,ii] = R.tex[11].value
        fluxgrid_355M[jj,ii] = TI[11].value
        columngrid[jj,ii] = R.column.value

pl.figure(1)
pl.clf()
extent = [densities.min(),densities.max(),temperatures.min(),temperatures.max()]
for kk,(grid,freq,label) in enumerate(zip([taugrid_145,taugrid_355M,texgrid_145,texgrid_355M],['145','355M','145','355M'],[r'\tau',r'\tau','T_{ex}','T_{ex}'])):
    ax = pl.subplot(2,2,kk+1)
    #ax.imshow(grid, extent=extent)
    pc = ax.pcolormesh(np.log10(densities),temperatures,grid)
    ax.set_title('$%s$ p-H$_2$CO %s GHz' % (label,freq))
    ax.set_xlabel('log Density')
    ax.set_ylabel('Temperature')
    pl.colorbar(pc)

pl.figure(2)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,taugrid_145/taugrid_355M, vmax=1e4, vmin=0)
pl.colorbar(cax)
ax.set_title('$\\tau$ p-H$_2$CO 145 GHz/355MHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,texgrid_145/texgrid_355M)#, vmax=1.3, vmin=0.8)
pl.colorbar(cax)
ax.set_title('$T_{ex}$ p-H$_2$CO 145 GHz/355 MHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')

pl.figure(3)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,taugrid_71M/taugrid_355M, vmax=45, vmin=-2)
pl.colorbar(cax)
ax.set_title('$\\tau$ p-H$_2$CO 71/355 MHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,texgrid_71M/texgrid_355M, vmax=2, vmin=0.0)
pl.colorbar(cax)
ax.set_title('$T_{ex}$ p-H$_2$CO 71/355 MHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')


pl.figure(4)
pl.clf()
extent = [densities.min(),densities.max(),temperatures.min(),temperatures.max()]
for kk,(freq,name) in enumerate(zip([0.071,145,0.355],['71M','145','355M'])):
    ax = pl.subplot(2,2,kk+1)
    grid = eval('fluxgrid_%s' % name)
    #ax.imshow(grid, extent=extent)
    pc = ax.pcolormesh(np.log10(densities),temperatures,grid)
    ax.set_title('Flux p-H$_2$CO %i GHz' % (freq))
    ax.set_xlabel('log Density')
    ax.set_ylabel('Temperature')
    pl.colorbar(pc)

pl.figure(6)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,fluxgrid_71M/fluxgrid_355M)#, vmax=0.01, vmin=0)
pl.colorbar(cax)
ax.set_title('$S_{\\nu}$ p-H$_2$CO 71/355 MHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,fluxgrid_145/fluxgrid_355M)#, vmax=1.2, vmin=0.8)
pl.colorbar(cax)
ax.set_title('$S_{\\nu}$ p-H$_2$CO 145 GHz/355 MHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')


pl.show()

