"""
Create some simple grids for the low-frequency para-H2CO lines
"""
import pylab as pl
import numpy as np

if not 'texgrid_355M' in locals():
    execfile('ph2co_grid_computation.py')

pl.figure(1)
pl.clf()
extent = [densities.min(),densities.max(),temperatures.min(),temperatures.max()]
for kk,(grid,freq,label) in enumerate(zip([taugrid_145,taugrid_355M,texgrid_145,texgrid_355M],['145','355M','145','355M'],[r'\tau',r'\tau','T_{ex}','T_{ex}'])):
    ax = pl.subplot(2,2,kk+1)
    #ax.imshow(grid, extent=extent)
    pc = ax.pcolormesh(temperatures,np.log10(densities),grid)
    ax.set_title('$%s$ p-H$_2$CO %s GHz' % (label,freq))
    ax.set_ylabel('log Density')
    ax.set_xlabel('Temperature')
    pl.colorbar(pc)

pl.figure(7)
pl.clf()
extent = [densities.min(),densities.max(),temperatures.min(),temperatures.max()]
for kk,(grid,freq,label) in enumerate(zip([taugrid_71M,taugrid_355M,texgrid_71M,texgrid_355M],['71M','355M','71M','355M'],[r'\tau',r'\tau','T_{ex}','T_{ex}'])):
    ax = pl.subplot(2,2,kk+1)
    #ax.imshow(grid, extent=extent)
    pc = ax.pcolormesh(temperatures,np.log10(densities),grid, vmax=50 if 'ex' in label else None)
    ax.set_title('$%s$ p-H$_2$CO %sHz' % (label,freq))
    ax.set_ylabel('log Density')
    ax.set_xlabel('Temperature')
    pl.colorbar(pc)

pl.figure(2)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(temperatures,np.log10(densities),taugrid_145/taugrid_355M, vmax=1e4, vmin=0)
pl.colorbar(cax)
ax.set_title('$\\tau$ p-H$_2$CO 145 GHz/355MHz')
ax.set_ylabel('log Density')
ax.set_xlabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(temperatures,np.log10(densities),texgrid_145/texgrid_355M)#, vmax=1.3, vmin=0.8)
pl.colorbar(cax)
ax.set_title('$T_{ex}$ p-H$_2$CO 145 GHz/355 MHz')
ax.set_ylabel('log Density')
ax.set_xlabel('Temperature')

pl.figure(3)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(temperatures,np.log10(densities),taugrid_71M/taugrid_355M, vmax=45, vmin=-2)
pl.colorbar(cax)
ax.set_title('$\\tau$ p-H$_2$CO 71/355 MHz')
ax.set_ylabel('log Density')
ax.set_xlabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(temperatures,np.log10(densities),texgrid_71M/texgrid_355M, vmax=2, vmin=0.0)
pl.colorbar(cax)
ax.set_title('$T_{ex}$ p-H$_2$CO 71/355 MHz')
ax.set_ylabel('log Density')
ax.set_xlabel('Temperature')


pl.figure(4)
pl.clf()
extent = [densities.min(),densities.max(),temperatures.min(),temperatures.max()]
for kk,(freq,name) in enumerate(zip([0.071,145,0.355],['71M','145','355M'])):
    ax = pl.subplot(2,2,kk+1)
    grid = eval('fluxgrid_%s' % name)
    #ax.imshow(grid, extent=extent)
    pc = ax.pcolormesh(temperatures,np.log10(densities),grid)
    ax.set_title('Flux p-H$_2$CO %i GHz' % (freq))
    ax.set_ylabel('log Density')
    ax.set_xlabel('Temperature')
    pl.colorbar(pc)

pl.figure(6)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(temperatures,np.log10(densities),fluxgrid_71M/fluxgrid_355M)#, vmax=0.01, vmin=0)
pl.colorbar(cax)
ax.set_title('$S_{\\nu}$ p-H$_2$CO 71/355 MHz')
ax.set_ylabel('log Density')
ax.set_xlabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(temperatures,np.log10(densities),fluxgrid_145/fluxgrid_355M)#, vmax=1.2, vmin=0.8)
pl.colorbar(cax)
ax.set_title('$S_{\\nu}$ p-H$_2$CO 145 GHz/355 MHz')
ax.set_ylabel('log Density')
ax.set_xlabel('Temperature')


pl.show()

