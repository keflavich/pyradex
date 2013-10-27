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

taugrid_6 = np.empty([ndens,ntemp])
texgrid_6 = np.empty([ndens,ntemp])
fluxgrid_6 = np.empty([ndens,ntemp])
taugrid_140 = np.empty([ndens,ntemp])
texgrid_140 = np.empty([ndens,ntemp])
fluxgrid_140 = np.empty([ndens,ntemp])
taugrid_150 = np.empty([ndens,ntemp])
texgrid_150 = np.empty([ndens,ntemp])
fluxgrid_150 = np.empty([ndens,ntemp])
columngrid = np.empty([ndens,ntemp])

R = pyradex.Radex(species='oh2co-h2', abundance=abundance)
R.run_radex()

# get the table so we can look at the frequency grid
table = R.get_table()

# Target frequencies:
table[np.array([0,1,3])].pprint()

for ii,tt in enumerate(temperatures):
    R.temperature = tt
    for jj,dd in enumerate(densities):
        R.density = {'oH2':dd*fortho,'pH2':dd*(1-fortho)}
        R.abundance = abundance # reset column to the appropriate value
        R.run_radex(reuse_last=False, reload_molfile=True)

        TI = R.total_intensity
        taugrid_6[jj,ii] = R.tau[0]
        texgrid_6[jj,ii] = R.tex[0].value
        fluxgrid_6[jj,ii] = TI[0].value
        taugrid_140[jj,ii] = R.tau[1]
        texgrid_140[jj,ii] = R.tex[1].value
        fluxgrid_140[jj,ii] = TI[1].value
        taugrid_150[jj,ii] = R.tau[3]
        texgrid_150[jj,ii] = R.tex[3].value
        fluxgrid_150[jj,ii] = TI[3].value
        columngrid[jj,ii] = R.column.value

pl.figure(1)
pl.clf()
extent = [densities.min(),densities.max(),temperatures.min(),temperatures.max()]
for kk,(grid,freq,label) in enumerate(zip([taugrid_140,taugrid_150,texgrid_140,texgrid_150],['140','150','140','150'],[r'\tau',r'\tau','T_{ex}','T_{ex}'])):
    ax = pl.subplot(2,2,kk+1)
    #ax.imshow(grid, extent=extent)
    ax.pcolormesh(np.log10(densities),temperatures,grid)
    ax.set_title('$%s$ o-H$_2$CO %s GHz' % (label,freq))
    ax.set_xlabel('log Density')
    ax.set_ylabel('Temperature')

pl.figure(2)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,taugrid_140/taugrid_150, vmax=1.3, vmin=0.8)
pl.colorbar(cax)
ax.set_title('$\\tau$ o-H$_2$CO 140/150 GHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,texgrid_140/texgrid_150, vmax=1.3, vmin=0.8)
pl.colorbar(cax)
ax.set_title('$T_{ex}$ o-H$_2$CO 140/150 GHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')

pl.figure(3)
pl.clf()
ax = pl.subplot(2,1,1)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,taugrid_6/taugrid_150, vmax=0.5, vmin=0.01)
pl.colorbar(cax)
ax.set_title('$\\tau$ o-H$_2$CO 6/150 GHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')

ax = pl.subplot(2,1,2)
#ax.imshow(grid, extent=extent)
cax = ax.pcolormesh(np.log10(densities),temperatures,texgrid_6/texgrid_150, vmax=2, vmin=0.1)
pl.colorbar(cax)
ax.set_title('$T_{ex}$ o-H$_2$CO 6/150 GHz')
ax.set_xlabel('log Density')
ax.set_ylabel('Temperature')


pl.figure(4)
pl.clf()
extent = [densities.min(),densities.max(),temperatures.min(),temperatures.max()]
for kk,freq in enumerate([6,140,150]):
    ax = pl.subplot(2,2,kk+1)
    grid = eval('fluxgrid_%i' % freq)
    #ax.imshow(grid, extent=extent)
    ax.pcolormesh(np.log10(densities),temperatures,grid)
    ax.set_title('Flux o-H$_2$CO %s GHz' % (label,freq))
    ax.set_xlabel('log Density')
    ax.set_ylabel('Temperature')

pl.show()
