import os
import pyradex
import pylab as pl
import numpy as np
import matplotlib

# for pretty on-screen plots
if os.path.exists('/Users/adam/.matplotlib/ggplotrc'):
    matplotlib.rc_file('/Users/adam/.matplotlib/ggplotrc')

ndens = 20

#temperatures = np.linspace(10,50,nabund)
temperature = 50
densities = np.logspace(1.5,7,ndens)
abundances = 10**np.array([-9,-8.5])
opr = 0.01 # assume primarily para
fortho = opr/(1+opr)

nabund = abundances.size

taugrid_6cm = np.empty([ndens,nabund])
texgrid_6cm = np.empty([ndens,nabund])
fluxgrid_6cm = np.empty([ndens,nabund])
taugrid_2cm = np.empty([ndens,nabund])
texgrid_2cm = np.empty([ndens,nabund])
fluxgrid_2cm = np.empty([ndens,nabund])
columngrid = np.empty([ndens,nabund])

if not os.path.exists('oh2co-h2.dat'):
    import urllib
    urllib.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/oh2co-h2.dat')

R = pyradex.Radex(species='oh2co-h2', abundance=abundances[0], h2column=1e21)
R.run_radex()

# get the table so we can look at the frequency grid
table = R.get_table()

# Target frequencies:
table[np.array([0,2])].pprint()

R.temperature = temperature
for ii,aa in enumerate(abundances):
    for jj,dd in enumerate(densities):
        R.density = {'oH2':dd*fortho,'pH2':dd*(1-fortho)}
        R.abundance = aa
        R.run_radex(reuse_last=False, reload_molfile=True)

        TI = R.source_line_brightness_temperature
        taugrid_6cm[jj,ii] = R.tau[0]
        texgrid_6cm[jj,ii] = R.tex[0].value
        fluxgrid_6cm[jj,ii] = TI[0].value
        taugrid_2cm[jj,ii] = R.tau[1]
        texgrid_2cm[jj,ii] = R.tex[1].value
        fluxgrid_2cm[jj,ii] = TI[1].value
        columngrid[jj,ii] = R.column.value


pl.figure(1)
pl.clf()
ax1 = pl.subplot(2,1,1)
ax1.loglog(densities,taugrid_6cm[:,0])
ax1.loglog(densities,taugrid_2cm[:,0])
ax1.loglog(densities,taugrid_6cm[:,1],linestyle='--')
ax1.loglog(densities,taugrid_2cm[:,1],linestyle='--')
ax1.set_xticks([])
ax1.set_ylabel("$\\tau$")
ax2 = pl.subplot(2,1,2)
ax2.semilogx(densities,taugrid_6cm[:,0]/taugrid_2cm[:,0])
ax2.semilogx(densities,taugrid_6cm[:,1]/taugrid_2cm[:,1],linestyle='--')
ax2.set_xlabel("log n(H$_2$) [cm$^{-3}$]")
ax2.set_ylabel("Ratio")
pl.subplots_adjust(hspace=0)
pl.show()
