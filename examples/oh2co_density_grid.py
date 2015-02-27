import os
import pyradex
import pyradex.fjdu
import pylab as pl
import numpy as np
import matplotlib

# for pretty on-screen plots
if os.path.exists('/Users/adam/.matplotlib/ggplotrc'):
    matplotlib.rc_file('/Users/adam/.matplotlib/ggplotrc')

ndens = 100

#temperatures = np.linspace(10,50,nabund)
temperature = 50.
densities = np.logspace(1.5,7.,ndens)
abundances = np.array([-10, -9,-8],dtype='float')
opr = 1. # assume an even mix
fortho = opr/(1.+opr)
# background is... not obviously working
#background = 10.0 # instead of 2.73

nabund = abundances.size

taugrid_6cm = np.empty([ndens,nabund])
texgrid_6cm = np.empty([ndens,nabund])
tbgrid_6cm = np.empty([ndens,nabund])
fluxgrid_6cm = np.empty([ndens,nabund])

taugrid_2cm = np.empty([ndens,nabund])
texgrid_2cm = np.empty([ndens,nabund])
tbgrid_2cm = np.empty([ndens,nabund])
fluxgrid_2cm = np.empty([ndens,nabund])

taugrid_1cm = np.empty([ndens,nabund])
texgrid_1cm = np.empty([ndens,nabund])
tbgrid_1cm = np.empty([ndens,nabund])
fluxgrid_1cm = np.empty([ndens,nabund])

columngrid = np.empty([ndens,nabund])

if not os.path.exists('oh2co-h2.dat'):
    import urllib
    urllib.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/oh2co-h2.dat')

Radex = pyradex.fjdu.Fjdu
Radex = pyradex.Radex

R = Radex(species='oh2co-h2', collider_densities={'H2':1e6},
          column=1e12, temperature=temperature)
R.run_radex()

# get the table so we can look at the frequency grid
table = R.get_table()

# Target frequencies:
table[np.array([0,2,5])].pprint()

R.temperature = temperature
for ii,aa in enumerate(10.**abundances):
    for jj,dd in enumerate(densities):
        R.density = {'oH2':dd*fortho,'pH2':dd*(1-fortho)}
        R.abundance = aa
        R.run_radex(reuse_last=False, reload_molfile=True)

        TI = R.source_line_brightness_temperature
        taugrid_6cm[jj,ii] = R.tau[0]
        texgrid_6cm[jj,ii] = R.tex[0].value
        tbgrid_6cm[jj,ii] = R.T_B[0].value
        fluxgrid_6cm[jj,ii] = TI[0].value

        taugrid_2cm[jj,ii] = R.tau[2]
        texgrid_2cm[jj,ii] = R.tex[2].value
        tbgrid_2cm[jj,ii] = R.T_B[2].value
        fluxgrid_2cm[jj,ii] = TI[2].value

        taugrid_1cm[jj,ii] = R.tau[5]
        texgrid_1cm[jj,ii] = R.tex[5].value
        tbgrid_1cm[jj,ii] = R.T_B[5].value
        fluxgrid_1cm[jj,ii] = TI[5].value

        columngrid[jj,ii] = R.column.value


pl.figure(1)
pl.clf()
ax1 = pl.subplot(2,1,1)
ax1.loglog(densities,taugrid_6cm[:,0],label='$X=10^{%i}$' % abundances[0])
ax1.loglog(densities,taugrid_2cm[:,0])
ax1.loglog(densities,taugrid_6cm[:,1],linestyle='--',label='$X=10^{%i}$' % abundances[1])
ax1.loglog(densities,taugrid_2cm[:,1],linestyle='--')
ax1.loglog(densities,taugrid_6cm[:,2],linestyle=':',label='$X=10^{%i}$' % abundances[2])
ax1.loglog(densities,taugrid_2cm[:,2],linestyle=':')
ax1.set_xticks([])
ax1.set_ylabel("$\\tau$")
pl.legend(loc='best',fontsize=14)
ax2 = pl.subplot(2,1,2)
ax2.semilogx(densities,taugrid_6cm[:,0]/taugrid_2cm[:,0])
ax2.semilogx(densities,taugrid_6cm[:,1]/taugrid_2cm[:,1],linestyle='--')
ax2.set_xlabel("log n(H$_2$) [cm$^{-3}$]")
ax2.set_ylabel("Ratio")
pl.subplots_adjust(hspace=0)
pl.show()

pl.figure(2)
pl.clf()
ax1 = pl.subplot(2,1,1)
ax1.semilogx(densities,texgrid_6cm[:,0],label='$X=10^{%i}$' % abundances[0])
ax1.semilogx(densities,texgrid_2cm[:,0])
ax1.semilogx(densities,texgrid_6cm[:,1],linestyle='--',label='$X=10^{%i}$' % abundances[1])
ax1.semilogx(densities,texgrid_2cm[:,1],linestyle='--')
ax1.semilogx(densities,texgrid_6cm[:,2],linestyle=':',label='$X=10^{%i}$' % abundances[2])
ax1.semilogx(densities,texgrid_2cm[:,2],linestyle=':')
ax1.set_xticks([])
ax1.set_ylabel("$T_{ex}$")
ax1.hlines(2.73,10,1e7,color='k', linewidth=3, alpha=0.3, zorder=-1)
ax1.set_ylim(0,4)
ax1.set_xlim(densities.min(),densities.max())
pl.legend(loc='best', fontsize=14)
ax2 = pl.subplot(2,1,2)
ax2.semilogx(densities,texgrid_6cm[:,0]/texgrid_2cm[:,0])
ax2.semilogx(densities,texgrid_6cm[:,1]/texgrid_2cm[:,1],linestyle='--')
ax2.set_xlabel("log n(H$_2$) [cm$^{-3}$]")
ax2.set_ylabel("Ratio")
ax2.set_ylim(0,1)
pl.subplots_adjust(hspace=0)
pl.show()

pl.figure(3)
pl.clf()
ax1 = pl.subplot(2,1,1)
ax1.semilogx(densities,tbgrid_6cm[:,0],label='$X=10^{%i}$' % abundances[0])
ax1.semilogx(densities,tbgrid_2cm[:,0])
ax1.semilogx(densities,tbgrid_6cm[:,1],linestyle='--',label='$X=10^{%i}$' % abundances[1])
ax1.semilogx(densities,tbgrid_2cm[:,1],linestyle='--')
ax1.semilogx(densities,tbgrid_6cm[:,2],linestyle=':',label='$X=10^{%i}$' % abundances[2])
ax1.semilogx(densities,tbgrid_2cm[:,2],linestyle=':')
ax1.set_xticks([])
ax1.set_ylabel("$T_{B}$")
ax1.hlines(2.73,10,1e7,color='k', linewidth=3, alpha=0.3, zorder=-1)
ax1.set_ylim(-3,10)
ax1.set_xlim(densities.min(),densities.max())
pl.legend(loc='best', fontsize=14)
ax2 = pl.subplot(2,1,2)
ax2.semilogx(densities,tbgrid_6cm[:,0]/tbgrid_2cm[:,0])
ax2.semilogx(densities,tbgrid_6cm[:,1]/tbgrid_2cm[:,1],linestyle='--')
ax2.set_xlabel("log n(H$_2$) [cm$^{-3}$]")
ax2.set_ylabel("Ratio")
ax2.set_ylim(-1,10)
pl.subplots_adjust(hspace=0)
pl.show()


fig4 = pl.figure(4)
fig4.clf()
ax1 = pl.subplot(2,1,1)
cm2, = ax1.loglog(densities,taugrid_2cm[:,0],label='$X=10^{%i}$' % abundances[0])
cm1, = ax1.loglog(densities,taugrid_1cm[:,0])
ax1.loglog(densities,taugrid_2cm[:,1],linestyle='--',label='$X=10^{%i}$' % abundances[1],
           color=cm2.get_color())
ax1.loglog(densities,taugrid_1cm[:,1],linestyle='--',
           color=cm1.get_color())
ax1.loglog(densities,taugrid_2cm[:,2],linestyle=':',label='$X=10^{%i}$' % abundances[2],
           color=cm2.get_color())
ax1.loglog(densities,taugrid_1cm[:,2],linestyle=':',
           color=cm1.get_color())
ax1.set_xticks([])
ax1.set_ylabel("$\\tau$")
ax1.set_ylim(1e-5, 1e2)
pl.legend(loc='best',fontsize=14)
ax2 = pl.subplot(2,1,2)
ax2.semilogx(densities,1./(taugrid_2cm[:,0]/taugrid_1cm[:,0]))
ax2.semilogx(densities,1./(taugrid_2cm[:,1]/taugrid_1cm[:,1]),linestyle='--')
ax2.semilogx(densities,1./(taugrid_2cm[:,2]/taugrid_1cm[:,2]),linestyle=':')
ax2.set_xlabel("log n(H$_2$) [cm$^{-3}$]")
ax2.set_ylabel("Ratio")
pl.subplots_adjust(hspace=0)
pl.show()

pl.figure(6)
pl.clf()
ax1 = pl.subplot(2,1,1)
ax1.semilogx(densities,tbgrid_2cm[:,0],label='$X=10^{%i}$' % abundances[0])
ax1.semilogx(densities,tbgrid_1cm[:,0])
ax1.semilogx(densities,tbgrid_2cm[:,1],linestyle='--',label='$X=10^{%i}$' % abundances[1])
ax1.semilogx(densities,tbgrid_1cm[:,1],linestyle='--')
ax1.semilogx(densities,tbgrid_2cm[:,2],linestyle=':',label='$X=10^{%i}$' % abundances[2])
ax1.semilogx(densities,tbgrid_1cm[:,2],linestyle=':')
ax1.set_xticks([])
ax1.set_ylabel("$T_{B}$")
ax1.hlines(2.73,10,1e7,color='k', linewidth=3, alpha=0.3, zorder=-1)
ax1.set_ylim(-3,10)
ax1.set_xlim(densities.min(),densities.max())
pl.legend(loc='best', fontsize=14)
ax2 = pl.subplot(2,1,2)
ax2.semilogx(densities,tbgrid_2cm[:,0]/tbgrid_1cm[:,0])
ax2.semilogx(densities,tbgrid_2cm[:,1]/tbgrid_1cm[:,1],linestyle='--')
ax2.semilogx(densities,tbgrid_2cm[:,2]/tbgrid_1cm[:,2],linestyle=':')
ax2.set_xlabel("log n(H$_2$) [cm$^{-3}$]")
ax2.set_ylabel("Ratio")
ax2.set_ylim(-1,10)
pl.subplots_adjust(hspace=0)
pl.show()
