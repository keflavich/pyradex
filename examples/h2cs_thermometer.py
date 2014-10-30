"""
Simple proposal-writing experiment: For a given signal-to-noise in a line, what
signal-to-noise do you get in a derived parameter (e.g., temperature)?
"""
import pyradex
import pylab as pl
import numpy as np
import matplotlib
import astropy.units as u
import os

# for pretty on-screen plots
if os.path.exists('/Users/adam/.matplotlib/ggplotrc'):
    matplotlib.rc_file('/Users/adam/.matplotlib/ggplotrc')

# Download the data file if it's not here already
if not os.path.exists('ph2cs.dat'):
    import urllib
    import shutil
    fn,msg = urllib.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/ph2cs.dat')
    shutil.move(fn,'ph2cs.dat')

# Formatting tool
def latex_float(f):
    float_str = "{0:.1g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

# Define the grid parameters
ntemp = 25
temperatures = np.linspace(10,100,ntemp)

density = 1e4

abundance = 10**-8.5 # 10**-8.5
for abundance in (10**-8.5,10**-9):

        R = pyradex.Radex(species='ph2cs',
                          abundance=abundance,
                          collider_densities={'H2':density},
                          column=None,
                          temperature=temperatures[0])

        pl.figure(1)
        pl.clf()

        for temperature in [10,50,100]:
            R.temperature = temperature
            R.run_radex()
            wcs = np.linspace(103.0, 103.1, 1000)*u.GHz
            S = pyradex.synthspec.SyntheticSpectrum.from_RADEX(wcs, R,
                                                               linewidth=10*u.km/u.s)
            S.plot(label='%i K' % temperature)

        pl.legend(loc='best')

        # create a small grid...
        densities = [10**x for x in xrange(4,7)]
        ratio = {d:[] for d in densities}
        f1 = {d:[] for d in densities}
        f2 = {d:[] for d in densities}

        for density in densities:
            R.density = {'H2': density}
            for temperature in temperatures:
                R.temperature = temperature
                print R.run_radex(),

                F1 = R.T_B[2]
                F2 = R.T_B[12]

                ratio[density].append(F2/F1)
                f2[density].append(F2)
                f1[density].append(F1)
            print

        f1 = {d:np.array([x.value for x in f1[d]]) for d in densities}
        f2 = {d:np.array([x.value for x in f2[d]]) for d in densities}
        ratio = {d:np.array(ratio[d]) for d in densities}

        pl.figure(2)
        pl.clf()

        for d in densities:
            pl.plot(ratio[d],temperatures,label='$n=10^{%i}$' % (np.log10(d)))

            m = 1/((ratio[d][15]-ratio[d][5])/(temperatures[15]-temperatures[5]))
            b = temperatures[5]-ratio[d][5]*m
            line=(m,b)
            print d,m,b
            pl.plot(ratio[d],ratio[d]*line[0]+line[1],'--')

        pl.ylabel("Temperature")
        pl.xlabel("$S(3_{2,1}-2_{2,0})/S(3_{0,3}-2_{0,2})$")
        pl.legend(loc='best',fontsize=14)
        pl.title("X(p-H$_2$CS)$=10^{%0.1f}$" % (np.log10(abundance)))

        pl.axis([0,0.5,10,200,])

        pl.savefig("pH2CS_ratio_vs_temperature_X=%0.1e.pdf" % (abundance),bbox_inches='tight')

        pl.figure(3)
        pl.clf()
        for d in densities:
            L, = pl.plot(temperatures,f2[d],label='$n=10^{%i}$' % (np.log10(d)))
            pl.plot(temperatures,f1[d],'--',color=L.get_color())
        pl.xlabel("Temperature")
        pl.ylabel("$T_B(3_{2,1}-2_{2,0})$ (solid), $T_B(3_{0,3}-2_{0,2})$ (dashed)")
        ax = pl.gca()
        #pl.plot(ax.get_xlim(),[1.5,1.5],'k--',label='S/N$=5$, $\sigma=0.3$ K')
        #pl.plot(ax.get_xlim(),[2.5,2.5],'k:',label='S/N$=5$, $\sigma=0.5$ K')
        #pl.plot(ax.get_xlim(),[0.25,0.25],'k-.',label='S/N$=5$, $\sigma=0.05$ K')
        ax.axis([10,100,0,3.2])
        pl.legend(loc='best',fontsize=14)
        pl.title("X(p-H$_2$CS)$=10^{%0.1f}$" % (np.log10(abundance)))
        pl.savefig("pH2CS_321-220_vs_temperature_X=%0.1e.pdf" % (abundance),bbox_inches='tight')
