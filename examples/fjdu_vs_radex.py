import pyradex
import pyradex.fjdu
import pylab as pl
import numpy as np
from matplotlib.lines import Line2D


for density in ({'oH2': 100, 'pH2': 900, },
                {'oH2': 1000, 'pH2': 9000, },
                {'oH2': 10000, 'pH2': 90000, },):

    RR = pyradex.Radex(species='co', column=1e10, density=density, temperature=20)
    FF = pyradex.fjdu.Fjdu(species='co', column=1e10, density=density, temperature=20)

    fig1 = pl.figure(1)
    fig1.clf()
    fig2 = pl.figure(2)
    fig2.clf()
    fig1,(ax1a,ax1b,ax1c) = pl.subplots(nrows=3, num=1)
    fig2,(ax2a,ax2b,ax2c) = pl.subplots(nrows=3, num=2)

    for temperature in (20,500,1000):
        for column in (1e12, 1e14, 1e16):

            rtbl = RR(temperature=temperature, column=column, density=density)
            ftbl = FF(temperature=temperature, column=column, density=density)

            L1, = ax1a.plot(rtbl['upperlevel'], rtbl['upperlevelpop'], '-',
                            alpha=0.5, linewidth=2, label='Radex')
            L2, = ax1a.plot(ftbl['upperlevel'], ftbl['upperlevelpop'], '--',
                            alpha=0.5, linewidth=2, label="Fujun Du's MyRadex",
                            color=L1.get_color())
            ax1a.set_ylabel('Upper Energy Level Population')

            ax1b.plot(rtbl['upperlevel'], rtbl['T_B'], '-', alpha=0.5, linewidth=2,
                      label='Radex', color=L1.get_color())
            ax1b.semilogy(ftbl['upperlevel'], ftbl['T_B'], '--', alpha=0.5,
                          linewidth=2, label="Fujun Du's MyRadex",
                          color=L2.get_color())
            ax1b.set_ylabel('$T_B$')
            ax1b.set_ylim(1e-3,2e1)

            ax1c.set_ylabel('tau')
            ax1c.semilogy(rtbl['upperlevel'], rtbl['tau'], '-', alpha=0.5,
                          linewidth=2, label='Radex', color=L1.get_color())
            ax1c.plot(ftbl['upperlevel'], ftbl['tau'], '--', alpha=0.5,
                      linewidth=2, label="Fujun Du's MyRadex",
                      color=L2.get_color())
            ax1c.set_xlabel("Upper Energy Level Quantum Number")
            ax1c.set_ylim(1e-5,5e0)

            ax2a.plot(rtbl['upperlevel'],
                      rtbl['upperlevelpop']-ftbl['upperlevelpop'], '-', alpha=0.5,
                      linewidth=2, color=L1.get_color())
            ax2a.set_ylabel('Upper Energy Level Population')

            ax2b.plot(rtbl['upperlevel'], rtbl['T_B']-ftbl['T_B'], '-', alpha=0.5,
                      linewidth=2, label='Radex', color=L1.get_color())
            ax2b.set_ylabel('$T_B$')
            #ax2b.set_ylim(1e-3,2e1)

            ax2c.set_ylabel('tau')
            ax2c.plot(rtbl['upperlevel'], rtbl['tau']-ftbl['tau'], '-',
                          alpha=0.5, linewidth=2, color=L1.get_color())
            ax2c.set_xlabel("Upper Energy Level Quantum Number")
            #ax2c.set_ylim(1e-5,5e0)


    for ax in (ax1a,ax1b,ax1c,ax2a,ax2b,ax2c):
        ax.set_xlim(0,10)

    ax1a.set_title("Density = $10^{{{0}}}$ cm$^{{-3}}$".format(int(np.log10(FF.total_density.value))))
    ax2a.set_title("Density = $10^{{{0}}}$ cm$^{{-3}}$".format(int(np.log10(FF.total_density.value))))
    ax1a.legend((Line2D([0],[0], linewidth=2, color='k', linestyle='-'),
                Line2D([0],[0], linewidth=2, color='k', linestyle='--')),
               ('Radex', "Fujun Du's MyRadex"), loc='best')
    fig1.savefig("fjdu_vs_radex_CO_n{0}.png".format(int(np.log10(FF.total_density.value))))
    pl.draw()
    pl.show()
