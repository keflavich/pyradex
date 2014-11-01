import pyradex
import pyradex.fjdu
import pylab as pl

density = {'oH2': 100,
           'pH2': 900,
          }

RR = pyradex.Radex(species='co', column=1e10, density=density, temperature=20)
FF = pyradex.fjdu.Fjdu(species='co', column=1e10, density=density, temperature=20)

fig = pl.figure(1)
fig.clf()
ax1 = pl.subplot(3,1,1)
ax2 = pl.subplot(3,1,2)
ax3 = pl.subplot(3,1,3)

for temperature in (20,100,300):
    for column in (1e12, 1e14, 1e16):

        rtbl = RR(temperature=temperature, column=column, density=density)
        ftbl = FF(temperature=temperature, column=column, density=density)

        ax1.plot(rtbl['upperlevel'], rtbl['upperlevelpop'], '-', alpha=0.5,
                label='Radex')
        ax1.plot(ftbl['upperlevel'], ftbl['upperlevelpop'], '--', alpha=0.5,
                label="Fujun Du's MyRadex")
        ax1.set_ylabel('Upper Energy Level Population')

        ax2.plot(rtbl['upperlevel'], rtbl['T_B'], '-', alpha=0.5,
                label='Radex')
        ax2.plot(ftbl['upperlevel'], ftbl['T_B'], '--', alpha=0.5,
                label="Fujun Du's MyRadex")
        ax2.set_ylabel('$T_B$')

        ax3.set_ylabel('tau')
        ax3.plot(rtbl['upperlevel'], rtbl['tau'], '-', alpha=0.5,
                label='Radex')
        ax3.plot(ftbl['upperlevel'], ftbl['tau'], '--', alpha=0.5,
                label="Fujun Du's MyRadex")
        ax3.set_xlabel("Upper Energy Level Quantum Number")

for ax in (ax1,ax2,ax3):
    ax.set_xlim(0,10)

#ax3.legend(loc='best')
pl.draw()
pl.show()
