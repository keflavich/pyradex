"""
Create some simple grids for the 1mm para-H2CO lines
"""
import pyradex
import numpy as np

ntemp,ndens = 50,20

temperatures = np.linspace(10,500,ntemp)
densities = np.logspace(2.5,7,ndens)
abundance = 10**-8.5
abundance = 1.2e-9 # Johnston / Ao
opr = 0.01 # assume primarily para
opr = 3
fortho = opr/(1+opr)

import os
if not os.path.exists('ph2co-h2.dat'):
    import urllib
    urllib.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/ph2co-h2.dat')

# 1e23 from Johnston 2014
R = pyradex.Radex(species='ph2co-h2', abundance=abundance, h2column=1e23,
                  temperature=50,
                  collider_densities={'oH2':2e4*fortho,'pH2':2e4*(1-fortho)})
print R.escapeProbGeom # DEBUG
R.run_radex()

# get the table so we can look at the frequency grid
table = R.get_table()

# Target frequencies:
table[np.array([6,1,11])].pprint()

key_303 = np.where((table['upperlevel'] == '3_0_3') &
                   (table['frequency'] > 218) &
                   (table['frequency'] < 220))[0]
key_321 = np.where((table['upperlevel'] == '3_2_1') &
                   (table['frequency'] > 218) &
                   (table['frequency'] < 220))[0]
key_322 = np.where((table['upperlevel'] == '3_2_2') &
                   (table['frequency'] > 218) &
                   (table['frequency'] < 220))[0]


def compute_grid(densities=densities, temperatures=temperatures, fortho=fortho,
                 columnperbin=abundance*1e23, deltav=1.0, escapeProbGeom='lvg',
                 R=R):

    R.escapeProbGeom=escapeProbGeom

    ndens = len(densities)
    ntemp = len(temperatures)

    taugrid_303 = np.empty([ndens,ntemp])
    texgrid_303 = np.empty([ndens,ntemp])
    fluxgrid_303 = np.empty([ndens,ntemp])
    taugrid_321 = np.empty([ndens,ntemp])
    texgrid_321 = np.empty([ndens,ntemp])
    fluxgrid_321 = np.empty([ndens,ntemp])
    taugrid_322 = np.empty([ndens,ntemp])
    texgrid_322 = np.empty([ndens,ntemp])
    fluxgrid_322 = np.empty([ndens,ntemp])
    columngrid = np.empty([ndens,ntemp])

    for ii,tt in enumerate(temperatures):
        R.temperature = tt
        for jj,dd in enumerate(densities):
            R.density = {'oH2':dd*fortho,'pH2':dd*(1-fortho)}
            #R.abundance = abundance # reset column to the appropriate value
            R.column_per_bin = columnperbin
            R.deltav = deltav
            R.run_radex(reuse_last=False, reload_molfile=True)

            TI = R.source_line_surfbrightness
            taugrid_303[jj,ii] = R.tau[key_303]
            texgrid_303[jj,ii] = R.tex[key_303].value
            fluxgrid_303[jj,ii] = TI[key_303].value
            taugrid_321[jj,ii] = R.tau[key_321]
            texgrid_321[jj,ii] = R.tex[key_321].value
            fluxgrid_321[jj,ii] = TI[key_321].value
            taugrid_322[jj,ii] = R.tau[key_322]
            texgrid_322[jj,ii] = R.tex[key_322].value
            fluxgrid_322[jj,ii] = TI[key_322].value
            columngrid[jj,ii] = R.column.value

    return (TI,
            taugrid_303,texgrid_303,fluxgrid_303,
            taugrid_321,texgrid_321,fluxgrid_321,
            taugrid_322,texgrid_322,fluxgrid_322,
            columngrid)
(TI, taugrid_303,texgrid_303,fluxgrid_303,
 taugrid_321,texgrid_321,fluxgrid_321,
 taugrid_322,texgrid_322,fluxgrid_322,
 columngrid) = compute_grid()

if __name__ == "__main__":
    import pylab as pl
    from matplotlib.ticker import FuncFormatter
    import matplotlib


    johnstonpars = [
                    (1e13, 10, ),
                    (1e14, 5, ),
                    (1e14, 10, ),
                    (1e14, 20, ),
                    (1e15, 10, ),
        ]

    for geom in ('lvg','sphere'):
        for opr in (0.01, 3):
            for colh2co, dv in johnstonpars:
                (TI, taugrid_303,texgrid_303,fluxgrid_303,
                 taugrid_321,texgrid_321,fluxgrid_321,
                 taugrid_322,texgrid_322,fluxgrid_322,
                 columngrid) = compute_grid(columnperbin=colh2co, deltav=dv, escapeProbGeom=geom)

                fig = pl.figure(1)
                fig.clf()
                ax = fig.gca()

                im = ax.imshow((fluxgrid_303/fluxgrid_321).T, vmin=1.0, vmax=10, aspect=fluxgrid_303.shape[0]/float(fluxgrid_303.shape[1]),
                               norm=matplotlib.colors.LogNorm())
                fig.colorbar(im)
                c = ax.contour((fluxgrid_303/fluxgrid_321).T, levels=[1.4,1.5,1.6,1.7,1.8], colors=['w']*3)
                ax.clabel(c, inline=1, fontsize=10, color='w') 
                ax.xaxis.set_major_formatter(FuncFormatter(lambda x,y: "{0:0.2f}".format(np.log10(densities[int(x)])) if x<len(densities) else ""))
                ax.yaxis.set_major_formatter(FuncFormatter(lambda x,y: str(temperatures[int(x)]) if x<len(temperatures) else ""))
                ax.set_ylabel("Temperature (K)")
                ax.set_xlabel("Log Density (cm$^{-3}$)")
                ax.set_title("N(p-H$_2$CO) = $10^{{{0}}}$, geom={1}, dv={2}, opr={3}".format(np.log10(R.column_per_bin.value),
                                                                                             R.escapeProbGeom,
                                                                                             R.deltav.value,
                                                                                             opr))
                pl.savefig("paraH2CO_thermometry_303to321_geom{1}_col{0}_deltav{2}_opr{3}_denstemgrid.png".format(np.log10(R.column_per_bin.value),
                                                                                                         R.escapeProbGeom,
                                                                                                         R.deltav.value,
                                                                                                         opr),
                           bbox_inches='tight')

                fig = pl.figure(2)
                fig.clf()
                ax = fig.gca()

                im = ax.imshow((fluxgrid_303/fluxgrid_322).T, vmin=1.0, vmax=10, aspect=fluxgrid_303.shape[0]/float(fluxgrid_303.shape[1]),
                               norm=matplotlib.colors.LogNorm())
                fig.colorbar(im)
                c = ax.contour((fluxgrid_303/fluxgrid_322).T, levels=[1.4,1.5,1.6,1.7,1.8], colors=['w']*3)
                ax.clabel(c, inline=1, fontsize=10, color='w') 
                ax.xaxis.set_major_formatter(FuncFormatter(lambda x,y: "{0:0.2f}".format(np.log10(densities[int(x)])) if x<len(densities) else ""))
                ax.yaxis.set_major_formatter(FuncFormatter(lambda x,y: str(temperatures[int(x)]) if x<len(temperatures) else ""))
                ax.set_ylabel("Temperature (K)")
                ax.set_xlabel("Log Density (cm$^{-3}$)")
                ax.set_title("N(p-H$_2$CO) = $10^{{{0}}}$, geom={1}, dv={2}, opr={3}".format(np.log10(R.column_per_bin.value),
                                                                                             R.escapeProbGeom,
                                                                                             R.deltav.value,
                                                                                             opr))
                fig.savefig("paraH2CO_thermometry_303to322_geom{1}_col{0}_deltav{2}_opr{3}_denstemgrid.png".format(np.log10(R.column_per_bin.value),
                                                                                                         R.escapeProbGeom,
                                                                                                         R.deltav.value,
                                                                                                         opr),
                           bbox_inches='tight')
                pl.show()
