"""
Create some simple grids for the 1mm para-H2CO lines
"""
import pyradex
import numpy as np

ntemp,ndens = 50,20

temperatures = np.linspace(10,500,ntemp)
densities = np.logspace(2.5,7,ndens)
abundance = 10**-8.5
opr = 0.01 # assume primarily para
fortho = opr/(1+opr)

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

import os
if not os.path.exists('ph2co-h2.dat'):
    import urllib
    urllib.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/ph2co-h2.dat')

R = pyradex.Radex(species='ph2co-h2', abundance=abundance, h2column=1e23,
                  temperature=50,
                  collider_densities={'oH2':2e4*fortho,'pH2':2e4*(1-fortho)})
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

if False:

    for ii,tt in enumerate(temperatures):
        R.temperature = tt
        for jj,dd in enumerate(densities):
            R.density = {'oH2':dd*fortho,'pH2':dd*(1-fortho)}
            R.abundance = abundance # reset column to the appropriate value
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
