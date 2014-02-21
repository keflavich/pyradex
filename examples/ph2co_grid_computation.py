"""
Create some simple grids for the low-frequency para-H2CO lines
"""
import pyradex
import numpy as np

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

R = pyradex.Radex(species='ph2co-h2', abundance=abundance, h2column=1e23)
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

        TI = R.source_line_surfbrightness
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
