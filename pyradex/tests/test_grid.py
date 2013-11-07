from .. import Radex
import numpy as np

def buildgrid(densities=np.logspace(0,8,20), abundance=10**-8.5, fortho=1e-3, **kwargs):

    tau = np.zeros_like(densities)

    R = Radex(species='co', abundance=abundance, collider_densities={'oH2':10,'pH2':10})
    for jj,dd in enumerate(densities):
        R.density = {'oH2':dd*fortho,'pH2':dd*(1-fortho)}
        R.abundance = abundance # reset column to the appropriate value
        R.run_radex(**kwargs)
        tau[jj] = R.tau[0]

    return tau

def try_variants():

    taud = {}
    for reuse_last in (True,False):
        for reload_molfile in (True,False):
            s = ('T' if reuse_last else 'F') + ('T' if reload_molfile else 'F')
            taud[s] = buildgrid(reuse_last=reuse_last, reload_molfile=reload_molfile)

    return taud
