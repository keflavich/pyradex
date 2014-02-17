import pyradex
import pylab as pl
from astropy import units as u
from pyradex.utils import united,uvalue,get_datafile
import numpy as np
import os

def test_co():
    mydir = os.path.dirname(os.path.realpath(__file__))

    print mydir
    datapath,speciesdat = get_datafile('co', savedir=os.path.join(mydir,'data/'))
    print datapath,speciesdat

    R = pyradex.Radex(h2column=1e21,species='co',abundance=1e-5,column=None,datapath=datapath)
    D = pyradex.despotic_interface.Despotic(hcolumn=2e21,species='co',abundance=1e-5/2,datapath=datapath)

    # make sure we're dealing with identical qtys
    assert uvalue(D.temperature,u.K) == uvalue(R.temperature,u.K)
    assert R.density == D.density
    np.testing.assert_approx_equal(uvalue(D.cloud.colDen/2,u.cm**-2), uvalue(R.h2column,u.cm**-2))

    # make sure RADEX converged
    assert R.run_radex() < 200

    print "RADEX N/dv/1.064: ",R.radex.cphys.cdmol/(R.radex.cphys.deltav/1e5)/1.064
    print "DESPOTIC xs NH / dvdr",D.cloud.emitters['co'].abundance * D.cloud.colDen / (D.cloud.dVdr*3.08e18/1e5)

    TR = R.get_table()
    TD = D.get_table(noClump=True)

    R.source_area = 1*u.arcsec**2

    print(R.line_flux_density)
    print(R.line_brightness_temperature(1*u.arcsec**2))
    print(R.source_line_brightness_temperature)
    print(R.T_B)

    print(TR[:5])
    print(TD[:5])
    print(np.array(TR['tau'][:5]))
    print(np.array(TD['tau'][:5]))
    print(np.array(TR['tau'][:5])/np.array(TD['tau'][:5]))

