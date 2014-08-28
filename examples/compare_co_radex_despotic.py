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

    R = pyradex.Radex(species='co', abundance=1e-4, column=None,
                      temperature=10,
                      collider_densities={'ph2':990,'oh2':10}, deltav=2.0,
                      datapath=datapath, escapeProbGeom='lvg')
    D = pyradex.despotic_interface.Despotic(hcolumn=2e21, species='co',
                                            temperature=10,
                                            abundance=1e-4,
                                            collider_densities={'ph2':990,'oh2':10},
                                            deltav=2.0,
                                            datapath=datapath,
                                            escapeProbGeom='lvg')

    # make sure we're dealing with identical qtys
    assert uvalue(D.temperature,u.K) == uvalue(R.temperature,u.K)
    assert R.density == D.density
    #np.testing.assert_approx_equal(uvalue(D.cloud.colDen/2,u.cm**-2), uvalue(R.h2column,u.cm**-2))

    # make sure RADEX converged
    assert R.run_radex() < 200

    print "RADEX N/dv/1.0645: ",R.radex.cphys.cdmol/(R.radex.cphys.deltav/1e5)/1.0645
    print "RADEX N/dv/1.27: ",R.radex.cphys.cdmol/(R.radex.cphys.deltav/1e5)/1.27
    print "DESPOTIC xs NH / dvdr",D.cloud.emitters['co'].abundance * D.cloud.colDen / (D.cloud.dVdr*3.08e18/1e5)

    TR = R.get_table()
    TD = D.get_table(noClump=True)

    R.source_area = 1*u.arcsec**2

    print("DV: D: {0}  R: {1}".format(D.deltav, R.deltav))
    print(R.line_flux_density[:5])
    print(R.line_brightness_temperature(1*u.arcsec**2)[:5])
    print(R.source_line_brightness_temperature[:5])
    print(R.T_B[:5])

    TR[:5]['frequency','Tex','tau','upperlevelpop','lowerlevelpop'].pprint()
    TD[:5]['frequency','Tex','tau','upperlevelpop','lowerlevelpop'].pprint()
    print(np.array(TR['tau'][:5]))
    print(np.array(TD['tau'][:5]))
    print(np.array(TR['tau'][:5])/np.array(TD['tau'][:5]))

    return D,R

def compare_levpops(D,R):
    uplevs = {'R':[], 'D': []}
    betas = {'R':[], 'D': []}

    for dv in range(1,6):
        R.deltav = dv
        # 1.27 gets them to very nearly agree, but still not quite, and not consistent.
        # I think I'm no closer to an answer.
        D.deltav = dv
        R.run_radex()
        D.recompute()
        uplevs['R'].append(R.upperlevelpop)
        uplevs['D'].append(D.upperlevelpop)
        betas['R'].append(R.beta[:len(D.beta)])
        betas['D'].append(D.beta)

    uplevs['R'] = np.array(uplevs['R'])
    uplevs['D'] = np.array(uplevs['D'])
    betas['R'] = np.array(betas['R'])
    betas['D'] = np.array(betas['D'])

    pl.plot((uplevs['R'][:,:6]-uplevs['D'][:,:6])/uplevs['D'][:,:6])
    #pl.plot(uplevs['D'][:,:4], linestyle=':')
    pl.figure(3)
    pl.clf()
    pl.plot((betas['R'][:,:5]-betas['D'][:,:5]).T, '-')

    return uplevs

if __name__ == '__main__':
    D,R = test_co()
    pl.figure(1)
    uplevs = compare_levpops(D,R)
    pl.figure(2)
    pl.plot(D.beta)
    pl.plot(R.beta[:len(D.beta)])
    pl.show()

