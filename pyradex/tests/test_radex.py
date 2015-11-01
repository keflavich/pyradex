import pytest
import os
import distutils.spawn
import numpy as np
from astropy import units as u
from astropy import log

from ..core import parse_outfile,pyradex,Radex

exepath = 'Radex/bin/radex'
#if os.path.isfile(exepath) and os.access(exepath, os.X_OK):
#    exepath = exepath
if distutils.spawn.find_executable('radex'):
    exepath = distutils.spawn.find_executable('radex')
else:
    exepath = 'radex'

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

def test_parse_example():
    data = parse_outfile(data_path('example.out'))
    data.pprint(show_unit=True)

@pytest.mark.skipif(exepath=='radex',reason='radex not installed')
def test_call():
    data = pyradex(executable=exepath,species='Radex/data/hco+',minfreq=50)
    data.pprint(show_unit=True)

@pytest.mark.skipif(exepath=='radex',reason='radex not installed')
@pytest.mark.parametrize('molecule,',('co','13co','c18o','o-h2co','p-nh3',))
def test_molecules(molecule):
    if os.path.isfile('examples/%s.dat' % molecule):
        molecule = 'examples/%s' % molecule
    elif os.path.isfile('Radex/data/%s.dat' % molecule):
        molecule = 'Radex/data/%s' % molecule
    else:
        return

    data = pyradex(executable=exepath,species=molecule,minfreq=1,maxfreq=250)
    data.pprint(show_unit=True)

def test_radex_class():
    R = Radex(datapath='examples/',species='co',abundance=1e-4,
                      column=1e15, collider_densities=None, temperature=20)
    assert hasattr(R,'radex')

def test_change_abundance():
    R = Radex(datapath='examples/',species='co',abundance=1e-4,
                      column=1e15, collider_densities=None, temperature=20)
    totdens = R.total_density
    R.abundance = 1e-6
    assert totdens == R.total_density

def test_consistent_abund():
    with pytest.raises(ValueError):
        R = Radex(datapath='examples/', species='co', abundance=1e-4,
                          column=1e15, density=1e3)
    with pytest.raises(ValueError):
        R = Radex(datapath='examples/', species='co', abundance=1e-4,
                          column=1e15, collider_densities={'H2':1e3})
    with pytest.raises(ValueError):
        R = Radex(datapath='examples/', species='co', abundance=1e-4,
                          column_per_bin=1e15)
    with pytest.raises(ValueError):
        R = Radex(datapath='examples/', species='co', abundance=None,
                          column=None)

def test_selfconsistent_density():
    rdx = Radex(species='co', collider_densities={'H2':1e3},
                        column_per_bin=1e13, temperature=20)
    np.testing.assert_almost_equal(rdx.total_density.value, 1e3)
    rdx.temperature = 30
    np.testing.assert_almost_equal(rdx.total_density.value, 1e3)
    rdx.density = rdx.density
    np.testing.assert_almost_equal(rdx.total_density.value, 1e3)
    rdx.density = {'H2':1e3}
    np.testing.assert_almost_equal(rdx.total_density.value, 1e3)
    rdx.density = {'oH2':990,'pH2':10}
    np.testing.assert_almost_equal(rdx.total_density.value, 1e3)

def test_consistent_parchanges():
    rdx = Radex(species='co', collider_densities={'H2':1e3},
                        column_per_bin=1e13, temperature=20)
    np.testing.assert_almost_equal(rdx.abundance, 1e13/(1e3*(u.pc.to(u.cm))))
    assert rdx.locked_parameter == 'column'
    rdx.abundance=1e-9
    assert rdx.locked_parameter == 'abundance'
    np.testing.assert_almost_equal(rdx.total_density.to(u.cm**-3).value, 1e13/1e-9/u.pc.to(u.cm))
    rdx.density = 1e3
    rdx.column_per_bin = 1e13
    np.testing.assert_almost_equal(rdx.abundance, 1e13/(1e3*(u.pc.to(u.cm))))

def test_radex_results():
    # default parameters for radex online
    rdx = Radex(species='co', collider_densities={'H2':1e4}, column_per_bin=1e14, deltav=1.0,
                        temperature=30, tbackground=2.73)
    rdx.run_radex()
    assert rdx.temperature.value == 30.0 # no approximates allowed
    assert rdx.column.value == 1e14
    #       LINE         E_UP       FREQ        WAVEL     T_EX      TAU        T_R       POP        POP       FLUX        FLUX
    #                    (K)        (GHz)       (um)      (K)                  (K)        UP        LOW      (K*km/s) (erg/cm2/s)
    # 1      -- 0          5.5    115.2712   2600.7576   56.131  1.786E-03  9.378E-02  3.640E-01  1.339E-01  9.983E-02  1.969E-09
    # RADEX online:
    # Transition        Frequency       Tex     tau             TR
    # 1   --   0        115.2712        54.863  1.824E-03       9.349E-02
    np.testing.assert_approx_equal(rdx.tex[0].value, 56.131, 5)
    np.testing.assert_approx_equal(rdx.tau[0], 1.786E-03, 4)
    np.testing.assert_approx_equal(rdx.upperlevelpop[0], 3.640E-01, 4)
    np.testing.assert_approx_equal(rdx.lowerlevelpop[0], 1.339E-01, 4)

def test_consistent_init():
    rdx = Radex(species='co', collider_densities={'H2':1e4}, column_per_bin=1e14, deltav=1.0,
                        temperature=30, tbackground=2.73)
    log.debug('crate at init: {0}'.format(rdx.radex.collie.crate[0,1]))
    log.debug(str((rdx.density, rdx.temperature, rdx.column, rdx._cddv)))
    rdx.run_radex()
    log.debug('crate at run: {0}'.format(rdx.radex.collie.crate[0,1]))
    tex0 = rdx.tex[0].value
    crate0 = np.copy(rdx.radex.collie.crate)
    rdx = Radex(species='co', collider_densities={'H2':1e4}, column_per_bin=1e14, deltav=1.0,
                        temperature=25, tbackground=2.73)
    log.debug('crate at temchange: {0}'.format(rdx.radex.collie.crate[0,1]))
    rdx = Radex(species='co', collider_densities={'H2':1e4}, column_per_bin=1e14, deltav=1.0,
                        temperature=30, tbackground=2.73)
    log.debug('crate at init2: {0}'.format(rdx.radex.collie.crate[0,1]))
    log.debug(str((rdx.density, rdx.temperature, rdx.column, rdx._cddv)))
    rdx.run_radex()
    log.debug('crate at run2: {0}'.format(rdx.radex.collie.crate[0,1]))
    tex1 = rdx.tex[0].value
    crate1 = np.copy(rdx.radex.collie.crate)
    assert tex0==tex1
    assert np.all(crate0 == crate1)

def test_thermal_opr():
    # Check if H2 is specified as total H2, the thermal fraction of O/P H2 is used
    rdx = Radex(species='co', collider_densities={'H2':1e4}, column_per_bin=1e14, deltav=1.0,
                        temperature=30, tbackground=2.73)
    opr = 9.0*np.exp(-170.6/30)
    fortho = opr/(1+opr)
    np.testing.assert_almost_equal(rdx.density['oH2'].value, fortho*1e4)
    np.testing.assert_almost_equal(rdx.density['pH2'].value, (1-fortho)*1e4)

    rdx.temperature = 50
    opr = 9.0*np.exp(-170.6/50)
    fortho = opr/(1+opr)
    np.testing.assert_almost_equal(rdx.density['oH2'].value, fortho*1e4)
    np.testing.assert_almost_equal(rdx.density['pH2'].value, (1-fortho)*1e4)

    # Check that if ortho is specified, density remains unchanged
    rdx = Radex(species='co', collider_densities={'oH2':1e4, 'pH2':0}, column_per_bin=1e14, deltav=1.0,
                        temperature=30, tbackground=2.73)
    assert rdx.density['oH2'].value == 1e4
    rdx.temperature = 50
    assert rdx.density['oH2'].value == 1e4

def test_mod_params():

    RR = Radex(datapath='examples/', species='co', column=1e15,
                       density=1e3, temperature=20)

    tbl = RR()

    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.69274406690759, decimal=2)

    RR.column = 1e14
    tbl = RR()

    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.0986662583317646, decimal=2)

    RR.density=1e4
    tbl = RR()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 25.381267019506591, decimal=1)

    RR.temperature=25
    tbl = RR()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 37.88, decimal=1)

    RR.deltav = 5 * u.km/u.s
    np.testing.assert_almost_equal(RR.deltav.to(u.km/u.s).value, 5)
    tbl = RR()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 37.83, decimal=1)

if __name__ == "__main__":
    test_call()
    test_parse_example()
    for mol in ['co','13co','c18o','o-h2co','p-nh3']:
        test_molecules(mol)
