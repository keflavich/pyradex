import pyradex
import pytest
import os
import distutils.spawn
import numpy as np
from astropy import units as u

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
    data = pyradex.parse_outfile(data_path('example.out'))
    data.pprint(show_unit=True)

@pytest.mark.skipif(exepath=='radex',reason='radex not installed')
def test_call():
    data = pyradex.pyradex(executable=exepath,species='Radex/data/hco+',minfreq=50)
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

    data = pyradex.pyradex(executable=exepath,species=molecule,minfreq=1,maxfreq=250)
    data.pprint(show_unit=True)

def test_radex_class():
    R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column=1e15, collider_densities=None)
    assert hasattr(R,'radex')

def test_change_abundance():
    R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column=1e15, collider_densities=None)
    totdens = R.total_density
    R.abundance = 1e-6
    assert totdens == R.total_density

def test_consistent_abund():
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column=1e15,density=1e3)
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column=1e15,collider_densities={'H2':1e3})
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column_per_bin=1e15)
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=None,column=None)

def test_selfconsistent_density():
    rdx = pyradex.Radex(species='hco+', collider_densities={'H2':1e3}, column_per_bin=1e13)
    assert rdx.total_density.value == 1e3
    rdx.temperature = 30
    assert rdx.total_density.value == 1e3
    rdx.density = rdx.density
    assert rdx.total_density.value == 1e3
    rdx.density = {'H2':1e3}
    assert rdx.total_density.value == 1e3
    rdx.density = {'oH2':990,'pH2':10}
    assert rdx.total_density.value == 1e3

def test_consistent_parchanges():
    rdx = pyradex.Radex(species='hco+', collider_densities={'H2':1e3}, column_per_bin=1e13)
    np.testing.assert_almost_equal(rdx.abundance, 1e13/(1e3*(u.pc.to(u.cm))))
    assert rdx.locked_parameter == 'column'
    rdx.abundance=1e-9
    assert rdx.locked_parameter == 'abundance'
    np.testing.assert_almost_equal(rdx.total_density.to(u.cm**-3).value, 1e13/1e-9/u.pc.to(u.cm))
    rdx.density = 1e3
    rdx.column_per_bin = 1e13
    np.testing.assert_almost_equal(rdx.abundance, 1e13/(1e3*(u.pc.to(u.cm))))


if __name__ == "__main__":
    test_call()
    test_parse_example()
    for mol in ['co','13co','c18o','o-h2co','p-nh3']:
        test_molecules(mol)
