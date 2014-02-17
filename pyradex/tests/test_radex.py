import pyradex
import pytest
import os
import distutils.spawn

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
    R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column=1e15)
    assert hasattr(R,'radex')

def test_change_abundance():
    R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column=1e15)
    totdens = R.total_density
    R.abundance = 1e-6
    assert totdens == R.total_density

def test_consistent_abund():
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4)
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,column=1e15,h2column=1e20)
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=1e-4,h2column=None,column=None)
    with pytest.raises(ValueError):
        R = pyradex.Radex(datapath='examples/',species='co',abundance=None,h2column=None,column=None)

if __name__ == "__main__":
    test_call()
    test_parse_example()
    for mol in ['co','13co','c18o','o-h2co','p-nh3']:
        test_molecules(mol)
