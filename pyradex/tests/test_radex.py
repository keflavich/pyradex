import pyradex
import pytest
import os

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

def test_parse_example():
    data = pyradex.parse_outfile(data_path('example.out'))
    data.pprint(show_unit=True)

def test_call():
    data = pyradex.pyradex()
    data.pprint(show_unit=True)

@pytest.mark.parametrize(('molecule',),('co','13co','c18o','o-h2co','p-nh3',))
def test_molecules(molecule):
    data = pyradex.pyradex(species=molecule,minfreq=1,maxfreq=250)
    data.pprint(show_unit=True)

if __name__ == "__main__":
    test_call()
    test_parse_example()
    for mol in ['co','13co','c18o','o-h2co','p-nh3']:
        test_molecules(mol)
