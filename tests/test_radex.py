import pyradex
import pytest

def test_parse_example():
    data = pyradex.parse_outfile('data/example.out')
    data.pprint(show_units=True)

def test_call():
    data = pyradex.radex()
    data.pprint(show_units=True)

@pytest.mark.parametrize(('molecule',),zip(['co','13co','c18o','o-h2co','p-nh3']))
def test_molecules(molecule):
    data = pyradex.radex(molecule=molecule,flow=1,fhigh=250)
    data.pprint(show_units=True)

if __name__ == "__main__":
    test_call()
    test_parse_example()
    for mol in ['co','13co','c18o','o-h2co','p-nh3']:
        test_molecules(mol)
