import pyradex

def test_parse_example():
    data = pyradex.parse_outfile('example.out')
    data.pprint(show_units=True)

def test_call():
    data = pyradex.radex()
    data.pprint(show_units=True)

if __name__ == "__main__":
    test_call()
    test_parse_example()
