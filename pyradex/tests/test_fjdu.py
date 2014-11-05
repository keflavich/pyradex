import pyradex.fjdu
import numpy as np

def test_simple():

    FF = pyradex.fjdu.Fjdu(datapath='examples/', species='co', column=1e15,
                           density=1e3, temperature=20)

    assert FF.params['ncol_x_cgs'] == 1e15
    assert FF.params['dens_x_cgs'] == 1e3
    assert FF.params['h2_density_cgs'] == 1e3
    assert FF.params['tkin'] == 20

    tbl = FF()

    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.69274406690759)
