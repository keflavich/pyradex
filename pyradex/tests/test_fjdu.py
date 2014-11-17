import pyradex.fjdu
import numpy as np
from astropy import units as u

def test_simple():

    FF = pyradex.fjdu.Fjdu(datapath='examples/', species='co', column=1e15,
                           density=1e3, temperature=20)

    assert FF.params['ncol_x_cgs'] == 1e15
    assert FF.params['dens_x_cgs'] == 1e3
    assert FF.params['h2_density_cgs'] == 1e3
    assert FF.params['tkin'] == 20

    tbl = FF()

    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.69274406690759)

def test_mod_params():

    FF = pyradex.fjdu.Fjdu(datapath='examples/', species='co', column=1e15,
                           density=1e3, temperature=20)

    tbl = FF()

    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.69274406690759)

    FF.column = 1e14
    tbl = FF()

    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.0986662583317646)

    FF.density=1e4
    tbl = FF()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 25.381267019506591)

    FF.temperature=25
    tbl = FF()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 37.455354325813083)

    FF.deltav = 5 * u.km/u.s
    np.testing.assert_almost_equal(FF.deltav.to(u.km/u.s).value, 5)
    tbl = FF()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 37.752430110617119)
