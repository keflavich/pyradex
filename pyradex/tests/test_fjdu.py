from .. import fjdu
import numpy as np
from astropy import units as u

def test_simple():

    FF = fjdu.Fjdu(datapath='examples/', species='co', column=1e15,
                   density={'pH2':1e3,'oH2':0}, temperature=20)

    assert FF.params['ncol_x_cgs'] == 1e15
    np.testing.assert_almost_equal(FF.params['dens_x_cgs'], 1e3)
    #np.testing.assert_almost_equal(FF.params['oh2_density_cgs'], 1.7739422658722102)
    #np.testing.assert_almost_equal(FF.params['ph2_density_cgs'], 998.22605773412772)
    np.testing.assert_almost_equal(FF.params['ph2_density_cgs'], 1e3)
    assert FF.params['tkin'] == 20

    tbl = FF()

    #np.testing.assert_almost_equal(tbl[0]['Tex'], 8.69274406690759)
    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.6897105103500127)

def test_mod_params():

    FF = fjdu.Fjdu(datapath='examples/', species='co', column=1e15,
                   density={'pH2':1e3,'oH2':0}, temperature=20)

    tbl = FF()

    # Before para/ortho h2 setting, was this
    #np.testing.assert_almost_equal(tbl[0]['Tex'], 8.69274406690759)
    #np.testing.assert_almost_equal(tbl[0]['Tex'], 8.6935272872891236)
    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.6897105103500127)

    FF.column = 1e14
    tbl = FF()

    #np.testing.assert_almost_equal(tbl[0]['Tex'], 8.0986662583317646)
    #np.testing.assert_almost_equal(tbl[0]['Tex'], 8.0994405565362371)
    np.testing.assert_almost_equal(tbl[0]['Tex'], 8.0956672866767292)

    FF.density=1e4
    tbl = FF()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 25.382518594741391)

    FF.temperature=25
    tbl = FF()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 37.463006941695028)

    FF.deltav = 5 * u.km/u.s
    np.testing.assert_almost_equal(FF.deltav.to(u.km/u.s).value, 5)
    tbl = FF()
    np.testing.assert_almost_equal(tbl[0]['Tex'], 37.760227295047343)
