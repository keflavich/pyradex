from __future__ import print_function
import numpy as np
import os
from astropy import constants
from astropy import log
import astropy.table
from .. import base_class
from ..utils import ImmutableDict,unitless,lower_keys

import astropy.units as u

class Fjdu(base_class.RadiativeTransferApproximator):
    def __init__(self, **kwargs):
        self.set_default_params()
        from myradex import myradex_wrapper
        self._myradex = myradex_wrapper

    def __call__(self, return_table=True, **kwargs):

        niter = self.run_radex(**kwargs)

        if return_table:
            return self.get_table()
        else:
            return niter

    def load_datfile(self, filename=None, verbose=False):
        filename = filename or self.molpath
        if filename is not None:
            self.dpath = os.path.dirname(filename) or self.dpath
            self.fname = os.path.basename(filename)
        nlevels, nitems, ntrans = self._myradex.config_basic(self.dpath,
                                                             self.fname,
                                                             verbose)
        self.params['n_levels'] = nlevels
        self.params['n_items'] = nitems
        self.params['n_transitions'] = ntrans

    _default_params = (('tkin', 1e3),
                       ('dv_CGS', 1e5),
                       ('dens_X_CGS', 1e6),
                       ('Ncol_X_CGS', 1e15),
                       ('H2_density_CGS', 1e9),
                       ('HI_density_CGS', 1e1),
                       ('oH2_density_CGS', 0.0),
                       ('pH2_densty_CGS', 0.0),
                       ('HII_density_CGS', 0.0),
                       ('Electron_density_CGS', 1e6),
                       ('n_levels', 0),
                       ('n_item', 0),
                       ('n_transitions', 0))

    def set_default_params(self):
        self.params = dict(self._default_params)

    def set_params(self, **kwargs):
        self.params.update(kwargs)

    @property
    def params(self):
        return lower_keys(self._params)

    @params.setter
    def params(self, value):
        if not isinstance(value, dict):
            raise TypeError('Parameters must be a dictionary.')
        default = lower_keys(dict(self._default_params))
        self._params = default
        for k in value:
            if k.lower() not in default:
                raise ValueError("{0} is not a valid key.".format(k))
            else:
                self._params[k] = value[k]

    @property
    def density(self):

        dd = {'H2': u.Quantity(self.params['h2_density_cgs'], self._u_cc),
              'OH2': u.Quantity(self.params['oh2_density_cgs'], self._u_cc),
              'PH2': u.Quantity(self.params['ph2_density_cgs'], self._u_cc),
              'E': u.Quantity(self.params['Electron_density_cgs'], self._u_cc),
              'H+': u.Quantity(self.params['HII_density_cgs'], self._u_cc),
              'H': u.Quantity(self.params['HI_density_cgs'], self._u_cc),
              'He': u.Quantity(0, self._u_cc),}
        return ImmutableDict(dd)

    @density.setter
    def density(self, value):
        raise NotImplementedError

    @property
    def temperature(self):
        return u.Quantity(self.params[tkin], u.K)

    @temperature.setter
    def temperature(self, tkin):
        if hasattr(tkin,'to'):
            tkin = unitless(u.Quantity(tkin, u.K))
        if tkin <= 0 or tkin > 1e4:
            raise ValueError('Must have kinetic temperature > 0 and < 10^4 K')
        self.params['tkin'] = tkin

    @property
    def column_per_bin(self):
        return u.Quantity(self.params['ncol_x_cgs'], self._u_sc)

    @column_per_bin.setter
    def column_per_bin(self, col):
        if hasattr(col,'to'):
            col = unitless(u.Quantity(col, self._u_sc))
        if col < 1e5 or col > 1e25:
            raise ValueError("Extremely low or extremely high column.")
        self.params['ncol_x_cgs'] = col

        col = u.Quantity(col, self._u_sc)
        if not self._is_locked:
            self._is_locked = True
            if self.locked_parameter == 'density':
                ab = (col/(self.total_density * self.length))
                if hasattr(ab, 'decompose'):
                    self.abundance = ab.decompose().value
                else:
                    self.abundance = ab / (self._u_cc*u.pc).to(self._u_sc)
            elif self.locked_parameter == 'abundance':
                self.density = col / self.length / self.abundance
            self._lock_param('column')
            self._is_locked = False

    @property
    def abundance(self):
        return self._abundance

    @abundance.setter
    def abundance(self, abund):
        self._abundance = abund
        if not self._is_locked:
            self._is_locked = True
            if self.locked_parameter == 'column':
                dens = self.column_per_bin / self.length / abund
                self.density = dens
            elif self.locked_parameter == 'density':
                col = self.total_density*self.length*abund
                self.column_per_bin = u.Quantity(col, u.cm**-2)
            self._lock_param('abundance')
            self._is_locked=False

    @property
    def deltav(self):
        return u.Quantity(self.params['dv_cgs']/self._kms_to_cms, self._u_kms)

    _kms_to_cms = 1e-5
    _u_cms = u.cm/u.s

    @deltav.setter
    def deltav(self, dv):
        self.params['dv_cgs'] = unitless(u.Quantity(dv*self._kms_to_cms,
                                                    self._u_cms))

    @property
    def molpath(self):
        if hasattr(self,'_molpath'):
            return self._molpath

    @molpath.setter
    def molpath(self, molfile):
        if "~" in molfile:
            molfile = os.path.expandpath(molfile)
        utils.verify_collisionratefile(molfile)
        self._molpath = molfile

    @property
    def escprobProbGeom(self):
        return 'lvg'
        
    def run_radex(self, **kwargs):
        self.set_params(**kwargs)
        energies, f_occupations, data_transitions, cooling_rate = \
                self._myradex.run_one_params(**self.params)
        self._energies = u.Quantity(energies, u.K) # excitation temperature
        self._data_dict = cast_into_dic("".join(self._myradex.column_names),
                                        data_transitions)
        self._level_population = f_occupations

    _um_to_ghz = u.um.to(u.GHz, equivalencies=u.spectral())

    @property
    def frequency(self):
        return u.Quantity(self._data_dict[lam]*_um_to_ghz, unit=u.GHz)

    @property
    def level_population(self):
        return self._level_population

    @property
    def tex(self):
        return u.Quantity(self._data_dict['Tex'], u.K)

    Tex = tex

    @property
    def tau(self):
        return self._data_dict['tau']

    @property
    def upperstateenergy(self):
        return u.Quantity(self._data_dict['Eup'], u.K)

    @property
    def upperlevel(self):
        return self._data_dict['iup']

    @property
    def lowerlevel(self):
        return self._data_dict['ilow']

    @property
    def upperlevelpop(self):
        return self._data_dict['fup']

    @property
    def lowerlevelpop(self):
        return self._data_dict['flow']

    @property
    def source_line_brightness_temperature(self):
        return u.Quantity(self._data_dict['flux_K'], u.K)

    @property
    def source_line_brightness(self):
        return u.Quantity(self._data_dict['flux'], self._u_brightness)

    @property
    def beta(self):
        return self._data_dict['beta']


def cast_into_dic(col_names, arr):
    '''col_names is column_info, and arr is data_transitions'''
    names = col_names.split()
    return {names[i]: arr[i,:] for i in xrange(len(names))}
