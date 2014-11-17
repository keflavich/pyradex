from __future__ import print_function
import numpy as np
import os
from astropy import constants
from astropy import log
import astropy.table
from .. import base_class
from ..utils import ImmutableDict,unitless,lower_keys
from .. import utils

import astropy.units as u

class Fjdu(base_class.RadiativeTransferApproximator):
    def __init__(self, datapath=None, species='co',
                 density=None,
                 temperature=None,
                 tbg=2.73,
                 column=None,
                 escapeProbGeom='lvg',
                 **kwargs):

        if os.getenv('RADEX_DATAPATH') and datapath is None:
            datapath = os.getenv('RADEX_DATAPATH')

        self.datapath = os.path.dirname(datapath)
        self.species = species

        self.set_default_params()
        self.set_params(temperature=temperature, density=density,
                        column=column, geotype=escapeProbGeom, **kwargs)
        self.tbg = tbg
        from pyradex.fjdu import wrapper_my_radex
        myradex_wrapper = wrapper_my_radex.myradex_wrapper
        self._myradex = myradex_wrapper
        self._is_locked = False
        self._locked_parameter = None

    def __call__(self, return_table=True, **kwargs):

        niter = self.run_radex(**kwargs)

        if return_table:
            return self.get_table()
        else:
            return niter

    def load_datafile(self, filename=None, verbose=False):
        filename = filename or self.molpath
        self.datapath = (os.path.dirname(filename) or self.datapath)+"/"
        self.fname = os.path.basename(filename)

        nlevels, nitems, ntrans = self._myradex.config_basic(self.datapath,
                                                             self.fname,
                                                             unitless(self.tbg),
                                                             verbose)
        self.set_params(**{'n_levels': nlevels,
                           'n_item': nitems,
                           'n_transitions': ntrans})

    def run_radex(self, **kwargs):
        self.set_params(**kwargs)
        self.load_datafile()
        energies, f_occupations, data_transitions, cooling_rate = \
                self._myradex.run_one_params(**self.params)
        self._energies = u.Quantity(energies, u.K) # excitation temperature
        self._data_dict = cast_into_dic("".join(self._myradex.column_names),
                                        data_transitions)
        self._level_population = f_occupations


    _default_params = (('tkin', 0.0),
                       ('dv_CGS', 1e5),
                       ('dens_X_CGS', 0.0),
                       ('Ncol_X_CGS', 0.0),
                       ('H2_density_CGS', 0.0),
                       ('HI_density_CGS', 0.0),
                       ('oH2_density_CGS', 0.0),
                       ('pH2_density_CGS', 0.0),
                       ('HII_density_CGS', 0.0),
                       ('Electron_density_CGS', 0.0),
                       ('n_levels', 0),
                       ('n_item', 0),
                       ('n_transitions', 0),
                       ('geotype', 'lvg'),
                      )

    _keyword_map = {'temperature': 'tkin',
                    'deltav': 'dv_cgs',
                    'column': 'ncol_x_cgs'}

    _density_keyword_map = {'h2': 'h2_density_cgs',
                            'h': 'hi_density_cgs',
                            'oh2': 'oh2_density_cgs',
                            'ph2': 'ph2_density_cgs',
                            'hii': 'hii_density_cgs',
                            'e': 'electron_density_cgs',
                           }

    def set_default_params(self):
        self._params = lower_keys(dict(self._default_params))

    def set_params(self, **kwargs):
        default = lower_keys(dict(self._default_params))
        for k in kwargs:
            if kwargs[k] is None:
                continue
            if k == 'deltav':
                # deltav requires unit conversion
                self.deltav = kwargs[k]
            elif k.lower() in self._keyword_map:
                self._params[self._keyword_map[k]] = kwargs[k]
            elif k.lower() in ('density','collider_densities'):
                self.density = kwargs[k]
            elif k == 'species':
                self.species = kwargs[k]
            elif k.lower() not in default:
                raise ValueError("{0} is not a valid key.".format(k))
            else:
                self._params[k] = kwargs[k]

    @property
    def params(self):
        return lower_keys(self._params)

    @params.setter
    def params(self, value):
        if not isinstance(value, dict):
            raise TypeError('Parameters must be a dictionary.')
        self.set_params(**value)

    @property
    def density(self):

        dd = {'H2': u.Quantity(self.params['h2_density_cgs'], self._u_cc),
              'OH2': u.Quantity(self.params['oh2_density_cgs'], self._u_cc),
              'PH2': u.Quantity(self.params['ph2_density_cgs'], self._u_cc),
              'E': u.Quantity(self.params['electron_density_cgs'], self._u_cc),
              'H+': u.Quantity(self.params['hii_density_cgs'], self._u_cc),
              'H': u.Quantity(self.params['hi_density_cgs'], self._u_cc),
              'He': u.Quantity(0, self._u_cc),}
        return ImmutableDict(dd)

    @density.setter
    def density(self, value):
        if isinstance(value, dict):
            # dictionary of collider densities
            for k in value:
                if k.lower() in self._density_keyword_map:
                    key = self._density_keyword_map[k.lower()]
                    self._params[key.lower()] = value[k]
                elif k.lower() in self._density_keyword_map.values():
                    self._params[k.lower()] = value[k]
                else:
                    raise KeyError("Collider {0} not recognized.".format(k))
            self._params['dens_x_cgs'] = self.total_density.value
        else:
            self._params['dens_x_cgs'] = value
            for k in self._density_keyword_map.values():
                self._params[k] = 0.0
            self._params['h2_density_cgs'] = value

    @property
    def temperature(self):
        return u.Quantity(self.params['tkin'], u.K)

    @temperature.setter
    def temperature(self, tkin):
        if hasattr(tkin,'to'):
            tkin = unitless(u.Quantity(tkin, u.K))
        if tkin <= 0 or tkin > 1e4:
            raise ValueError('Must have kinetic temperature > 0 and < 10^4 K')
        self.set_params(tkin=tkin)

    @property
    def column_per_bin(self):
        return u.Quantity(self.params['ncol_x_cgs'], self._u_sc)

    @column_per_bin.setter
    def column_per_bin(self, col):
        if hasattr(col,'to'):
            col = unitless(u.Quantity(col, self._u_sc))
        if col < 1e5 or col > 1e25:
            raise ValueError("Extremely low or extremely high column.")
        self.set_params(ncol_x_cgs=col)

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
    def tbg(self):
        return u.Quantity(self._tbg, u.K)

    @tbg.setter
    def tbg(self, tbg):
        if hasattr(tbg, 'value'):
            self._tbg = unitless(u.Quantity(tbg, u.K))
        else:
            self._tbg = tbg

    @property
    def deltav(self):
        return u.Quantity(self.params['dv_cgs']*self._kms_to_cms, self._u_kms)

    _kms_to_cms = 1e-5
    _u_cms = u.cm/u.s

    @deltav.setter
    def deltav(self, dv):
        self.set_params(deltav=unitless(u.Quantity(dv/self._kms_to_cms,
                                                   self._u_cms)))

    @property
    def molpath(self):
        if hasattr(self,'_molpath'):
            return self._molpath

    @molpath.setter
    def molpath(self, molfile):
        if "~" in molfile:
            molfile = os.path.expanduser(molfile)
        utils.verify_collisionratefile(molfile)
        self._molpath = molfile

    @property
    def datapath(self):
        return self._datapath

    @datapath.setter
    def datapath(self, datapath):
        self._datapath = datapath

    @property
    def escapeprobProbGeom(self):
        return self._params['geotype']

    @escapeprobProbGeom.setter
    def escapeprobProbGeom(self, value):
        if value in ('lvg','spherical','slab'):
            self._params['geotype'] = value
        else:
            raise ValueError("Geometry must be spherical, slab, or lvg")
        
    _um_to_ghz = u.um.to(u.GHz, equivalencies=u.spectral())

    @property
    def frequency(self):
        return u.Quantity(self._um_to_ghz/self._data_dict['lam'], unit=u.GHz)

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
    def upperlevelnumber(self):
        return self._data_dict['iup']

    @property
    def lowerlevelnumber(self):
        return self._data_dict['ilow']

    @property
    def upperlevelpop(self):
        return self._data_dict['fup']

    @property
    def lowerlevelpop(self):
        return self._data_dict['flow']

    @property
    def source_line_brightness_temperature(self):
        return u.Quantity(self._data_dict['Tr'], u.K)
    
    #@property
    #def source_line_surfbrightness(self):
    #    return u.Quantity(self._data_dict['flux'], self._u_brightness)

    @property
    def source_brightness(self):
        return u.Quantity(self._data_dict['flux_dens'], self._u_brightness)

    @property
    def background_brightness(self):
        return u.Quantity(self._data_dict['Jback'], self._u_brightness)
    #    return self.tbg.to(self._u_brightness)

    @property
    def beta(self):
        return self._data_dict['beta']

    @property
    def statistical_weight(self):
        return self._data_dict['gup']


def cast_into_dic(col_names, arr):
    '''col_names is column_info, and arr is data_transitions'''
    names = col_names.split()
    return {names[i]: arr[i,:] for i in xrange(len(names))}
