from __future__ import print_function
import numpy as np
import astropy.units as u
_quantity = u.Quantity
import os

from .utils import QuantityOff,ImmutableDict,unitless,grouper
from . import utils

from astropy import units as u
from astropy import constants
from astropy import log
import astropy.table

# maybe an ABC?
class RadiativeTransferApproximator(object):
    _u_brightness = (u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1 * u.sr**-1)
    _u_sc = u.cm**-2
    _u_cc = u.cm**-3

    _u_gradient = u.cm**-2 / (u.km/u.s) / u.pc
    _u_kms = u.km/u.s
    _u_cms = u.cm/u.s

    @property
    def locked_parameter(self):
        return self._locked_parameter

    def _lock_param(self, parname):
        self._locked_parameter = parname

    _all_valid_colliders = {'H2':'H2',
                            'PH2':'pH2',
                            'OH2':'oH2',
                            'E':'e',
                            'H':'H',
                            'HE':'He',
                            'H+':'H+'}

    @property
    def density(self):
        raise NotImplementedError

    @density.setter
    def density(self, collider_density):
        raise NotImplementedError


    @property
    def valid_colliders(self):
        return self._valid_colliders
    
    @property
    def total_density(self):
        """
        The total density *by number of particles* 
        The *mass density* can be dramatically different!
        """
        vc = [x.lower() for x in self.valid_colliders]
        if 'h2' in vc:
            useh2 = 1
            useoph2 = 0
        elif 'oh2' in vc or 'ph2' in vc:
            useh2 = 0
            useoph2 = 1
        else:
            # Weird case: no H2 colliders at all
            useoph2 = 0
            useh2 = 0

        weights = {'H2': useh2,
                   'PH2': useoph2,
                   'OH2': useoph2,
                   'E': 1,
                   'H': 1,
                   'He': 1,
                   'H+': 1,}

        return u.Quantity([self.density[k]*weights[k] for k in self.density]).sum()


    @property
    def mass_density(self):

        vc = [x.lower() for x in self.valid_colliders]
        if 'h2' in vc:
            useh2 = 1
            useoph2 = 0
        elif 'oh2' in vc or 'ph2' in vc:
            useh2 = 0
            useoph2 = 1
        else:
            # Weird case: no H2 colliders at all
            useoph2 = 0
            useh2 = 0

        weights = {'H2': 2*useh2,
                   'PH2': 2*useoph2,
                   'OH2': 2*useoph2,
                   'E': 1/1836.,
                   'H': 1,
                   'He': 4,
                   'H+': 1,}

        return np.sum( (self.density[k]*weights[k] for k in self.density)
                     )*constants.m_p


    @property
    def opr(self):
        return self.density['OH2']/self.density['PH2']

    @property
    def oprh2(self):
        return self.opr

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, species):
        if hasattr(self,'_species') and self._species == species:
            return
        self._species = species
        try:
            self.molpath = os.path.join(self.datapath,species+'.dat')
        except IOError:
            log.warn("Did not find data file for species {0} "
                     "in path {1}.  Downloading it.".format(species,
                                                            self.datapath))
            utils.get_datafile(species, self.datapath)
            self.molpath = os.path.join(self.datapath,species+'.dat')

        self._valid_colliders = utils.get_colliders(self.molpath)
        vc = [x.lower() for x in self._valid_colliders]
        if 'h2' in vc and ('oh2' in vc or 'ph2' in vc):
            log.warn("oH2/pH2 and h2 are both in the datafile: "
                     "The resulting density/total density are invalid.")

    @property
    def molpath(self):
        raise NotImplementedError

    @molpath.setter
    def molpath(self, molfile):
        raise NotImplementedError


    @property
    def datapath(self):
        return os.path.dirname(self.molpath)

    @datapath.setter
    def datapath(self, radat):
        raise NotImplementedError

    @property
    def escapeProbGeom(self):
        raise NotImplementedError

    @escapeProbGeom.setter
    def escapeProbGeom(self, escapeProbGeom):
        raise NotImplementedError
        

    @property
    def column(self):
        return self.column_per_bin

    @column.setter
    def column(self, value):
        self.column_per_bin = value


    @property
    def column_per_kms_perpc(self):
        return self.column_per_bin / self.deltav


    @column_per_kms_perpc.setter
    def column_per_kms_perpc(self, cddv):

        cddv = u.Quantity(cddv, self._u_gradient)

        self.column_per_bin = cddv * u.Quantity(self.deltav, self._u_kms) * self.length()

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
        return self._deltav


    @deltav.setter
    def deltav(self, dv):
        self._deltav = u.Quantity(dv, self._u_kms)

    @property
    def length(self):
        """ Hard-coded, assumed length-scale """
        return u.Quantity(1, u.pc)

    @property
    def tbg(self):
        raise NotImplementedError

    def _validate_colliders(self):
        """
        Check whether the density of at least one collider in the associated
        LAMDA data file is nonzero
        """
        valid_colliders = [x.lower() for x in self.valid_colliders]

        density = self.density

        OK = False
        matched_colliders = []
        for collider in valid_colliders:
            if unitless(density[self._all_valid_colliders[collider.upper()]]) > 0:
                OK = True
                matched_colliders.append(collider.lower())

        if not OK:
            raise ValueError("The colliders in the data file {0} ".format(self.molpath)
                             + "have density 0.")

        bad_colliders = []
        for collider in density:
            if (unitless(density[collider]) > 0
                and (collider.lower() not in valid_colliders)):
                if (collider.lower() in ('oh2','ph2') and 'h2' in
                    matched_colliders):
                    # All is OK: we're allowed to have mismatches of this sort
                    continue
                elif (collider.lower() == 'h2' and ('oh2' in matched_colliders
                                                    or 'ph2' in
                                                    matched_colliders)):
                    # again, all OK
                    continue
                bad_colliders.append(collider)
                OK = False

        if not OK:
            raise ValueError("There are colliders with specified densities >0 "
                             "that do not have corresponding collision rates."
                             "  The bad colliders are {0}".format(bad_colliders))


    @property
    def source_area(self):
        if hasattr(self, '_source_area'):
            return self._source_area

    @source_area.setter
    def source_area(self, source_area):
        self._source_area = source_area

    @property
    def source_line_surfbrightness(self):
        return self.source_brightness - self.background_brightness

    def line_brightness_temperature(self,beamsize):
        """
        Return the line surface brightness in kelvins for a given beam area
        (Assumes the frequencies are rest frequencies)
        """
        #return (self.line_flux * beamsize)
        # because each line has a different frequency, have to loop it
        try:
            return u.Quantity([x.to(u.K, u.brightness_temperature(beamsize, f)).value
                               for x,f in zip(self.line_flux_density,self.frequency)
                               ],
                              unit=u.K)
        except AttributeError as ex:
            raise NotImplementedError("line brightness temperature is not implemented "
                                      "without reference to astropy units yet")

    @property
    def source_line_brightness_temperature(self):
        """
        The surface brightness of the source assuming it is observed with a
        beam matched to its size and it has ff=1

        (this is consistent with the online RADEX calculator)
        """
        #return (self.line_flux * beamsize)
        # because each line has a different frequency, have to loop it
        return ((self.source_line_surfbrightness*u.sr).
                 to(u.K, u.brightness_temperature(1*u.sr,
                                                  self.frequency)))

    @property
    def T_B(self):
        return self.source_line_brightness_temperature


    @property
    def background_brightness(self):
        raise NotImplementedError

    @property
    def flux_density(self):
        """
        Convert the source surface brightness to a flux density by specifying
        the emitting area of the source (in steradian-equivalent units)

        This is the non-background-subtracted version
        """

        if not self.source_area:
            raise AttributeError("Need to specify a source area in order to compute the flux density")

        return self.source_brightness * self.source_area

    @property
    def line_flux_density(self):
        """
        Background-subtracted version of flux_density
        """

        if not self.source_area:
            raise AttributeError("Need to specify a source area in order to compute the flux density")

        return self.source_line_surfbrightness * self.source_area


    @property
    def source_brightness(self):
        """
        RADEX compat?  (check)
        """

        raise NotImplementedError

    @property
    def source_brightness_beta(self):

        raise NotImplementedError

    @property
    def beta(self):
        raise NotImplementedError

    def get_table(self):
        columns = [
            astropy.table.Column(name='Tex', data=self.tex, unit=u.K),
            astropy.table.Column(name='tau', data=self.tau, unit=''),
            astropy.table.Column(name='frequency', data=self.frequency,
                                 unit=u.GHz),
            astropy.table.Column(name='upperstateenergy',
                                 data=self.upperstateenergy, unit=u.K),
            astropy.table.Column(name='upperlevel',
                                 data=self.upperlevelnumber,
                                 unit=''),
            astropy.table.Column(name='lowerlevel',
                                 data=self.lowerlevelnumber,
                                 unit=''),
            astropy.table.Column(name='upperlevelpop',
                                 data=self.upperlevelpop,
                                 unit=''),
            astropy.table.Column(name='lowerlevelpop',
                                 data=self.lowerlevelpop,
                                 unit=''),
            astropy.table.Column(name='brightness',
                                 data=self.source_line_surfbrightness),
            astropy.table.Column(name='T_B', data=self.T_B), # T_B is pre-masked
        ]
        if self.source_area:
            columns.append(astropy.table.Column(name='flux',data=self.line_flux_density[mask]))

        T = astropy.table.Table(columns)

        return T

    def get_synthspec(self, fmin, fmax, npts=1000, **kwargs):
        """
        Generate a synthetic spectrum of the selected molecule over the
        specified frequency range.  This task is good for quick-looks but has a
        lot of overhead for generating models and should not be used for
        fitting (unless you have a conveniently small amount of data)

        Parameters
        ----------
        fmin : `~astropy.units.Quantity`
        fmax : `~astropy.units.Quantity`
            Frequency-equivalent quantity
        """
        wcs = synthspec.FrequencyArray(fmin, fmax, npts)
        S = synthspec.SyntheticSpectrum.from_RADEX(wcs, self, **kwargs)

        return S
