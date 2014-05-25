"""
Tools to generate synthetic spectra given a table of line strengths
"""
import numpy as np

from astropy import logging
from astropy.modeling import models
from astropy import units as u
from astropy import constants as c

try:
    import specutils
    Spectrum = specutils.Spectrum1D
    Spectrum1DLinearWCS = specutils.wcs.Spectrum1DLinearWCS

    class FrequencyArray(Spectrum1DLinearWCS):
        """
        Generate a frequency array given a minimum frequency, maximum frequency,
        and # of points

        Parameters
        ----------
        minfreq: u.Quantity
        maxfreq: u.Quantity
            Minimum and maximum frequency to plot
        npts: int
            The number of spectral points to include
        """
        def __init__(self, minfreq, maxfreq, npts):
            self.minfreq = minfreq
            self.maxfreq = maxfreq
            self.npts = npts
            self.delta = (maxfreq-minfreq)/npts
            self.unit = u.Hz
            self.refpix = 0
            self.refval = minfreq
            #super(FrequencyArray, self).__init__(minfreq, (maxfreq-minfreq)/npts, 0, u.Hz)

        def __call__(self, pixel_indices):
            return (self.refval + self.delta * (pixel_indices -
                                                self.refpix))

        @property
        def lookuptable(self):
            if not hasattr(self, '_lookuptable'):
                self._lookuptable = self(np.arange(npts))
            return self._lookuptable


    class SyntheticSpectrum(Spectrum):
        """
        Synthetic Spectrum class - neato!
        """

        def __init__(self, wcs, species, linewidth):
            self.wcs = wcs
            self.species = species
            self.linewidth = linewidth

        @classmethod
        def from_table(cls, wcs, table, species,
                       linewidth=1.0*u.km/u.s,
                       profile_function=models.Gaussian1D):
            """
            Create a synthetic spectrum from a RADEX (or DESPOTIC, eventually)
            output

            Parameters
            ----------
            wcs: SpectralWCS
                A spectral world coordinate system.  You can generate one with
                FrequencyArray or specutils.wcs.Spectrum1DLookupWCS
            table: astropy.Table
                Result of the RADEX query (from R.get_table())
            linewidth: u.Quantity (km/s)
                The width of the line to plot
            npts: int
                The number of spectral points to include
            profile_function: astropy.modeling.model
                The model function to use.  Must accept, in order:
                
                 * flux (peak)
                 * frequency center (Hz)
                 * frequency width (Hz)
                 

            Examples
            --------
            >>> from pyradex import Radex,synthspec
            >>> R = Radex(species='ch3cn')
            >>> R.run_radex()
            >>> wcs = synthspec.FrequencyArray(91.95*u.GHz, 92*u.GHz, npts=1000)
            >>> S = synthspec.SyntheticSpectrum.from_table(wcs, R.get_Table(),
            ...                                 species='ch3cn')
            >>> S.plot()
            """

            self = cls(wcs, species, linewidth)

            self.profile_function = profile_function

            self.species = species
            self.linewidth = linewidth
            self.wcs = wcs
            if hasattr(wcs,'minfreq'):
                self.minfreq,self.maxfreq = wcs.minfreq,wcs.maxfreq
            else:
                self.minfreq,self.maxfreq = wcs.min(),wcs.max()

            linefreqs = u.Quantity(table['frequency'],
                                   unit=u.Unit(table['frequency'].unit))
            self.table = table[(linefreqs>self.minfreq) & (linefreqs<self.maxfreq)]
            self.linefreqs = linefreqs
            self.width_frequency = (linewidth/c.c *
                                    u.Quantity(self.table['frequency'],
                                               unit=u.Unit(self.table['frequency'].unit)))
            self.T_B = self.table['T_B']

            data = self.get_profile()

            super(Spectrum,self).__init__(data=data, wcs=self.wcs,
                                          unit=u.Unit(table['T_B'].unit))

            return self

        @classmethod
        def from_RADEX(cls, wcs, rad,
                       linewidth=1.0*u.km/u.s,
                       profile_function=models.Gaussian1D):
            """
            Create a synthetic spectrum from a RADEX class

            Parameters
            ----------
            wcs: SpectralWCS
                A spectral world coordinate system.  You can generate one with
                FrequencyArray or specutils.wcs.Spectrum1DLookupWCS
            rad: pyradex.Radex instance
                Result of the RADEX query
            linewidth: u.Quantity (km/s)
                The width of the line to plot
            npts: int
                The number of spectral points to include
            profile_function: astropy.modeling.model
                The model function to use.  Must accept, in order:
                
                 * flux (peak)
                 * frequency center (Hz)
                 * frequency width (Hz)
                 

            Examples
            --------
            >>> from pyradex import Radex,synthspec
            >>> R = Radex(species='ch3cn')
            >>> R.run_radex()
            >>> wcs = synthspec.FrequencyArray(91.95*u.GHz, 92*u.GHz, npts=1000)
            >>> S = synthspec.SyntheticSpectrum.from_RADEX(wcs, R)
            >>> S.plot()
            """

            self = cls(wcs, rad.species, linewidth)

            self.profile_function = profile_function

            self.wcs = wcs
            if hasattr(wcs,'minfreq'):
                self.minfreq,self.maxfreq = wcs.minfreq,wcs.maxfreq
            else:
                self.minfreq,self.maxfreq = wcs.min(),wcs.max()

            self.rad = rad
            linefreqs = rad.frequency
            linefreq_mask = (linefreqs>self.minfreq) & (linefreqs<self.maxfreq)
            included_frequencies_mask = linefreq_mask[rad.inds_frequencies_included]
            self.linefreqs = linefreqs[linefreq_mask]
            self.T_B = rad.T_B[included_frequencies_mask]
            self.width_frequency = (linewidth/c.c * self.linefreqs)

            data = self.get_profile()

            super(Spectrum,self).__init__(data=data, wcs=self.wcs,
                                          unit=u.Unit(rad.T_B.unit))

            return self

        def get_profile(self, velocity_offset=0*u.km/u.s):
            
            def model(xpts):
                if isinstance(xpts,u.Quantity):
                    xpts = xpts.to(u.Hz).value
                M = np.zeros_like(xpts)
                freqs = self.linefreqs + (self.linefreqs*velocity_offset/c.c)
                for freq,flux,width in zip(freqs,
                                           self.T_B,
                                           self.width_frequency):
                    M += self.profile_function(flux.value, freq.to(u.Hz).value,
                                               width.to(u.Hz).value)(xpts)

                return M

            X = self.wcs(np.arange(self.wcs.npts))

            return model(X)

        def plot(self, update_data=False, *args, **kwargs):
            import pylab as pl

            if update_data:
                self.data = self.get_profile()

            pl.gca().set_xlabel(self.dispersion.unit.to_string())
            if hasattr(self.data,'unit'):
                pl.gca().set_ylabel(self.data.unit.to_string())
                data = self.data.value
            else:
                data = self.data

            return pl.plot(self.dispersion.value, data, *args, **kwargs)

        def __call__(self, linewidth=None, velocity_offset=0*u.km/u.s, **kwargs):
            """
            Return a synthetic spectrum created by calling RADEX.  Parameters
            are passed to pyradex.Radex (except linewidth)

            Parameters
            ----------
            linewidth: u.Quantity (km/s)
                The width of the line to plot

            Examples
            --------
            >>> from pyradex import Radex,synthspec
            >>> radex_pars = dict(temperature=20, column=1e13,
            ...                   abundance=10**-8.5,
            ...                   collider_densities={'H2':1e4})
            >>> R = Radex(species='oh2co-h2', **radex_pars)
            >>> R.run_radex()
            >>> wcs = synthspec.FrequencyArray(4.828*u.GHz, 4.830*u.GHz, npts=1000)
            >>> S = synthspec.SyntheticSpectrum.from_RADEX(wcs, R)
            >>> S.plot()
            >>> radex_pars['temperature'] = 50
            >>> S2 = S(velocity_offset=2*u.km/u.s, **radex_pars)
            >>> S2.plot()
            """

            from .core import Radex

            rad = Radex(species=self.species, **kwargs)

            if linewidth is None:
                linewidth = self.linewidth
            else:
                self.linewidth = linewidth

            self.rad = rad
            linefreqs = rad.frequency
            linefreq_mask = (linefreqs>self.minfreq) & (linefreqs<self.maxfreq)
            self.linefreqs = linefreqs[linefreq_mask]
            self.T_B = rad.T_B[linefreq_mask]
            self.width_frequency = (linewidth/c.c * self.linefreqs)

            self.data = self.get_profile()

            return self
except ImportError as exc:
    logging.warn("Failed to import synthetic spectrum package: "+exc.args[0])
