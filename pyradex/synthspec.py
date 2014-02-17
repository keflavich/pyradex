"""
Tools to generate synthetic spectra given a table of line strengths
"""

import specutils
Spectrum = specutils.Spectrum1D
Spectrum1DLinearWCS = specutils.wcs.Spectrum1DLinearWCS

from astropy.modeling import models
from astropy import units as u
from astropy import constants as c
import numpy as np

class SyntheticSpectrum(Spectrum):
    """
    Synthetic Spectrum class - neato!
    """

    def __init__(self, minfreq, maxfreq, table,
                 linewidth=1.0*u.km/u.s,
                 npts=1000,
                 profile_function=models.Gaussian1D):
        """
        Create a synthetic spectrum from a RADEX (or DESPOTIC, eventually)
        output

        Parameters
        ----------
        minfreq: u.Quantity
        maxfreq: u.Quantity
            Minimum and maximum frequency to plot
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
        >>> import pyradex
        >>> R = pyradex.Radex(species='ch3cn')
        >>> R.run_radex()
        >>> S = pyradex.synthspec.SyntheticSpectrum(91.95*u.GHz,92*u.GHz,R.get_table())
        >>> S.plot()
        """

        self.minfreq = minfreq
        self.maxfreq = maxfreq
        self.profile_function = profile_function
        self.npts = npts

        linefreqs = u.Quantity(table['frequency'],unit=u.Unit(table['frequency'].unit))
        self.table = table[(linefreqs>minfreq) & (linefreqs<maxfreq)]
        self.width_frequency = linewidth/c.c * u.Quantity(self.table['frequency'], unit=u.Unit(self.table['frequency'].unit))

        self.wcs = Spectrum1DLinearWCS(minfreq, (maxfreq-minfreq)/npts, 0, u.Hz)

        data = self.get_profile()

        super(Spectrum,self).__init__(data=data, wcs=self.wcs,
                                      unit=u.Unit(table['T_B'].unit))


    def get_profile(self):
        
        def model(xpts):
            if isinstance(xpts,u.Quantity):
                xpts = xpts.to(u.Hz).value
            M = np.zeros_like(xpts)
            for freq,flux,width in zip(u.Quantity(self.table['frequency'],unit=u.Unit(self.table['frequency'].unit)),
                                       u.Quantity(self.table['T_B'],unit=u.Unit(self.table['T_B'].unit)),
                                       self.width_frequency):
                M += self.profile_function(flux.value, freq.to(u.Hz).value, width.to(u.Hz).value)(xpts)

            return M

        X = self.wcs(np.arange(self.npts))

        return model(X)

    def plot(self, *args, **kwargs):
        import pylab as pl

        pl.gca().set_xlabel(self.dispersion.unit.to_string())
        if hasattr(self.data,'unit'):
            pl.gca().set_ylabel(self.data.unit.to_string())
            data = self.data.value
        else:
            data = self.data

        return pl.plot(self.dispersion.value, data, *args, **kwargs)
