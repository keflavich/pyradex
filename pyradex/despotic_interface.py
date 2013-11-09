import despotic
import numpy as np
import os
from collections import defaultdict
from astropy import units as u

class Despotic(object):
    """
    A class meant to be similar to RADEX in terms of how it is called but to
    use Despotic instead of RADEX as the backend
    """
    def __call__(self, **kwargs):
        # reset the parameters appropriately
        self.__init__(**kwargs)
        return self.lineLum(self.species)

    def __init__(self,
                 collider_densities={'ph2':990,'oh2':10},
                 temperature=30,
                 species='co',
                 datapath=None,
                 hcolumn=1e21,
                 abundance=1e-5,
                 #column=1e13,
                 tbackground=2.7315,
                 deltav=1.0,
                 length=3.085677581467192e+18,
                 method='lvg',
                 outfile='radex.out',
                 logfile='radex.log',
                 debug=False,
                 ):
        """
        Interface to DESPOTIC

        Parameters
        ----------
        collider_densities: dict
            Dictionary giving the volume densities of the collider(s) in units of
            cm^-3.  Valid entries are h2,oh2,ph2,e,He,H,H+.  The keys are
            case-insensitive.
        temperature: float
            Local gas temperature in K
        species: str
            A string specifying a valid chemical species.  This is used to look
            up the specified molecule
        hcolumn: float
            The total column density of hydrogen.
        abundance: float
            The molecule's abundance relative to H.
        tbackground: float
            Background radiation temperature (e.g., CMB)
        deltav: float
            The FWHM line width (really, the single-zone velocity width to
            scale the column density by: this is most sensibly interpreted as a
            velocity gradient (dv/length))
        length: float
            For abundance calculations, the line-of-sight size scale to
            associate with a velocity bin.  Typically, assume 1 km/s/pc,
            and thus set deltav=1, length=3.08e18
        datapath: str
            Path to the molecular data files
        """
        
        self.cloud = despotic.cloud()
        self.cloud.nH = float(np.sum([collider_densities[k]*2 if 'h2' in k.lower()
                                      else collider_densities[k]
                                      for k in collider_densities]))

        for k in collider_densities.keys():
            collider_densities[k.lower()] = collider_densities[k]

        if 'ph2' in collider_densities:
            self.cloud.comp.xpH2 = collider_densities['ph2'] / self.cloud.nH
        if 'oh2' in collider_densities:
            self.cloud.comp.xoH2 = collider_densities['oh2'] / self.cloud.nH

        self.cloud.Td = temperature
        self.cloud.Tg = temperature
        self.cloud.dVdr = (deltav*u.km/u.s / (length*u.pc)).to(1/u.s).value
        self.cloud.colDen = hcolumn


        if tbackground > 2.7315:
            self.cloud.rad.TradDust = tbackground

        self.species = species
        if datapath is None:
            emitterFile = species+'.dat'
        else:
            emitterFile = os.path.expanduser(os.path.join(datapath, species+'.dat'))
        self.cloud.addEmitter(species, abundance, emitterFile=emitterFile)

    def lineLum(self):
        return self.cloud.lineLum(self.species)

    @property
    def density(self):
        d = {'H2':self.cloud.comp.xH2*self.nH,
             'oH2':self.cloud.comp.xoH2*self.nH,
             'pH2':self.cloud.comp.xpH2*self.nH,
             'e':self.cloud.comp.xe*self.nH,
             'H':self.cloud.comp.xHI*self.nH,
             'He':self.cloud.comp.xHe*self.nH,
             'H+':self.cloud.comp.xHplus*self.nH}
        if u:
            for k in d:
                d[k] = d[k] * u.cm**-3
        return d

    @property
    def nH(self):
        return self.cloud.nH

    @density.setter
    def density(self, collider_density):
        collider_densities = defaultdict(lambda: 0)
        for k in collider_density:
            collider_densities[k.upper()] = collider_density[k]

        if 'OH2' in collider_densities:
            if not 'PH2' in collider_densities:
                raise ValueError("If o-H2 density is specified, p-H2 must also be.")
            self.cloud.comp.xH2 = (collider_densities['OH2'] + collider_densities['PH2'])/self.nH
            self.cloud.comp.xoH2 = (collider_densities['OH2'])/self.nH
            self.cloud.comp.xpH2 = (collider_densities['PH2'])/self.nH
        elif 'H2' in collider_densities:
            warnings.warn("Using a default ortho-to-para ratio (which "
                          "will only affect species for which independent "
                          "ortho & para collision rates are given)")
            self.cloud.comp.xH2 = collider_densities['H2']/self.nH

            T = self.temperature.value if hasattr(self.temperature,'value') else self.temperature
            if T > 0:
                opr = min(3.0,9.0*np.exp(-170.6/T))
            else:
                opr = 3.0
            fortho = opr/(1+opr)
            self.cloud.comp.xoH2 = collider_densities['H2']*fortho/self.nH
            self.cloud.comp.xpH2 = collider_densities['H2']*(1-fortho)/self.nH

        self.cloud.comp.xe = collider_densities['E']/self.nH
        self.cloud.comp.xHI = collider_densities['H']/self.nH
        self.cloud.comp.xHe = collider_densities['HE']/self.nH
        self.cloud.comp.xHplus = collider_densities['H+']/self.nH

        # skip H2 when computing by assuming OPR correctly distributes ortho & para
        self.radex.cphys.totdens = self.radex.cphys.density[1:].sum()

    @property
    def temperature(self):
        return self.cloud.Tg
