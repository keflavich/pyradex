try:
    import despotic
except ImportError:
    despotic=False
import numpy as np
import warnings
import os
from collections import defaultdict
from .utils import united,uvalue

try:
    from astropy import units as u
    from astropy import constants
    import astropy.table
except ImportError:
    u = False

class Despotic(object):
    """
    A class meant to be similar to RADEX in terms of how it is called but to
    use Despotic instead of RADEX as the backend
    """

    # escapeprobabilty geometry names
    _epgdict = {'lvg':'LVG','sphere':'sphere','slab':'slab'}

    def __call__(self, **kwargs):
        # reset the parameters appropriately
        self.__init__(**kwargs)
        return self.lineLum()

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
                 escapeProbGeom='lvg',
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
            The molecule's abundance relative to H (NOT H2 as is normally done!).
        tbackground: float
            Background radiation temperature (e.g., CMB)
        deltav: float
            The FWHM line width (really, the single-zone velocity width to
            scale the column density by: this is most sensibly interpreted as a
            velocity gradient (dv/dR))
        sigmaNT: float
            Nonthermal velocity dispersion
            (this is strictly ignored - deltav IS sigmant)
        datapath: str
            Path to the molecular data files
        """

        if not despotic:
            raise ImportError("Despotic could not be imported.  Please check that it is installed.")
        
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

        self.cloud.Td = uvalue(temperature,u.K)
        self.cloud.Tg = uvalue(temperature,u.K)
        self.cloud.dust.sigma10 = 0.0

        self.cloud.colDen = uvalue(hcolumn,u.cm**-2)


        if uvalue(tbackground,u.K) > 2.7315:
            self.cloud.rad.TradDust = uvalue(tbackground,u.K)

        self.species = species
        if datapath is None:
            emitterFile = species+'.dat'
        else:
            emitterFile = os.path.expanduser(os.path.join(datapath, species+'.dat'))
        self.cloud.addEmitter(species, abundance, emitterFile=emitterFile)


        self.cloud.comp.computeDerived(self.cloud.nH)

        self.escapeProbGeom = escapeProbGeom

        self.deltav = deltav

    @property
    def deltav(self):
        try:
            return united(self._dv,u.km/u.s)
        except u.UnitsError:
            return united(self._dv,u.km/u.s/u.pc)

    @deltav.setter
    def deltav(self, deltav):

        if self.escapeProbGeom == 'LVG':
            self._dv = united(self.cloud.dVdr,u.s**-1)
            # See notes.rst: DESPOTIC must have a different dVdR to get the same results as RADEX
            # 1.0645 / sqrt(8*log(2)) = sqrt(2 * pi) / (8*log(2))
            self.cloud.dVdr = uvalue(united(deltav,u.km/u.s/u.pc),u.s**-1) * np.sqrt(8*np.log(2)) * 2
            # (2*np.pi)**0.5/(8*np.log(2))

        else:
            self._dv = united(self.cloud.sigmaNT,u.km/u.s)

            FWHM = united(deltav,u.km/u.s)
            self.sigmaTot = FWHM/np.sqrt(8.0*np.log(2))
            self.cloud.sigmaNT = uvalue(np.sqrt(self.sigmaTot**2 -
                                                self.cs**2/self.cloud.emitters[self.species].data.molWgt),
                                        u.km/u.s)

    @property
    def cs(self):
        return np.sqrt(constants.k_B*united(self.cloud.Tg,u.K)/(self.cloud.comp.mu*constants.m_p)).to(u.km/u.s)

    def lineLum(self, **kwargs):
        if 'escapeProbGeom' not in kwargs:
            kwargs['escapeProbGeom'] = self.escapeProbGeom
        return self.cloud.lineLum(self.species, **kwargs)

    @property
    def escapeProbGeom(self):
        return self._epg

    @escapeProbGeom.setter
    def escapeProbGeom(self, escapeProbGeom):
        mdict = self._epgdict
        if escapeProbGeom.lower() not in mdict:
            raise ValueError("Invalid escapeProbGeom, must be one of "+",".join(mdict.values()))
        self._epg = mdict[escapeProbGeom]

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

    @nH.setter
    def nH(self, nh):
        self.cloud.nH = nh

    @property
    def nH2(self):
        return self.cloud.nH/2.

    @nH2.setter
    def nH2(self, nh2):
        self.nh = nh2*2.


    @density.setter
    def density(self, collider_density):
        collider_densities = defaultdict(lambda: 0)
        for k in collider_density:
            collider_densities[k.upper()] = collider_density[k]

        if 'OH2' in collider_densities:
            if not 'PH2' in collider_densities:
                raise ValueError("If o-H2 density is specified, p-H2 must also be.")
            collider_densities['H2'] = (collider_densities['OH2'] + collider_densities['PH2'])
        elif 'H2' in collider_densities:
            warnings.warn("Using a default ortho-to-para ratio (which "
                          "will only affect species for which independent "
                          "ortho & para collision rates are given)")

            T = self.temperature.value if hasattr(self.temperature,'value') else self.temperature
            if T > 0:
                opr = min(3.0,9.0*np.exp(-170.6/T))
            else:
                opr = 3.0
            fortho = opr/(1+opr)
            collider_densities['OH2'] = collider_densities['H2']*fortho
            collider_densities['PH2'] = collider_densities['H2']*(1-fortho)

        total_density = np.sum([collider_densities[x] * (2 if '2' in x else 1)
                                for x in (['OH2','PH2','H','E','HE','H+'])])
        self.nH = total_density

        self.cloud.comp.xH2 = collider_densities['H2']/self.nH
        self.cloud.comp.xoH2 = (collider_densities['OH2'])/self.nH
        self.cloud.comp.xpH2 = (collider_densities['PH2'])/self.nH

        self.cloud.comp.xe = collider_densities['E']/self.nH
        self.cloud.comp.xHI = collider_densities['H']/self.nH
        self.cloud.comp.xHe = collider_densities['HE']/self.nH
        self.cloud.comp.xHplus = collider_densities['H+']/self.nH

        self.cloud.comp.computeDerived()

    @property
    def temperature(self):
        return self.cloud.Tg

    def get_table(self,**kwargs):

        D = self.lineLum(**kwargs)

        names = D[0].keys()
        T = astropy.table.Table(names=names,dtypes=[type(D[0][k]) for k in names])
    
        for row in D:
            T.add_row([row[k] for k in names])

        return T
