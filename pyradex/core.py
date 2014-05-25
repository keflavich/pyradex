from __future__ import print_function
import subprocess
import tempfile
import numpy as np
import warnings
import astropy.units as u
from collections import defaultdict
import itertools
import os

from . import utils
from . import synthspec

try:
    from astropy import units as u
    from astropy import constants
    import astropy.table
except ImportError:
    u = False

__all__ = ['pyradex', 'write_input', 'parse_outfile', 'call_radex', 'Radex']


# silly tool needed for fortran misrepresentation of strings
# http://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks
def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)
        

def pyradex(executable='radex', minfreq=100, maxfreq=130,
            collider_densities={'H2':1}, debug=False, delete_tempfile=True,
            **kwargs):
    """
    Get the radex results for a set of input parameters


    Parameters
    ----------
    executable : str
        Full path to the RADEX executable
    minfreq : float
        Lowest frequency line to store, in GHz
        (note: any astropy.unit spectroscopic unit is also allowed)
    maxfreq : float
        Highest frequency line to store
    collider_densities : dict
        Collider names and their number densities
        If the molecule specified has both o-H2 and p-H2, you will get a
        WARNING if you specify 'H2'
        An ortho/para example:
        collider_densities = {'oH2':900, 'pH2':100} 
        which will yield H2 = 1000

    See write_input for additional parameters

    Returns
    -------
    An astropy table containing the RADEX returns

    .. WARNING:: If RADEX spits out *******, it will be replaced with -999
    """
    warnings.warn("pyradex is deprecated: it uses very slow hard disk file i/o.  Use pyradex.Radex instead if you can.")

    infile,outfile = write_input(minfreq=minfreq, maxfreq=maxfreq,
            delete_tempfile=delete_tempfile,
            collider_densities=collider_densities, **kwargs)

    logfile = call_radex(executable, infile.name, debug=debug,
                         delete_tempfile=delete_tempfile)

    check_logfile(logfile.name)

    data = parse_outfile(outfile.name)

    if debug:
        with open(infile.name,'r') as inf:
            print("Input:")
            print(inf.read())
        with open(outfile.name,'r') as out:
            print("Output:")
            print(out.read())

    infile.close()
    outfile.close()
    logfile.close()

    return data

def check_logfile(logfilename):
    with open(logfilename,'r') as f:
        if "Warning: Assuming thermal o/p ratio" in f.read():
            warnings.warn("Assumed thermal o/p ratio since only H2 was given but collider file has o- and p- H2")

def write_input(temperature=10, column=1e12, collider_densities={'H2':1},
        bw=0.01, tbg=2.73, species='co', velocity_gradient=1.0, minfreq=1,
        maxfreq=10, delete_tempfile=True):
    """
    Write radex.inp file parameters

    Parameters
    ----------
    temperature : float
        Kinetic temperature (K)
    collider_densities : dict
        Collider names and their number densities
    column : float
        column density of the molecule
    species : str
        Name of the molecule (specifically, the prefix for the file name, e.g.
        for "co.dat", species='co').  Case sensitive!
    tbg : float
        Temperature of the background radiation (e.g. CMB)
    velocity_gradient : float
        Velocity gradient per pc in km/s
    """

    if hasattr(minfreq, 'unit'):
        minfreq = minfreq.to('GHz',u.spectral()).value
    if hasattr(maxfreq, 'unit'):
        maxfreq = maxfreq.to('GHz',u.spectral()).value

    infile = tempfile.NamedTemporaryFile(mode='w', delete=delete_tempfile)
    outfile = tempfile.NamedTemporaryFile(mode='w', delete=delete_tempfile)
    infile.write(species+'.dat\n')
    infile.write(outfile.name+'\n')
    infile.write(str(minfreq)+' '+str(maxfreq)+'\n')
    infile.write(str(temperature)+'\n')

    # RADEX doesn't allow densities < 1e-3
    for k in collider_densities.keys():
        if collider_densities[k] < 1e-3:
            collider_densities.pop(k)

    infile.write('%s\n' % len(collider_densities))
    for name,dens in collider_densities.iteritems():
        infile.write('%s\n' % name)
        infile.write(str(dens)+'\n')
    infile.write(str(tbg)+'\n')
    infile.write(str(column)+'\n')
    infile.write(str(velocity_gradient)+'\n')
    # end the input file
    infile.write('0\n')
    infile.flush()
    return infile,outfile

def call_radex(executable, inpfilename, debug=False, delete_tempfile=True):

    logfile = tempfile.NamedTemporaryFile(mode='w', delete=delete_tempfile)
    cmd = '{radex} < {inpfile} > {logfile}'.format(
        radex=executable,
        inpfile=inpfilename,
        logfile=logfile.name)
    if debug:
        print("Command:",cmd)

    result = subprocess.call(cmd, shell=True)
    if result != 0:
        print("RADEX returned error code %i" % result)
        with open(logfile.name,'r') as f:
            print(f.read())

    return logfile


header_names = ['J_up','J_low','E_UP','FREQ', 'WAVE', 'T_EX', 'TAU', 'T_R', 'POP_UP', 'POP_LOW', 'FLUX_Kkms',   'FLUX_Inu']
header_units = [None,       None, u.K,   u.GHz,  u.um,   u.K,    None,  u.K,   None,     None,     u.K*u.km/u.s, u.erg/u.cm**2/u.s]
dtypes       = [str,   str,    float, float,  float,  float,  float, float,  float,    float,     float,         float]

def parse_outfile(filename):
    with open(filename,'r') as f:
        alllines = f.readlines()
        header = {L.split(":")[0][2:].strip():L.split(":")[1].strip()
                for L in alllines
                if L[0]=='*'}
        lines = [L.replace("--","  ") for L in alllines 
                if (L[0] != '*' 
                    and 'iterat' not in L 
                    and 'GHz' not in L 
                    and 'TAU' not in L)]
    data_list = [[x if '*' not in x else '-999' for x in L.split()] for L in lines]
    if len(data_list) == 0:
        raise ValueError("No lines included?")
    data_in_columns = map(list,zip(*data_list))
    columns = [astropy.table.Column(data=C, name=name.lower(), units=units, dtype=dtype) 
            for C,name,units,dtype in zip(data_in_columns, header_names, header_units, dtypes)]
    data = astropy.table.Table(columns, meta=header)
    return data

class Radex(object):

    def __call__(self, return_table=True, **kwargs):
        # reset the parameters appropriately
        self.__init__(**kwargs)
        niter = self.run_radex()

        if return_table:
            return self.get_table()
        else:
            return niter

    def __init__(self,
                 collider_densities={'ph2':990,'oh2':10},
                 temperature=30,
                 species='co',
                 column=None,
                 h2column=None,
                 tbackground=2.7315,
                 deltav=1.0,
                 abundance=None,
                 datapath=None,
                 escapeProbGeom='lvg',
                 outfile='radex.out',
                 logfile='radex.log',
                 debug=False,
                 mu=2.8,
                 source_area=None,
                 ):
        """
        Direct wrapper of the radex FORTRAN code

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
        column: float
            The column density of the molecule of interest.  If the column is
            specified, the abundance is equal to (column/(dv*total
            density*length)).
        h2column: float
            The column of h2. 
        abundance: float
            The molecule's abundance relative to the total collider density in
            each velocity bin, i.e. column = abundance * density * length * dv.
            If both abundance and column are specified, abundance is ignored.
        tbackground: float
            Background radiation temperature (e.g., CMB)
        deltav: float
            The FWHM line width (really, the single-zone velocity width to
            scale the column density by: this is most sensibly interpreted as a
            velocity gradient (dv/length))
        datapath: str
            Path to the molecular data files.  If it is not specified, defaults
            to the current directory, OR the shell variable RADEX_DATAPATH if
            it is specified.
        outfile: str
            Output file name
        logfile: str
            Log file name
        escapeProbGeom: 'lvg','sphere','slab'
            Which escape probability method to use
        mu: float
            Mean mass per particle in AMU.  Set to 2.8 for H2+Helium mix
        source_area: float / unit
            The emitting area of the source on the sky in steradians
        """
        from pyradex.radex import radex
        self.radex = radex

        self.mu = mu

        if os.getenv('RADEX_DATAPATH') and datapath is None:
            datapath = os.getenv('RADEX_DATAPATH')

        if datapath is not None:
            self.datapath = datapath
            if self.datapath != os.path.expanduser(datapath):
                raise ValueError("Data path %s was not successfully stored; instead %s was." % (datapath,self.datapath))
        self.molpath = os.path.join(self.datapath,species+'.dat')
        if self.molpath == '':
            raise ValueError("Must set a species name.")
        if not os.path.exists(self.molpath):
            raise ValueError("Must specify a valid path to a molecular data file else RADEX will crash."
                             "  Current path is {0}".format(self.molpath))

        self.density = collider_densities
        self.species = species

        self.outfile = outfile
        self.logfile = logfile
        self.escapeProbGeom = escapeProbGeom

        self.deltav = deltav
        self._set_parameters()

        if None not in (column,abundance,h2column) and column/h2column != abundance:
            raise ValueError("Cannot specify column, h2column, and abundance unless they are consistent.")

        if sum(x is None for x in (column,abundance,h2column)) >= 2:
            raise ValueError("Must specify at least two of: h2column, column, and abundance")

        if abundance:
            self._abundance = abundance
        else:
            self._abundance = column/h2column
        
        if h2column:
            self.h2column = h2column

        if column:
            self.column = column

        self.temperature = temperature
        self.tbg = tbackground

        self.debug = debug

        self.source_area = source_area

    def _set_parameters(self):

        #self.radex.cphys.cdmol = self.column
        #self.radex.cphys.tkin = self.temperature
        if u:
            self.radex.cphys.deltav = self.deltav.to(u.cm/u.s).value
        else:
            self.radex.cphys.deltav = self.deltav*1e5

        # these parameters are only used for outputs and therefore can be ignored
        self.radex.freq.fmin = 0
        self.radex.freq.fmax = 1e10

        self.miniter = 10
        self.maxiter = 200

    _valid_colliders = ['H2','PH2','OH2','E','H','HE','H+']

    @property
    def density(self):
        d = {'H2':self.radex.cphys.density[0],
             'pH2':self.radex.cphys.density[1],
             'hH2':self.radex.cphys.density[2],
             'e':self.radex.cphys.density[3],
             'H':self.radex.cphys.density[4],
             'He':self.radex.cphys.density[5],
             'H+':self.radex.cphys.density[6]}
        if u:
            for k in d:
                d[k] = d[k] * u.cm**-3
        return d

    @density.setter
    def density(self, collider_density):

        collider_ids = {'H2': 0,
                        'PH2': 1,
                        'OH2': 2,
                        'E': 3,
                        'H': 4,
                        'HE': 5,
                        'H+': 6}

        collider_densities = defaultdict(lambda: 0)
        for k in collider_density:
            collider_densities[k.upper()] = collider_density[k]
            if k.upper() not in self._valid_colliders:
                raise ValueError('Collider %s is not one of the valid colliders: %s' %
                                 (k,self._valid_colliders))

        if 'OH2' in collider_densities:
            if not 'PH2' in collider_densities:
                raise ValueError("If o-H2 density is specified, p-H2 must also be.")
            # TODO: look up whether RADEX uses density[0] if density[1] and [2] are specified
            # (it looks like the answer is "no" based on a quick test)
            self.radex.cphys.density[0] = 0 # collider_densities['OH2'] + collider_densities['PH2']
            # PARA is [1], ORTHO is [2]
            # See lines 91, 92 of io.f
            self.radex.cphys.density[1] = collider_densities['PH2']
            self.radex.cphys.density[2] = collider_densities['OH2']
        elif 'H2' in collider_densities:
            warnings.warn("Using a default ortho-to-para ratio (which "
                          "will only affect species for which independent "
                          "ortho & para collision rates are given)")
            self.radex.cphys.density[0] = collider_densities['H2']

            T = self.temperature.value if hasattr(self.temperature,'value') else self.temperature
            if T > 0:
                opr = min(3.0,9.0*np.exp(-170.6/T))
            else:
                opr = 3.0
            fortho = opr/(1+opr)
            self.radex.cphys.density[1] = collider_densities['H2']*(1-fortho)
            self.radex.cphys.density[2] = collider_densities['H2']*(fortho)

        self.radex.cphys.density[3] = collider_densities['E']
        self.radex.cphys.density[4] = collider_densities['H']
        self.radex.cphys.density[5] = collider_densities['HE']
        self.radex.cphys.density[6] = collider_densities['H+']

        # skip H2 when computing by assuming OPR correctly distributes ortho & para
        # It's not obvious that RADEX does this correctly in readdata.f
        self.radex.cphys.totdens = self.radex.cphys.density[1:].sum()
    
    @property
    def total_density(self):
        if u:
            return self.radex.cphys.totdens * u.cm**-3
        else:
            return self.radex.cphys.totdens

    @property
    def mass_density(self):
        d = {'H2':self.radex.cphys.density[0]*0, # ignore H2, use o/p H2
             'pH2':self.radex.cphys.density[1]*2,
             'oH2':self.radex.cphys.density[2]*2,
             'e':self.radex.cphys.density[3]/1836.,
             'H':self.radex.cphys.density[4]*2,
             'He':self.radex.cphys.density[5]*4,
             'H+':self.radex.cphys.density[6]}
        if u:
            return np.sum(d.values())*constants.m_p
        else:
            return np.sum(d.values())


    @property
    def opr(self):
        return self.radex.cphys.density[1]/self.radex.cphys.density[2]

    @property
    def molpath(self):
        return "".join(self.radex.impex.molfile).strip()

    @molpath.setter
    def molpath(self, molfile):
        if "~" in molfile:
            molfile = os.path.expandpath(molfile)
        self.radex.impex.molfile[:] = ""
        utils.verify_collisionratefile(molfile)
        self.radex.impex.molfile[:len(molfile)] = molfile

    @property
    def outfile(self):
        return self.radex.impex.outfile

    @outfile.setter
    def outfile(self, outfile):
        self.radex.impex.outfile[:len(self.outfile)] = outfile

    @property
    def logfile(self):
        return self.radex.setup.logfile

    @logfile.setter
    def logfile(self, logfile):
        self.radex.setup.logfile[:len(self.logfile)] = logfile

    @property
    def datapath(self):
        return os.path.expanduser("".join(self.radex.setup.radat).strip())

    @datapath.setter
    def datapath(self, radat):
        # self.radex data path not needed if molecule given as full path
        self.radex.setup.radat[:] = ""
        self.radex.setup.radat[:len(radat)] = radat


    @property
    def escapeProbGeom(self):
        return self.radex.setup.method

    @escapeProbGeom.setter
    def escapeProbGeom(self, escapeProbGeom):
        mdict = {'lvg':2,'sphere':1,'slab':3}
        if escapeProbGeom not in mdict:
            raise ValueError("Invalid escapeProbGeom, must be one of "+",".join(mdict))
        self.radex.setup.method = mdict[escapeProbGeom]
        

    @property
    def level_population(self):
        return self.radex.collie.xpop

    @property
    def tex(self):
        if u:
            return self.radex.radi.tex * u.K
        else:
            return self.radex.radi.tex

    Tex = tex

    @property
    def tau(self):
        # taul(iline) = cddv*(xpop(n)*gstat(m)/gstat(n)-xpop(m))
        #$         /(fgaus*xt/aeinst(iline))
        return self.radex.radi.taul

    @property
    def frequency(self):
        if u:
            return self.radex.radi.spfreq * u.GHz
        else:
            return self.radex.radi.spfreq

    @property
    def temperature(self):
        if u:
            return self.radex.cphys.tkin*u.K
        else:
            return self.radex.cphys.tkin

    @temperature.setter
    def temperature(self, tkin):
        if hasattr(tkin,'to'):
            tkin = tkin.to(u.K).value
        if tkin <= 0 or tkin > 1e4:
            raise ValueError('Must have kinetic temperature > 0 and < 10^4 K')
        self.radex.cphys.tkin = tkin
        
        if not os.path.exists(self.molpath):
            raise IOError("File not found: %s" % self.molpath)
        # must re-read molecular file and re-interpolate to new temperature
        self.radex.readdata()

    @property
    def column(self):
        if u:
            return self.radex.cphys.cdmol * u.cm**-2
        else:
            return self.radex.cphys.cdmol

    @column.setter
    def column(self, col):
        if hasattr(col,'to'):
            col = col.to('cm**-2').value
        if col < 1e5 or col > 1e25:
            raise ValueError("Extremely low or extremely high column.")
        self.radex.cphys.cdmol = col

    @property
    def column_per_bin(self):
        return self.column / self.deltav

    @column_per_bin.setter
    def column_per_bin(self, cddv):
        if u:
            self.column = cddv * self.deltav.to(u.km/u.s).value
        else:
            self.column = cddv * self.deltav

    @property
    def abundance(self):
        #abund = self.column / (self.total_density * self.length)
        #if u:
        #    return abund.decompose().value
        #else:
        #    return abund
        return self._abundance

    @abundance.setter
    def abundance(self, abund):
        col = self.h2column * abund
        self._abundance = abund
        #col = abund * self.total_density * self.length
        if u:
            # need to divide the column per km/s by
            self.column = col.to(u.cm**-2).value
        else:
            self.column = col

    @property
    def h2column(self):
        return self.column / self.abundance
        #return self.total_density * self.length

    @h2column.setter
    def h2column(self, nh2, unit='cm**-2'):
        if u:
            if not hasattr(nh2,'to'):
                nh2 = nh2*u.Unit(unit)
            self.column = nh2.to(u.cm**-2) * self.abundance
        else:
            self.column = nh2 * self.abundance

    @property
    def deltav(self):
        return self._deltav

    @deltav.setter
    def deltav(self, dv):
        if u:
            self._deltav = dv * u.km/u.s
        else:
            self._deltav = dv

    def length(self):
        """ Hard-coded, assumed length-scale """
        if u:
            return 1*u.pc
        else:
            return 1

    #@property
    #def length(self):
    #    if self.escapeProbGeom == 'lvg':
    #        return self.column_per_bin / self.total_density
    #    else:
    #        return self.h2column / self.total_density

    @property
    def debug(self):
        return self.radex.dbg.debug

    @debug.setter
    def debug(self, debug):
        self.radex.dbg.debug = debug

    @property
    def tbg(self):
        if u:
            return self.radex.cphys.tbg * u.K
        else:
            return self.radex.cphys.tbg

    @tbg.setter
    def tbg(self, tbg):
        #print("Set TBG=%f" % tbg)
        self.radex.cphys.tbg = tbg
        self.radex.backrad()

    def run_radex(self, silent=True, reuse_last=False, reload_molfile=True):
        """
        Run the iterative matrix solution using a python loop

        Parameters
        ----------
        silent: bool
            Print a message when iteration is done?
        reuse_last: bool
            If this is True, the matrix iterator will start at iteration 1
            rather than iteration 0, and it will therefore repopulate the rate
            matrix based on the radiative background alone.  In principle,
            setting this to True should result in a significantly faster
            convergence; in practice, it does not.
        reload_molfile: bool
            Re-read the molecular line file?  This is needed if the collision
            rates are different and have not been updated by, e.g., changing
            the temperature (which automatically runs the `readdata` function)
        """

        if reload_molfile or self.radex.collie.ctot.sum()==0:
            self.radex.readdata()

        #self.radex.backrad()
            
        self._iter_counter = 1 if reuse_last else 0
        
        converged = np.array(False)

        last = self.level_population.copy()

        while not converged:
            if self._iter_counter >= self.maxiter:
                if not silent:
                    print("Did not converge in %i iterations, stopping." % self.maxiter)
                break

            self.radex.matrix(self._iter_counter, converged)
            if np.abs(last-self.level_population).sum() < 1e-16 and self._iter_counter>self.miniter:
                if not silent:
                    print("Stopped changing after %i iterations" % self._iter_counter)
                break
            last = self.level_population.copy()
            self._iter_counter += 1

        if converged and not silent:
            print("Successfully converged after %i iterations" % self._iter_counter)

        return self._iter_counter

    @property
    def quantum_number(self):
        return np.array([("".join(x)).strip() for x in grouper(self.radex.quant.qnum.T.ravel().tolist(),6)])

    @property
    def upperlevelnumber(self):
        return self.radex.imolec.iupp

    @property
    def lowerlevelnumber(self):
        return self.radex.imolec.ilow

    @property
    def upperlevelindex(self):
        return self.radex.imolec.iupp-1

    @property
    def upperlevelpop(self):
        return self.level_population[self.upperlevelindex]

    @property
    def lowerlevelindex(self):
        return self.radex.imolec.ilow-1

    @property
    def lowerlevelpop(self):
        return self.level_population[self.lowerlevelindex]

    @property
    def upperstateenergy(self):
        return self.radex.rmolec.eup

    @property
    def source_area(self):
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
        if u:
            #return (self.line_flux * beamsize)
            # because each line has a different frequency, have to loop it
            OK_freqs = self.frequency != 0
            return u.Quantity([x.to(u.K, u.brightness_temperature(beamsize, f)).value
                               for x,f in zip(self.line_flux_density[OK_freqs],self.frequency[OK_freqs])
                               ],
                              unit=u.K)
        else:
            raise NotImplementedError("Astropy's units are required for this conversion.")

    @property
    def inds_frequencies_included(self):
        """
        The indices of the line frequencies fitted by RADEX
        (RADEX can hold up to 99999 frequencies, but usually uses ~100)
        """
        OK_freqs = self.frequency != 0
        return np.where(OK_freqs)[0]

    @property
    def source_line_brightness_temperature(self):
        """
        The surface brightness of the source assuming it is observed with a
        beam matched to its size and it has ff=1
        """
        if u:
            #return (self.line_flux * beamsize)
            # because each line has a different frequency, have to loop it
            OK_freqs = self.inds_frequencies_included
            return u.Quantity([(x*u.sr).to(u.K, u.brightness_temperature(1*u.sr, f)).value
                               for x,f in zip(self.source_line_surfbrightness[OK_freqs],
                                              self.frequency[OK_freqs])
                               ],
                              unit=u.K)
        else:
            raise NotImplementedError("Astropy's units are required for this conversion.")

    @property
    def T_B(self):
        return self.source_line_brightness_temperature

    @property
    def background_brightness(self):
        if u:
            return self.radex.radi.backi * u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1 * u.sr**-1
        else:
            return self.radex.radi.backi

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

        if u:
            thc = (2 * constants.h * constants.c).cgs / u.sr
            fk = (constants.h * constants.c / constants.k_B).cgs
        else:
            thc = 3.9728913665386055e-16
            fk = 1.4387769599838154

        ftau = np.exp(-self.tau)
        xt = self._xt
        xnu = self._xnu
        bnutex = thc*xt/(np.exp(fk*xnu/self.tex)-1.0)
        toti = self.background_brightness*ftau+bnutex*(1.0-ftau)

        return toti

    @property
    def source_brightness_beta(self):
        if u:
            thc = (2 * constants.h * constants.c).cgs / u.sr
            fk = (constants.h * constants.c / constants.k_B).cgs
        else:
            thc = 3.9728913665386055e-16
            fk = 1.4387769599838154
        xt = self._xt
        xnu = self._xnu
        bnutex = thc*xt/(np.exp(fk*xnu/self.tex)-1.0)
        toti = self.background_brightness*ftau+bnutex*(1-self.beta)
        return toti

    @property
    def beta(self):
        # this will probably be faster if vectorized (translated completely
        # from fortran to python)
        return np.array([self.radex.escprob(t) for t in self.tau])

    @property
    def _xnu(self):
        """
        Line frequency in inverse cm
        """
        if u:
            return self.radex.radi.xnu * u.cm**-1
        else:
            return self.radex.radi.xnu

    @property
    def _xt(self):
        # xt = xnu**3 # cm^-1 -> cm^-3
        return self._xnu**3

    @property
    def _cddv(self):
        return self.column / self.deltav

    @property
    def statistical_weight(self):
        return self.radex.rmolec.gstat

    def get_table(self):
        T = astropy.table.Table()
        mask = self.frequency.value != 0
        T.add_column(astropy.table.Column(name='Tex',data=self.tex[mask], unit='K'))
        T.add_column(astropy.table.Column(name='tau',data=self.tau[mask], unit=''))
        T.add_column(astropy.table.Column(name='frequency',data=self.frequency[mask], unit='GHz'))
        T.add_column(astropy.table.Column(name='upperstateenergy',data=self.upperstateenergy[mask], unit='K'))
        T.add_column(astropy.table.Column(name='upperlevel',data=self.quantum_number[self.upperlevelindex[mask]], unit=''))
        T.add_column(astropy.table.Column(name='lowerlevel',data=self.quantum_number[self.lowerlevelindex[mask]], unit=''))
        T.add_column(astropy.table.Column(name='upperlevelpop',data=self.level_population[self.upperlevelindex[mask]], unit=''))
        T.add_column(astropy.table.Column(name='lowerlevelpop',data=self.level_population[self.lowerlevelindex[mask]], unit=''))
        T.add_column(astropy.table.Column(name='brightness',data=self.source_line_surfbrightness[mask]))
        T.add_column(astropy.table.Column(name='T_B',data=self.T_B)) # T_B is pre-masked
        if self.source_area:
            T.add_column(astropy.table.Column(name='flux',data=self.line_flux_density[mask]))

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


def density_distribution(meandens, **kwargs):
    """
    Compute the LVG model for a single zone with an assumed density
    *distribution* but other properties fixed.
    """
    pass

def grid():
    pass
