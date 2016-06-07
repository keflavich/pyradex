from __future__ import print_function
import subprocess
import tempfile
import numpy as np
import warnings
import astropy.units as u
_quantity = u.Quantity
from collections import defaultdict
import os
import sys

from . import utils
from . import synthspec
from .utils import QuantityOff,ImmutableDict,unitless,grouper
from .base_class import RadiativeTransferApproximator

from astropy import units as u
from astropy import constants
from astropy import log
import astropy.table

PYVERSION = 3 if sys.version_info >= (3,0) else 2

__all__ = ['pyradex', 'write_input', 'parse_outfile', 'call_radex', 'Radex',
           'density_distribution']

        

def pyradex(executable='radex', minfreq=100, maxfreq=130,
            collider_densities={'H2':1}, debug=False, delete_tempfile=True,
            return_dict=False, **kwargs):
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
    warnings.warn("pyradex is deprecated: Use pyradex.Radex instead if you can.")

    infile,outfile = write_input(minfreq=minfreq, maxfreq=maxfreq,
            delete_tempfile=delete_tempfile,
            collider_densities=collider_densities, **kwargs)

    logfile = call_radex(executable, infile.name, debug=debug,
                         delete_tempfile=delete_tempfile)

    check_logfile(logfile.name)

    data = parse_outfile(outfile.name, return_dict=return_dict)

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
        minfreq = unitless(minfreq.to('GHz',u.spectral()))
    if hasattr(maxfreq, 'unit'):
        maxfreq = unitless(maxfreq.to('GHz',u.spectral()))

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

def parse_outfile(filename, return_dict=False):
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
        niter = [L.split(" ")[3]
                for L in alllines
                if 'iterat' in L]
    data_list = [[x if '*' not in x else '-999' for x in L.split()] for L in lines]
    if len(data_list) == 0:
        raise ValueError("No lines included?")
    data_in_columns = map(list,zip(*data_list))
    if return_dict:
        data = {name: C for C,name in zip(data_in_columns, header_names)}
        data['niter']=niter
        return data
    columns = [astropy.table.Column(data=C, name=name.lower(), unit=unit, dtype=dtype) 
            for C,name,unit,dtype in zip(data_in_columns, header_names, header_units, dtypes)]
    data = astropy.table.Table(columns, meta=header)
    return data

class Radex(RadiativeTransferApproximator):

    def __call__(self, return_table=True, **kwargs):
        # reset the parameters appropriately
        self.set_params(**kwargs)
        # No need to re-validate: it is already done when self.temperature is
        # set in __init__
        niter = self.run_radex(reload_molfile=False, validate_colliders=False)

        if return_table:
            return self.get_table()
        else:
            return niter

    def __init__(self,
                 collider_densities=None,
                 density=None,
                 total_density=None,
                 temperature=None,
                 species='co',
                 column=None,
                 column_per_bin=None,
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
            Dictionary giving the volume densities of the collider(s) in units
            of cm^-3.  Valid entries are h2,oh2,ph2,e,He,H,H+.  The keys are
            case-insensitive.
        density: float
        total_density: float
            (optional) Alternative to ``collider_densities``: can specify a
            single number indicating the total density of H2.  This should
            not be used when electrons or H atoms are the intended collider.
            These keywords are synonymous and therefore only one can be used.
        temperature: float
            Local gas temperature in K
        species: str
            A string specifying a valid chemical species.  This is used to look
            up the specified molecule
        column: float
        column_per_bin : float
            The column density of the molecule of interest per bin, where
            a bin is (deltav km/s * 1 pc). These keywords are synonymous and
            therefore only one can be specified.
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
                raise ValueError("Data path %s was not successfully stored;"
                                 " instead %s was." % (datapath,self.datapath))
        self.species = species
        if self.molpath == b'':
            raise ValueError("Must set a species name.")
        if not os.path.exists(self.molpath):
            raise ValueError("Must specify a valid path to a molecular data file "
                             "else RADEX will crash."
                             "  Current path is {0}".format(self.molpath))

        if sum(x is not None for x in (collider_densities,density,total_density)) > 1:
            raise ValueError("Can only specify one of density, total_density,"
                             " and collider_densities")

        if sum(x is not None for x in (column,column_per_bin)) > 1:
            raise ValueError("Can only specify one of column, column_per_bin.")

        n_specifications = sum(x is not None for x in (column, column_per_bin,
                                                       collider_densities,
                                                       density, total_density,
                                                       abundance)) 
        if (n_specifications > 2):
            raise ValueError("Can only specify two of column, density, and abundance.")
        if (n_specifications < 2):
            raise ValueError("Must specify two of column, density, and abundance.")

        self._locked_parameter = 'density'
        self._is_locked = True

        # This MUST happen before density is set, otherwise OPR will be 
        # incorrectly set.
        self.radex.cphys.tkin = unitless(temperature)

        if collider_densities:
            self.density = collider_densities
            self._is_locked = False
            if total_density:
                log.warn("`total_density` was specified, but `collider_densities` "
                         "was used instead.  Set `collider_densities=None` if you "
                         "want to use `total_density`.")
        elif total_density:
            self.density = total_density
            self._is_locked = False
        elif density:
            self.density = density
            self._is_locked = False
        else:
            self._locked_parameter = 'column'
            self._is_locked = True

        self.outfile = outfile
        self.logfile = logfile
        self.escapeProbGeom = escapeProbGeom

        self.deltav = deltav
        self._set_parameters()

        if column_per_bin is not None:
            self.column_per_bin = column_per_bin
        elif column is not None:
            self.column_per_bin = column
        else:
            self._locked_parameter = 'density'

        self._is_locked = False

        if abundance:
            self.abundance = abundance

        # This has to happen here, because the colliders are read in at
        # this point and rates interpolated
        self.temperature = temperature
        self.tbg = tbackground

        self.debug = debug

        self.source_area = source_area

    _u_brightness = (u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1 * u.sr**-1)
    _u_sc = u.cm**-2
    _u_cc = u.cm**-3

    _u_gradient = u.cm**-2 / (u.km/u.s) / u.pc
    _u_kms = u.km/u.s
    _u_cms = u.cm/u.s

    def set_params(self, density=None, collider_densities=None,
                   column=None, column_per_bin=None, temperature=None,
                   abundance=None, species=None, deltav=None, tbg=None,
                   escapeProbGeom=None):

        if species is not None:
            self.species = species

        if deltav is not None:
            self.deltav = deltav

        # This MUST happen before density is set, otherwise OPR will be 
        # incorrectly set.
        if temperature is not None:
            self.radex.cphys.tkin = unitless(temperature)

        if collider_densities is not None:
            self.density = collider_densities
        elif density is not None:
            if collider_densities is not None:
                raise ValueError('Can specify only one of density,'
                                 ' collider_densities')
            self.density = density

        if column is not None:
            self.column = column
        elif column_per_bin is not None:
            if column is not None:
                raise ValueError("Can specify only one of column,"
                                 "column per bin")
            self.column_per_bin = column_per_bin

        if temperature is not None:
            self.temperature = temperature

        if abundance is not None:
            self.abundance = abundance

        if tbg is not None:
            self.tbg = tbg

        if escapeProbGeom is not None:
            self.escapeProbGeom = escapeProbGeom

    @property
    def locked_parameter(self):
        return self._locked_parameter

    def _lock_param(self, parname):
        self._locked_parameter = parname

    def _set_parameters(self):

        #self.radex.cphys.cdmol = self.column
        #self.radex.cphys.tkin = self.temperature
        if hasattr(self.deltav, 'to'):
            self.radex.cphys.deltav = unitless(self.deltav.to(self._u_cms))
        else:
            self.radex.cphys.deltav = self.deltav * (self._u_cms.to(self._u_kms))

        # these parameters are only used for outputs and therefore can be ignored
        self.radex.freq.fmin = 0
        self.radex.freq.fmax = 1e10

        if not hasattr(self, 'miniter'):
            self.miniter = 10
        if not hasattr(self, 'maxiter'):
            self.maxiter = 200

    _all_valid_colliders = {'H2':'H2',
                            'PH2':'pH2',
                            'OH2':'oH2',
                            'E':'e',
                            'H':'H',
                            'HE':'He',
                            'H+':'H+'}

    @property
    def density(self):

        d = {'H2':self.radex.cphys.density[0],
             'pH2':self.radex.cphys.density[1],
             'oH2':self.radex.cphys.density[2],
             'e':self.radex.cphys.density[3],
             'H':self.radex.cphys.density[4],
             'He':self.radex.cphys.density[5],
             'H+':self.radex.cphys.density[6]}

        for k in d:
            d[k] = u.Quantity(d[k], self._u_cc)
        
        return ImmutableDict(d)

    @density.setter
    def density(self, collider_density):

        collider_ids = {'H2': 0,
                        'PH2': 1,
                        'OH2': 2,
                        'E': 3,
                        'H': 4,
                        'HE': 5,
                        'H+': 6}

        self._use_thermal_opr = False

        if isinstance(collider_density, (float,int,_quantity,np.ndarray)):
            log.warn("Assuming the density is n(H_2).")
            collider_density = {'H2': collider_density}

        collider_densities = defaultdict(lambda: 0)
        for k in collider_density:
            collider_densities[k.upper()] = unitless(u.Quantity(collider_density[k], self._u_cc))
            if k.upper() not in self._all_valid_colliders:
                raise ValueError('Collider %s is not one of the valid colliders: %s' %
                                 (k,self._all_valid_colliders))

        if (('OH2' in collider_densities and collider_densities['OH2'] !=0) or
            ('PH2' in collider_densities and collider_densities['PH2'] !=0)):

            # this is simply not true: NH3 has just ph2 as a collider
            #if not 'PH2' in collider_densities or not 'OH2' in collider_densities:
            #    raise ValueError("If o-H2 density is specified, p-H2 must also be.")
            # TODO: look up whether RADEX uses density[0] if density[1] and [2] are specified
            # (it looks like the answer is "no" based on a quick test)
            #self.radex.cphys.density[0] = 0 # collider_densities['OH2'] + collider_densities['PH2']
            # PARA is [1], ORTHO is [2]
            # See lines 91, 92 of io.f
            if 'PH2' in collider_densities:
                self.radex.cphys.density[1] = collider_densities['PH2']
            if 'OH2' in collider_densities:
                self.radex.cphys.density[2] = collider_densities['OH2']
            self._use_thermal_opr = False
        elif 'H2' in collider_densities:
            warnings.warn("Using a default ortho-to-para ratio (which "
                          "will only affect species for which independent "
                          "ortho & para collision rates are given)")
            self._use_thermal_opr = True
            #self.radex.cphys.density[0] = collider_densities['H2']

            T = unitless(self.temperature)
            if T > 0:
                # From Faure, private communication
                opr = min(3.0,9.0*np.exp(-170.6/T))
            else:
                opr = 3.0
            fortho = opr/(1+opr)
            log.debug("Set OPR to {0} and fortho to {1}".format(opr,fortho))
            self.radex.cphys.density[1] = collider_densities['H2']*(1-fortho)
            self.radex.cphys.density[2] = collider_densities['H2']*(fortho)

        # RADEX relies on n(H2) = n(oH2) + n(pH2)
        # We have set n(oH2) and n(pH2) above
        vc = [x.lower() for x in self.valid_colliders]
        if 'h2' in vc:
            self.radex.cphys.density[0] = self.radex.cphys.density[1:3].sum()
            self.radex.cphys.density[1] = 0
            self.radex.cphys.density[2] = 0
        elif 'oh2' in vc or 'ph2' in vc:
            self.radex.cphys.density[0] = 0

        self.radex.cphys.density[3] = collider_densities['E']
        self.radex.cphys.density[4] = collider_densities['H']
        self.radex.cphys.density[5] = collider_densities['HE']
        self.radex.cphys.density[6] = collider_densities['H+']

        # skip H2 when computing by assuming OPR correctly distributes ortho & para
        # It's not obvious that RADEX does this correctly in readdata.f
        self.radex.cphys.totdens = self.radex.cphys.density.sum()

        # Unfortunately,
        # must re-read molecular file and re-interpolate to new density
        self._validate_colliders()
        self.radex.readdata()

        if not self._is_locked:
            self._is_locked = True
            if self.locked_parameter == 'column':
                self.abundance = self.column_per_bin /(self.total_density*self.length)
            elif self.locked_parameter == 'abundance':
                self.column_per_bin = self.total_density * self.length * self.abundance
            self._lock_param('density')
            self._is_locked = False


    @property
    def valid_colliders(self):
        return self._valid_colliders
    
    @property
    def total_density(self):
        """
        The total density *by number of particles* 
        The *mass density* can be dramatically different!
        """
        return u.Quantity(self.radex.cphys.totdens, self._u_cc)


    @property
    def opr(self):
        return self.radex.cphys.density[1]/self.radex.cphys.density[2]

    @property
    def molpath(self):
        return b"".join(self.radex.impex.molfile).strip()

    @molpath.setter
    def molpath(self, molfile):
        if "~" in molfile:
            molfile = os.path.expanduser(molfile)
        if PYVERSION == 3:
            self.radex.impex.molfile[:] = np.bytes_([""]*len(self.radex.impex.molfile))
        else:
            self.radex.impex.molfile[:] = ""
        utils.verify_collisionratefile(molfile)
        self.radex.impex.molfile[:len(molfile)] = molfile

    @property
    def outfile(self):
        return self.radex.impex.outfile

    @outfile.setter
    def outfile(self, outfile):
        if PYVERSION == 3:
            self.radex.impex.outfile[:] = np.bytes_([""]*len(self.radex.impex.outfile))
        else:
            self.radex.impex.outfile[:] = ""
        self.radex.impex.outfile[:len(outfile)] = outfile

    @property
    def logfile(self):
        return self.radex.setup.logfile

    @logfile.setter
    def logfile(self, logfile):
        if PYVERSION == 3:
            self.radex.setup.logfile[:] = np.bytes_([""]*len(self.radex.setup.logfile))
        else:
            self.radex.setup.logfile[:] = ""
        self.radex.setup.logfile[:len(logfile)] = logfile

    @property
    def datapath(self):
        return os.path.expanduser(b"".join(self.radex.setup.radat).strip()).decode('utf-8')

    @datapath.setter
    def datapath(self, radat):
        # self.radex data path not needed if molecule given as full path
        if PYVERSION == 3:
            self.radex.setup.radat[:] = np.bytes_([""] * len(self.radex.setup.radat))
        else:
            self.radex.setup.radat[:] = ""
        # there is dangerous magic here: radat needs to be interpreted as an array,
        # but you can't make it an array of characters easily...
        self.radex.setup.radat[:len(radat)] = radat


    @property
    def escapeProbGeom(self):
        mdict = {2:'lvg',1:'sphere',3:'slab'}
        return mdict[int(self.radex.setup.method)]

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
        return u.Quantity(self.radex.radi.tex[self._mask], u.K)

    Tex = tex

    @property
    def tau(self):
        # taul(iline) = cddv*(xpop(n)*gstat(m)/gstat(n)-xpop(m))
        #$         /(fgaus*xt/aeinst(iline))
        return self.radex.radi.taul[self._mask]

    @property
    def frequency(self):
        return u.Quantity(self.radex.radi.spfreq[self._mask], u.GHz)

    @property
    def temperature(self):
        return u.Quantity(self.radex.cphys.tkin, u.K)

    @temperature.setter
    def temperature(self, tkin):
        if hasattr(tkin,'to'):
            tkin = unitless(u.Quantity(tkin, u.K))
        elif tkin is None:
            raise TypeError("Must specify tkin")

        if tkin <= 0 or tkin > 1e4:
            raise ValueError('Must have kinetic temperature > 0 and < 10^4 K')
        self.radex.cphys.tkin = tkin
        
        if not os.path.exists(self.molpath):
            raise IOError("File not found: %s" % self.molpath)
        # must re-read molecular file and re-interpolate to new temperature
        self._validate_colliders()
        #log.info("before DENS:"+str(self.radex.cphys.density))
        #log.info("before TOTDENS:"+str(self.radex.cphys.totdens))
        self.radex.readdata()
        #log.info("after DENS:"+str(self.radex.cphys.density))
        #log.info("after TOTDENS:"+str(self.radex.cphys.totdens))

        if self._use_thermal_opr:
            # Reset the density to a thermal value
            lp = self._locked_parameter
            self.density = (unitless(self.density['H2']) or
                            unitless(self.density['oH2']+self.density['pH2']))
            self._locked_parameter = lp

    @property
    def column(self):
        return self.column_per_bin

    @column.setter
    def column(self, value):
        self.column_per_bin = value

    @property
    def column_per_bin(self):
        return u.Quantity(self.radex.cphys.cdmol, self._u_sc)

    @column_per_bin.setter
    def column_per_bin(self, col):
        if hasattr(col,'to'):
            col = unitless(u.Quantity(col, self._u_sc))
        if col < 1e5 or col > 1e25:
            raise ValueError("Extremely low or extremely high column.")
        self.radex.cphys.cdmol = col

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
    def debug(self):
        return self.radex.dbg.debug

    @debug.setter
    def debug(self, debug):
        self.radex.dbg.debug = debug

    @property
    def tbg(self):
        return u.Quantity(self.radex.cphys.tbg, u.K)

    @tbg.setter
    def tbg(self, tbg):
        #print("Set TBG=%f" % tbg)
        if hasattr(tbg, 'value'):
            tbg = unitless(u.Quantity(tbg, u.K))
        self.radex.cphys.tbg = tbg
        self.radex.backrad()

    def run_radex(self, silent=True, reuse_last=False, reload_molfile=True,
                  abs_convergence_threshold=1e-16, rel_convergence_threshold=1e-8,
                  validate_colliders=True):
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
        validate_colliders: bool
            Validate the colliders before running the code.  This should always
            be done unless running in a grid, in which case it can cause a
            slowdown (~30%).
        """
        if validate_colliders:
            # 100 loops, best of 3: 7.48 ms per loop
            self._validate_colliders()

        if reload_molfile or self.radex.collie.ctot.sum()==0:
            # 100 loops, best of 3: 15.3 ms per loop
            self.radex.readdata()

        #self.radex.backrad()
        
        # Given the properties of *this* class, set the appropriate RADEX
        # fortran function values
        # 10000 loops, best of 3: 74 micros per loop
        self._set_parameters()
            
        self._iter_counter = 1 if reuse_last else 0
        
        converged = np.array(False)

        # 1000000 loops, best of 3: 1.79 micros per loop
        last = self.level_population.copy()

        while not converged:
            if self._iter_counter >= self.maxiter:
                if not silent:
                    print("Did not converge in %i iterations, stopping." % self.maxiter)
                break

            # 10000 loops, best of 3: 30.8 micros per loop
            self.radex.matrix(self._iter_counter, converged)
            level_diff = np.abs(last-self.level_population)
            frac_level_diff = level_diff/self.level_population
            if (((level_diff.sum() < abs_convergence_threshold) or
                 (frac_level_diff.sum() < rel_convergence_threshold)) and
                self._iter_counter>self.miniter):
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
        return np.array([(b"".join(x)).strip() for x in
                         grouper(self.radex.quant.qnum.T.ravel().tolist(),6)])

    @property
    def upperlevelnumber(self):
        # wrong return self.radex.imolec.iupp[self._mask]
        return self.quantum_number[self.upperlevelindex]

    @property
    def lowerlevelnumber(self):
        # wrong return self.radex.imolec.ilow[self._mask]
        return self.quantum_number[self.lowerlevelindex]

    @property
    def upperlevelindex(self):
        return self.radex.imolec.iupp[self._mask]-1

    @property
    def upperlevelpop(self):
        return self.level_population[self.upperlevelindex]

    @property
    def lowerlevelindex(self):
        return self.radex.imolec.ilow[self._mask]-1

    @property
    def lowerlevelpop(self):
        return self.level_population[self.lowerlevelindex]

    @property
    def upperstateenergy(self):
        return self.radex.rmolec.eup[self._mask]

    @property
    def inds_frequencies_included(self):
        """
        The indices of the line frequencies fitted by RADEX
        (RADEX can hold up to 99999 frequencies, but usually uses ~100)
        """
        return np.where(self._mask)[0]

    @property
    def background_brightness(self):
        return u.Quantity(self.radex.radi.backi[self._mask], self._u_brightness)

    _thc = (2 * constants.h * constants.c).cgs / u.sr
    _fk = (constants.h * constants.c / constants.k_B).cgs
    _thc_value = _thc.value
    _fk_value = _fk.value

    @property
    def source_brightness(self):
        """
        RADEX compat?  (check)
        """

        fk = self._fk_value
        thc = self._thc_value

        with QuantityOff():
            ftau = np.exp(-self.tau)
            xt = self._xt
            xnu = self._xnu
            earg = fk*xnu/self.tex
            bnutex = thc*xt/(np.exp(earg)-1.0)
            toti = self.background_brightness*ftau+bnutex*(1.0-ftau)

        return u.Quantity(toti, self._u_brightness)

    @property
    def source_brightness_beta(self):
        fk = self._fk_value
        thc = self._thc_value

        with QuantityOff():
            ftau = np.exp(-self.tau)
            xt = self._xt
            xnu = self._xnu
            earg = fk*xnu/self.tex
            bnutex = thc*xt/(np.exp(earg)-1.0)
            toti = self.background_brightness*ftau+bnutex*(1-self.beta)

        return u.Quantity(toti, self._u_brightness)

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
        return u.Quantity(self.radex.radi.xnu[self._mask], u.cm**-1)

    @property
    def _xt(self):
        # xt = xnu**3 # cm^-1 -> cm^-3
        return self._xnu**3

    @property
    def _cddv(self):
        return self.column / self.deltav

    @property
    def _statistical_weight(self):
        return self.radex.rmolec.gstat

    @property
    def upperlevel_statisticalweight(self):
        return self._statistical_weight[self.upperlevelindex]

    @property
    def lowerlevel_statisticalweight(self):
        return self._statistical_weight[self.lowerlevelindex]

    @property
    def _mask(self):
        return self.radex.radi.spfreq != 0

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

    def partition_function(self, temperature=None):
        """
        Equation 46 of Mangum & Shirley 2015:

            Q = Sum( g_i exp(-E_i / kT) )
        """
        warnings.warn("The partition function may be very inaccurate using "
                      "LAMDA files because they include a small fraction of"
                      " the total available states.")
        gi = self.upperlevel_statisticalweight
        Ei = u.Quantity(self.upperstateenergy, unit=u.K)
        if temperature is None:
            temperature = self.temperature
        if not hasattr(temperature, 'unit'):
            temperature = u.Quantity(temperature, unit=u.K)
        return (gi*np.exp(-Ei/(temperature))).sum()


def density_distribution(densarr, distr, moleculecolumn, tauthresh=0.8,
                         opr=None, line_ids=[], mincol=None, Radex=Radex,
                         **kwargs):
    """
    Compute the LVG model for a single zone with an assumed density
    *distribution* but other properties fixed.

    Parameters
    ----------
    dendarr : array
        Array of densities corresponding to the distribution function
    distr : array
        The density distribution corresponding to the density array
    moleculecolumn : quantity
        The total column density of the molecule in question.  It will be
        redistributed across the appropriate densities.  Units: cm^-2
        [this is wrong - each density will assume a too-low optical depth]
    """
    try:
        np.testing.assert_almost_equal(distr.sum(), 1)
    except AssertionError:
        raise ValueError("The distribution must be normalized.")

    if not line_ids:
        raise ValueError("Specify at least one line ID")

    meandens = (densarr*distr).mean()
    if opr is None:
        collider_densities = {'H2': meandens}
    else:
        fortho = opr/(1+opr)
        collider_densities = {'oH2':meandens*fortho,'pH2':meandens*(1-fortho)}


    # Test whether the multi-slab model is reasonable by checking:
    # if the column was all at the mean density, would any lines be
    # optically thick?
    R = Radex(collider_densities=collider_densities, column=moleculecolumn, **kwargs)
    R.run_radex()
    if np.any(R.tau > tauthresh):
        warnings.warn(("At least one line optical depth is >{tauthresh}.  "
                       "Smoothing may be invalid.").format(tauthresh=tauthresh))

    # set the optical depth from the *mean* density assuming the *total* column
    tau = R.tau
    print("Mean density: {0}  Optical Depth: {1}".format(meandens, tau[line_ids]))

    _thc = (2 * constants.h * constants.c).cgs / u.sr
    _fk = (constants.h * constants.c / constants.k_B).cgs
    _thc_value = _thc.value
    _fk_value = _fk.value
    _u_brightness = (u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1 * u.sr**-1)

    xnu = R.frequency.to(u.cm**-1, u.spectral()).value

    linestrengths = []
    texs = []
    for dens,prob in zip(densarr,distr):
        if opr is None:
            collider_densities = {'H2':dens}
        else:
            collider_densities = {'oH2':dens*fortho,'pH2':dens*(1-fortho)}
        R.density = collider_densities
        try:
            R.column = moleculecolumn * prob
            if mincol is not None and R.column < mincol:
                R.column = mincol
            R.run_radex()
        except ValueError as ex:
            if ex.args[0] == "Extremely low or extremely high column.":
                if R.column > u.Quantity(1e20, u.cm**-2):
                    raise ex
                else:
                    texs.append(np.zeros_like(line_ids)+2.73)
                    linestrengths.append(np.zeros_like(line_ids))
                    continue
            else:
                raise ex
        
        if hasattr(R, 'radex'):
            R.radex.radi.taul[:len(tau)] = tau
        elif hasattr(R, '_data_dict'):
            R._data_dict['tau'] = tau

        fk = _fk_value
        thc = _thc_value

        with QuantityOff():
            ftau = np.exp(-tau)
            xt = xnu**3
            earg = fk*xnu/R.tex
            bnutex = thc*xt/(np.exp(earg)-1.0)
            toti_nounit = R.background_brightness*ftau+bnutex*(1.0-ftau)

        toti = u.Quantity(toti_nounit, _u_brightness)
        totK = ((toti*u.sr).to(u.K, u.brightness_temperature(1*u.sr,
                                                             R.frequency)))

        linestrengths.append(totK[line_ids])
        texs.append(R.tex[line_ids])

    linestrengths = np.array(linestrengths)
    texs = np.array(texs)

    return R, linestrengths, linestrengths.sum(axis=0), texs, tau[line_ids]


def grid():
    pass
