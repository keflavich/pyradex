from __future__ import print_function
import subprocess
import tempfile
import astropy.io.ascii 
import astropy.table
import StringIO
import re
import numpy as np

from read_radex import read_radex

__all__ = ['radex','write_input','parse_outfile', 'call_radex']

def radex(executable='radex', flow=100, fhigh=130, collider_densities={'H2':1},
        debug=False, delete_tempfile=True, **kwargs):
    """
    Get the radex results for a set of input parameters


    Parameters
    ----------
    executable : str
        Full path to the RADEX executable
    flow : float
        Lowest frequency line to store
    fhigh : float
        Highest frequency line to store

    See write_input for additional parameters
    """

    infile,outfile = write_input(flow=flow, fhigh=fhigh,
            delete_tempfile=delete_tempfile,
            collider_densities=collider_densities, **kwargs)

    logfile = call_radex(executable, infile.name, debug=debug,
            delete_tempfile=delete_tempfile)

    data = parse_outfile(outfile.name)#, len(collider_densities), flow, fhigh)

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

def write_input(tkin=10, column_density=1e12, collider_densities={'H2':1},
        bw=0.01, tbg=2.73, molecule='co', velocity_gradient=1.0, flow=1,
        fhigh=10, delete_tempfile=True):
    """
    Write radex.inp file parameters

    Parameters
    ----------
    tkin : float
        Kinetic temperature (K)
    collider_densities : dict
        Collider names and their number densities
    column_density : float
        column density of the molecule
    molecule : str
        Name of the molecule (specifically, the prefix for the file name, e.g.
        for "co.dat", molecule='co').  Case sensitive!
    tbg : float
        Temperature of the background radiation (e.g. CMB)
    velocity_gradient : float
        Velocity gradient per pc in km/s
    """
    infile = tempfile.NamedTemporaryFile(mode='w', delete=delete_tempfile)
    outfile = tempfile.NamedTemporaryFile(mode='w', delete=delete_tempfile)
    infile.write(molecule+'.dat\n')
    infile.write(outfile.name+'\n')
    infile.write(str(flow)+' '+str(fhigh)+'\n')
    infile.write(str(tkin)+'\n')
    infile.write('%s\n' % len(collider_densities))
    for name,dens in collider_densities.iteritems():
        infile.write('%s\n' % name)
        infile.write(str(dens)+'\n')
    infile.write(str(tbg)+'\n')
    infile.write(str(column_density)+'\n')
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

    result = subprocess.call(cmd,
        shell=True)
    if result != 0:
        print("RADEX returned error code %i" % result)
        with open(logfile.name,'r') as f:
            print(f.read())

    return logfile


import astropy.units as u
header_names = ['Row#','Line#','E_UP','FREQ', 'WAVE', 'T_EX', 'TAU', 'T_R', 'POP_UP', 'POP_LOW', 'FLUX_Kkms',   'FLUX_Inu']
header_units = [None,       None, u.K,   u.GHz,  u.um,   u.K,    None,  u.K,   None,     None,     u.K*u.km/u.s, u.erg/u.cm**2/u.s]
dtypes       = [int,   int,    float, float,  float,  float,  float, float,  float,    float,     float,         float]

def parse_outfile(filename):
    with open(filename,'r') as f:
        lines = [L.replace("--","  ") for L in f 
                if (L[0] != '*' 
                    and 'iterat' not in L 
                    and 'GHz' not in L 
                    and 'TAU' not in L)]
    data_list = [L.split() for L in lines]
    data_in_columns = map(list,zip(*data_list))
    columns = [astropy.table.Column(data=C, name=name, units=units, dtype=dtype) 
            for C,name,units,dtype in zip(data_in_columns, header_names, header_units, dtypes)]
    data = astropy.table.Table(columns)
    return data
