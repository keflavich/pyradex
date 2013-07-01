import subprocess
import tempfile
from __future__ import print_function

def radex(executable='radex', flow=1, fhigh=10, **kwargs):
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

    infile,outfile = write_input(flow=flow, fhigh=fhigh, **kwargs)

    logfile = call_radex(executable, outfile.name)

    parse_outfile(outfile.name)

def write_input(tkin=10, column_density=1e12, collider_densities={'H2':1},
        bw=0.01, tbg=2.73, molecule='co', velocity_gradient=1.0, flow=1,
        fhigh=10):
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
    infile = tempfile.NamedTemporaryFile(mode='w', delete=True)
    outfile = tempfile.NamedTemporaryFile(mode='w', delete=True)
    infile.write(molecule+'.dat\n')
    infile.write(outfile.name+'\n')
    infile.write(str(flow)+' '+str(fhigh)+'\n')
    infile.write(str(tkin)+'\n')
    infile.write('%s\n' % len(colliders))
    for name,dens in colliders.iteritems():
        infile.write('%s\n' % name)
        infile.write(str(dens)+'\n')
    infile.write(str(tbg)+'\n')
    infile.write(str(column_density)+'\n')
    infile.write(str(velocity_gradient)+'\n')
    # end the input file
    infile.write('0\n')
    infile.flush()
    return infile,outfile

def call_radex(executable, inpfilename):

    logfile = tempfile.NamedTemporaryFile(mode='w', delete=True)
    result = subprocess.call('{radex} < {inpfile} > {logfile}'.format(
        radex=executable,
        inpfile=inpfilename,
        logfile=logfile.name),
        shell=True)
    if result != 0:
        print("RADEX returned error code %i" % result)

    return logfile
