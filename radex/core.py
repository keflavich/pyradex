import subprocess
import tempfile

def radex(inputs, executable='radex'):

    infile,outfile = write_input(*args)

    call_radex(executable, outfile)

def write_input(tkin=10, column_density=1e12, collider_densities={'H2':1}, 
        bw=0.01, tbg=2.73, molecule='co', velocity_gradient=1.0):
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
    infile.write(str(flow*(1-bw))+' '+str(fupp/(1-bw))+'\n')
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

def call_radex(executable):

    os.system('%s < radex.inp > /dev/null' % executable)
