from astropy import units as u
import os
import errno

def mkdir_p(path):
    """ mkdir -p equivalent [used by get_datafile]"""
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def united(qty, unit):
    if isinstance(qty,u.Quantity):
        return qty.to(unit)
    else:
        return qty*u.Unit(unit)

def uvalue(qty, unit):
    return united(qty, unit).value

def get_datafile(species, savedir='./'):
    """
    Load a molecular data file and save it into the specified directory
    """
    from astroquery import lamda

    datapath = os.path.join(savedir,species)

    species,suffix = os.path.splitext(species)
    if suffix == "":
        datapath += ".dat"
    elif suffix != ".dat":
        raise ValueError("Molecular data file must either be a species name or species.dat")

    if not os.path.isdir(savedir):
        mkdir_p(savedir)

    if not os.path.isfile(datapath):
        data = lamda.query(species, return_datafile=True)
        with open(datapath,'w') as out:
            out.writelines(data)

    return os.path.split(datapath)
