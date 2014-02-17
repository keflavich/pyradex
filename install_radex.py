from __future__ import print_function
import tarfile
import re
import os
import shutil
try:
    from astropy.utils.data import download_file
except ImportError:
    download_file = None


def install_radex(download=True, extract=True, patch=True, compile=True):
    if download:
        filename = download_radex()
    if extract:
        extract_radex(filename)
    if patch:
        patch_radex()
    if compile:
        compile_radex()

def download_radex(redownload=True,
                   url='http://www.sron.rug.nl/~vdtak/radex/radex_public.tar.gz'):

    filename = 'radex_public.tar.gz'

    if os.path.isfile(filename) and not redownload:
        return filename

    print("Downloading RADEX")

    try:
        filename = download_file(url, cache=True)
    except:
        import requests
        r = requests.get(url,
                         #data={'filename':filename},
                         stream=True,
                         verify=False)
        with open(filename,'wb') as f:
            for chunk in r.iter_content(1024):
                f.write(chunk)

    print("Download succeeded, or at least didn't obviously fail.")

    return filename

def extract_radex(filename='radex_public.tar.gz'):
    with tarfile.open(filename,mode='r:gz') as tar:
        # this will fail if readonly files are overwritten: tf.extractall()
        for file_ in tar:
            try:
                tar.extract(file_)
            except IOError as e:
                os.remove(file_.name)
                tar.extract(file_)
            finally:
                os.chmod(file_.name, file_.mode)

def patch_radex():
    # PATCHES
    radlines = []
    vers = re.compile("parameter\(version = '([a-z0-9]*)'\)")
    with open('Radex/src/radex.inc') as f:
        for line in f.readlines():
            if ('parameter(radat' in line or
                'parameter(version' in line or
                'parameter(logfile' in line or
                'parameter (method' in line):
                line = 'c'+line
            if vers.search(line):
                radlines.append("c      version = '%s'\n" % vers.search(line).groups()[0])
            if 'parameter(debug' in line:
                line = 'c'+line
                radlines.append('      common/dbg/debug\n')
            radlines.append(line)

            if 'parameter (method = 3)' in line:
                radlines.append('      common/setup/radat,method,version,logfile\n')

    with open('Radex/src/radex.inc','w') as f:
        for line in radlines:
            f.write(line)

    with open('Radex/src/background.f') as f:
        lines = f.readlines()

    with open('Radex/src/background.f','w') as f:
        for line in lines:
            if 'parameter(huge' in line:
                f.write('       parameter(huge=1.0e38)\n')
                f.write('c      ! highest allowed by f90 (fdvt 28apr06)\n')
            else:
                f.write(line)

"""
Works for hpc:
    PATH=/Users/adam/repos/hpc/bin/:/usr/bin:~/virtual-python/bin/:/bin FFLAGS='-m64 -fPIC' CFLAGS='-fno-strict-aliasing -fno-common -dynamic -m64 -g -O2' LDFLAGS='-m64 -undefined dynamic_lookup -bundle' python -c "import install_radex; install_radex.compile_radex(f77exec='/Users/adam/repos/hpc/bin/gfortran')"

Works for 4.2.3:
    FFLAGS='-arch i686 -arch x86_64 -fPIC' CFLAGS='-fno-strict-aliasing -fno-common -dynamic -arch i386 -arch x86_64 -g -O2' LDFLAGS='-arch i686 -arch x86_64 -undefined dynamic_lookup -bundle' python setup.py install_radex
"""

def compile_radex(fcompiler='gfortran',f77exec=None):
    r1 = os.system('f2py -h pyradex/radex/radex.pyf Radex/src/*.f --overwrite-signature > radex_build.log')
    if r1 != 0:
        raise SystemError("f2py failed with error %i" % r1)

    if f77exec is None:
        f77exec=''
    else:
        f77exec = '--f77exec=%s' % f77exec
    cmd = 'f2py -m radex -c Radex/src/*.f --fcompiler=%s %s >> radex_build.log' % (fcompiler,f77exec)
    r2 = os.system(cmd)
    if r2 != 0:
        raise SystemError("f2py failed with error %i" % r2)

    r3 = os.system('mv radex.so pyradex/radex/')
    if r3 != 0:
        raise SystemError("moving failed with error %i; radex.so was not created successfully" % r3)

def build_radex_executable(datapath='./'):
    filename = download_radex(redownload=False)
    # need to re-extract the RADEX source to get an un-patched version
    extract_radex(filename)
    compile_radex_source(datapath=datapath)


def compile_radex_source(datapath='./'):
    """
    Compile the source file in the Radex/ directory

    May be good to download & untar first
    """

    cwd = os.getcwd()

    # Set datapath to the 'examples' directory if not specified (since there is
    # at least co.dat there)
    if datapath is None:
        datapath = os.path.join(os.getcwd(), 'examples')
        print("Datapath was not defined.  Set to ",datapath)

    # make sure ~ gets expanded
    datapath = os.path.expanduser(datapath)

    try:
        os.link('Radex/data/hco+.dat',os.path.join(datapath,'hco+.dat'))
    except OSError:
        pass

    # I think fortran requires a trailing slash... can't hurt
    if datapath[-1] != '/':
        datapath = datapath+'/'

    os.chdir('Radex/src/')
    with open('radex.inc','r') as f:
        lines = [L.replace('/Users/floris/Radex/moldat/',datapath)
                 if 'radat' in L else L
                 for L in f.readlines()]
    with open('radex.inc','w') as of:
        of.writelines(lines)

    method_types = {1:'sphere',2:'lvg',3:'slab'}
    for method in method_types:
        radex_inc_method('./',method=method)
        r1 = os.system('make')
        if r1 != 0:
            raise SystemError("radex make failed with error %i" % r1)
        shutil.move('../bin/radex','../bin/radex_%s' % method_types[method])

    os.chdir(cwd)

def radex_inc_method(datapath, method=1):
    """
    Convert the radex.inc file to a method

    Parameters
    ----------
    datapath: path
        a directory path containing the target radex.inc
    method: 1,2,3
        1: sphere
        2: lvg
        3: slab
    """
    fn = os.path.join(datapath, 'radex.inc')
    with open(fn,'r') as f:
        lines = [L.replace('/Users/floris/Radex/moldat/',datapath)
                 if 'radat' in L else L
                 for L in f.readlines()]
    
    radlines = []
    for line in lines:
        if ('parameter (method' in line):
            line = 'c'+line
        if 'method = 3' in line:
            radlines.append('      parameter (method = %i)\n' % method)
        radlines.append(line)

    with open(fn,'w') as f:
        f.writelines(radlines)

