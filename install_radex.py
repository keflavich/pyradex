#raise "This was left in a terrible debug state may 25, 2021"
from __future__ import print_function
import tarfile
import sys
import re
import os
import shutil
import glob
from numpy import f2py
from subprocess import CompletedProcess
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
                   url='https://personal.sron.nl/~vdtak/radex/radex_public.tar.gz'):

    filename = 'radex_public.tar.gz'

    if os.path.isfile(filename) and not redownload:
        return filename

    print("Downloading RADEX")

    try:
        filename = download_file(url, cache=True)
        assert os.path.exists(filename)
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
    print("Extracting RADEX source from file {0}".format(filename))
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

    with open('Radex/src/readdata.f') as f:
        lines = f.readlines()

    with open('Radex/src/readdata.f','w') as f:
        # comment out the block dealing with ortho/para ratio: let python
        # handle that
        for ii,line in enumerate(lines):
            if ii <= 225 or ii > 235:
                f.write(line)
            else:
                f.write("c"+line[1:])
            #if 'density(3) = density(1)/(1.d0+1.d0/opr)' in line:
            #    f.write(line)
            #    f.write('c        For conservation of total density, set n(H2) = 0\n')
            #    f.write('         density(1) = 0.0\n')
            #else:
            #    f.write(line)



"""
Works for hpc:
    PATH=/Users/adam/repos/hpc/bin/:/usr/bin:~/virtual-python/bin/:/bin FFLAGS='-m64 -fPIC' CFLAGS='-fno-strict-aliasing -fno-automatic -fno-common -dynamic -m64 -g -O2' LDFLAGS='-m64 -undefined dynamic_lookup -bundle' python -c "import install_radex; install_radex.compile_radex(f77exec='/Users/adam/repos/hpc/bin/gfortran')"

Works for 4.2.3:
    FFLAGS='-arch i686 -arch x86_64 -fPIC' CFLAGS='-fno-strict-aliasing -fno-common -dynamic -arch i386 -arch x86_64 -g -O2' LDFLAGS='-arch i686 -arch x86_64 -undefined dynamic_lookup -bundle' python setup.py install_radex
"""

def compile_radex(fcompiler='gfortran',f77exec=None):
    #r1 = os.system('f2py -h pyradex/radex/radex.pyf Radex/src/*.f --overwrite-signature > radex_build.log')
    pwd = os.getcwd()
    os.chdir('Radex/src/')
    files = glob.glob('*.f')
    include_path = '--include-paths {0}'.format(os.getcwd())
    f2py.run_main(' -h radex.pyf --overwrite-signature'.split()+include_path.split()+files)

    if f77exec is None:
        f77exec=''
    else:
        f77exec = '--f77exec=%s' % f77exec
    #cmd = '-m radex -c %s --fcompiler=%s %s' % (" ".join(files), fcompiler, f77exec)
    #f2py.run_main(['-m','radex','-c','--fcompiler={0}'.format(fcompiler), f77exec,] + files)
    source_list = []
    for fn in files:
        with open(fn, 'r') as f:
            source_list.append(f.read())
    source = "\n".join(source_list)

    #with open("merged_source.f", 'wb') as fh:
    #    fh.write(source)

    include_path = '-I{0}'.format(os.getcwd())

    import platform
    mac_ver = platform.mac_ver()

    # not sure if this check is robust enough
    if int(mac_ver[0].split('.')[0]) >= 12:
        linker_path = '-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib'
    else:
        linker_path = ''

    extra_args = '--f77flags="-fno-automatic" --fcompiler={0} {1} {2} {3}'.format(fcompiler,
                                                                                  f77exec,
                                                                                  include_path,
                                                                                  linker_path)

    print(f"Running f2py with fcompiler={fcompiler}, f77exec={f77exec}, include_path={include_path}, linker_path={linker_path}")
    print(f"extra args = {extra_args}")
    print(f"Current directory = {os.getcwd()}")

    #f2py_path = os.path.join(sys.exec_prefix, 'bin', 'f2py')

    # Hack doesn't work.
    command = f'f2py -c -m radex {extra_args} *.f'
    # print(command)
    # import subprocess
    # r2 = subprocess.run(command.split(), stdout=subprocess.PIPE,
    #                     stderr=subprocess.PIPE, env=os.environ)

    r2 = f2py.compile(source="merged_source.f", modulename='radex',
                      verbose=True, full_output=False,
                      extra_args=extra_args,)
    print(f"Done running f2py in {os.getcwd()}.  r2={r2}")
    # 0 on success, or a subprocess.CompletedProcess if full_output=True
    if not isinstance(r2, CompletedProcess) and r2 != 0:
        r3 = f2py.compile(source="merged_source.f", modulename='radex',
                          verbose=True, full_output=True,
                          extra_args=extra_args,)
        raise SystemError(f"f2py failed with error {r2}\n"
                          "Try running the command: "
                          f"\n\ncd Radex/src/\n{command}\ncd -\nmv Radex/src/*so pyradex/radex/\n\n"
                          "See also Github issues 39 and 40")
    os.chdir(pwd)

    outfile = glob.glob("Radex/src/radex.*so")
    if len(outfile) != 1:
        print("outfile = {0}".format(outfile))
        raise OSError("Did not find the correct .so file(s)!  Compilation has failed.\n"
                      "Try running the command: "
                      f"\n\ncd Radex/src/\n{command}\ncd -\nmv Radex/src/*so pyradex/radex/\n\n"
                      "See also Github issues 39 and 40")
    sofile = outfile[0]
    r3 = shutil.move(sofile, 'pyradex/radex/radex.so')

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
