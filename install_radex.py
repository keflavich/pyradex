from __future__ import print_function
import tarfile
import re
import os
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

def download_radex(url='http://www.sron.rug.nl/~vdtak/radex/radex_public.tar.gz'):

    print("Downloading RADEX")

    try:
        filename = download_file(url, cache=True)
    except:
        import requests
        filename = 'radex_public.tar.gz'
        r = requests.get(url,
                         #data={'filename':filename},
                         stream=True,
                         verify=False)
        with open(filename,'wb') as f:
            for chunk in r.iter_content(1024):
                f.write(chunk)

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

def compile_radex():
    r1 = os.system('f2py -h pyradex/radex/radex.pyf Radex/src/*.f --overwrite-signature > radex_build.log')
    if r1 != 0:
        raise SystemError("f2py failed with error %i" % r1)
    r2 = os.system('f2py -m radex -c Radex/src/*.f --fcompiler=gfortran >> radex_build.log')
    if r2 != 0:
        raise SystemError("f2py failed with error %i" % r2)
    r3 = os.system('mv radex.so pyradex/radex/')
    if r3 != 0:
        raise SystemError("moving failed with error %i; radex.so was not created successfully" % r3)
