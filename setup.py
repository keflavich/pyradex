#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import builtins

# Ensure that astropy-helpers is available
import ah_bootstrap  # noqa

from setuptools import setup
from setuptools.config import read_configuration

from astropy_helpers.setup_helpers import register_commands, get_package_info
from astropy_helpers.version_helpers import generate_version_py

# Store the package name in a built-in variable so it's easy
# to get from other parts of the setup infrastructure
builtins._ASTROPY_PACKAGE_NAME_ = read_configuration('setup.cfg')['metadata']['name']

# Create a dictionary with setup command overrides. Note that this gets
# information about the package (name and version) from the setup.cfg file.
cmdclass = register_commands()

# Freeze build information in version.py. Note that this gets information
# about the package (name and version) from the setup.cfg file.
version = generate_version_py()

# Get configuration information from all of the various subpackages.
# See the docstring for setup_helpers.update_package_files for more
# details.
package_info = get_package_info()


import sys
import glob

if any((x in sys.argv for x in ('develop','bdist'))):
    # use setuptools for develop, but nothing else
    from setuptools import setup
    from distutils.core import Command
else:
    from distutils.core import setup, Command

with open('README.rst') as file:
    long_description = file.read()

#with open('CHANGES') as file:
#    long_description += file.read()

#import os
#if not os.path.exists('pyradex/radex/radex.so'):
#    import install_radex
#    install_radex.install_radex()

class InstallRadex(Command):
    """
    Compile the f2py radex modules needed for the python wrapper
    """

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import install_radex
        install_radex.install_radex()

class BuildRadexExecutable(Command):
    """
    Create the files radex_sphere, radex_lvg, radex_slab in the Radex/bin/
    directory
    """

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import install_radex
        install_radex.build_radex_executable()


class InstallFjdu(Command):
    """
    Compile Fujun Du's "myradex"
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        cwd = os.getcwd()
        os.chdir('myRadex')
        if os.path.exists('wrapper_my_radex.so'):
            os.remove('wrapper_my_radex.so')
        os.system('make clean')
        os.system('make wrapper')
        result = os.system('make sub_trivials.o')
        if result != 0:
            raise ValueError("Compilation has failed.  Check gfortran version?")
        os.chdir(cwd)
        for fn in glob.glob('myRadex/wrapper_my_radex*.so'):
            outpath = 'pyradex/fjdu/{0}'.format(os.path.basename(fn))
            if os.path.exists(outpath):
                os.remove(outpath)
            os.link(fn, outpath)
        os.chdir('myRadex')
        os.system('make clean')
        os.chdir(cwd)


import subprocess
import shutil
import os

class PyTest(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):

        if os.path.exists('build'):
            shutil.rmtree('build')
        #errno1 = subprocess.call(['py.test','--genscript=runtests.py'])
        errno2 = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno2)

#setup(version=version, cmdclass=cmdclass, **package_info)
cmdclass.update({
                'install_radex': InstallRadex,
                'build_radex_exe': BuildRadexExecutable,
                'install_myradex': InstallFjdu,
                'install_fjdu': InstallFjdu,
               })
package_info['package_data'].update(
    {'pyradex.radex':['radex.so'],
     'pyradex.fjdu':['wrapper_my_radex.so'],
     'pyradex.tests':['data/example.out']},)
package_info['packages'] += ['pyradex','pyradex.radex','pyradex.tests','pyradex.fjdu']

setup(name='pyradex',
      version=version,
      #description='Python-RADEX',
      #long_description=long_description,
      #author='Adam Ginsburg & Julia Kamenetzky',
      #author_email='adam.g.ginsburg@gmail.com',
      #url='http://github.com/keflavich/pyradex/',
      requires=['requests', 'astroquery', ],
      install_requires=['astropy>=0.4.1', 'requests>=2.4.1',],
      cmdclass=cmdclass,
      **package_info
      )

if os.getenv('RADEX_DATAPATH'):
    print("Installation has completed.  RADEX_DATAPATH={0}".format(os.getenv('RADEX_DATAPATH')))
else:
    print("Installation has completed.  Make sure to set your RADEX_DATAPATH variable!")
