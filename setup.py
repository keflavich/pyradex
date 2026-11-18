#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import glob

if any((x in sys.argv for x in ('develop', 'bdist'))):
    # use setuptools for develop, but nothing else
    from setuptools import setup
    from distutils.core import Command
else:
    from distutils.core import setup, Command

with open('README.rst') as file:
    long_description = file.read()

#with open('CHANGES') as file:
#    long_description += file.read()

__version__ = ""
with open("pyradex/version.py") as f:
    exec(f.read())

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
        # Use pytest directly instead of deprecated astropy test runners
        errno = subprocess.call([sys.executable, '-m', 'pytest', 'pyradex/tests/'])
        raise SystemExit(errno)

radex_shared_object_files = [os.path.basename(fn) for fn in glob.glob("pyradex/radex/*.so")]
myradex_shared_object_files = [os.path.basename(fn) for fn in glob.glob("pyradex/fjdu/*.so")]

print(f"Found shared object files={radex_shared_object_files} for RADEX.  (if that is a blank, it means radex didn't install successfully)")
print(f"Found shared object files={myradex_shared_object_files} for RADEX.  (if that is a blank, it means fjdu's myradex didn't install successfully)")

setup(name='pyradex',
      version=__version__,
      description='Python-RADEX',
      long_description=long_description,
      author='Adam Ginsburg & Julia Kamenetzky',
      author_email='adam.g.ginsburg@gmail.com',
      url='http://github.com/keflavich/pyradex/',
      packages=['pyradex','pyradex.radex','pyradex.tests','pyradex.fjdu'],
      package_data={'pyradex.radex':['radex.so'] + radex_shared_object_files,
                    'pyradex.fjdu':['wrapper_my_radex.so'] + myradex_shared_object_files,
                    'pyradex.tests':['data/example.out']},
      requires=['requests', 'astroquery', ],
      install_requires=['astropy>=0.4.1', 'requests>=2.4.1',],
      cmdclass={'test': PyTest,
                'install_radex': InstallRadex,
                'build_radex_exe': BuildRadexExecutable,
                'install_myradex': InstallFjdu,
                'install_fjdu': InstallFjdu,
               },
      #include_package_data=True,
      zip_safe=False
      )

if os.getenv('RADEX_DATAPATH'):
    print("Installation has completed.  RADEX_DATAPATH={0}".format(os.getenv('RADEX_DATAPATH')))
else:
    print("Installation has completed.  Make sure to set your RADEX_DATAPATH variable!")
