#!/usr/bin/env python
import sys

if 'develop' in sys.argv:
    # use setuptools for develop, but nothing else
    from setuptools import setup
    Command = object
else:
    from distutils.core import setup, Command

with open('README.rst') as file:
    long_description = file.read()

#with open('CHANGES') as file:
#    long_description += file.read()


execfile('pyradex/version.py')

#import os
#if not os.path.exists('pyradex/radex/radex.so'):
#    import install_radex
#    install_radex.install_radex()

class InstallRadex(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import install_radex
        install_radex.install_radex()

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

setup(name='pyradex',
      version=version,
      description='Python-RADEX',
      long_description=long_description,
      author='Adam Ginsburg & Julia Kamenetzky',
      author_email='adam.g.ginsburg@gmail.com',
      url='http://github.com/keflavich/pyradex/',
      packages=['pyradex','pyradex.radex','pyradex.tests'],
      package_data={'pyradex.radex':['radex.so'],
                    'pyradex.tests':['data/example.out']},
      cmdclass={'test': PyTest, 'install_radex': InstallRadex},
      #include_package_data=True,
      )
