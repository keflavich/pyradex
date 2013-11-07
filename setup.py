#!/usr/bin/env python
import sys

if 'develop' in sys.argv:
    # use setuptools for develop, but nothing else
    from setuptools import setup
else:
    from distutils.core import setup, Command

with open('README.rst') as file:
    long_description = file.read()

#with open('CHANGES') as file:
#    long_description += file.read()


version = "0.2"

import os
if not os.path.exists('pyradex/radex/radex.so'):
    import install_radex
    install_radex.install_radex()

import subprocess

class PyTest(Command):

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        errno = subprocess.call(['py.test'])
        raise SystemExit(errno)

setup(name='pyradex',
      version=version,
      description='Python-RADEX',
      long_description=long_description,
      author='Adam Ginsburg & Julia Kamenetzky',
      author_email='adam.g.ginsburg@gmail.com',
      url='http://github.com/keflavich/pyradex/',
      packages=['pyradex','pyradex.radex'],
      package_data={'pyradex.radex':['radex.so']},
      cmdclass={'test': PyTest},
      #include_package_data=True,
      )
