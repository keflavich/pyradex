# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

def get_package_data():
    paths_test = [os.path.join('data', '*.out')]
    return {'pyradex.tests': paths_test}
