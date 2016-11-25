#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

NAME = 'create_timelines'
VERSION = '0.3.0'
DESCRIPTION = 'Create timelines using TOAST'


# Utility function to read the README file.
def read(fname):
    import os
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='', top_path=None):
    config = Configuration(NAME, parent_package, top_path,
                           version=VERSION,
                           description=DESCRIPTION,
                           author='Maurizio Tomasi',
                           author_email='maurizio.tomasi@unimi.it',
                           license='MIT',
                           url='https://github.com/ziotom78/create_timelines',
                           long_description=read('README.md'),
                           requires=['numpy', 'healpy', 'astropy', 'toast', 'quaternionarray'])
    return config


if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
