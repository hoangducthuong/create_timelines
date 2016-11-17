#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

NAME = 'create_timelines'
VERSION = '0.1.0'
DESCRIPTION = 'Create timelines using TOAST'
FORTRAN2003_FLAG = '-std=f2003'


# Utility function to read the README file.
def read(fname):
    import os
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def configuration(parent_package='', top_path=None):
    extensions = [Extension('quat_to_pointings',
                            sources=['quat_to_pointings.f90'],
                            extra_f90_compile_args=[FORTRAN2003_FLAG])]

    config = Configuration(NAME, parent_package, top_path,
                           version=VERSION,
                           description=DESCRIPTION,
                           author='Maurizio Tomasi',
                           author_email='maurizio.tomasi@unimi.it',
                           license='MIT',
                           url='https://github.com/ziotom78/create_timelines',
                           long_description=read('README.md'),
                           requires=['numpy', 'healpy', 'astropy', 'toast'],
                           ext_modules=extensions)
    return config


if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
