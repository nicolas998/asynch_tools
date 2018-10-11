#!/usr/bin/env python
import os
from numpy.distutils.core import setup, Extension

setup(
    name='asynch_tools',
    version='0.0.1',
    author='Nicolas Velasquez G',
    author_email='nicolas.velasquezgiron@gmail.com',    
    packages=['asynch_tools'],
    package_data={'asynch_tools':['inputs.py','outputs.py']},
    url='https://github.com/nicolas998/asynch_tools.git',
    license='LICENSE.txt',
    description='Tools to work with inputs and outputs of asynch in python',
    long_description=open('README.txt').read(),
    install_requires=[ ],
    #ext_modules=[ext1, ext2],
	)
