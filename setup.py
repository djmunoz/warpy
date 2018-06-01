from __future__ import print_function
from distutils.core import setup
import numpy as np
import os
import sys


setup(
    name='warpy',
    version="0.0.1",
    author='Diego J. Munoz',
    author_email = 'diego.munoz.anguita@gmail.com',
    url='https://github.com/',
    packages=['warpy'],
    description='Solver linear warp equations in the bending wave regime',
    install_requires = ['numpy','scipy'],
)
