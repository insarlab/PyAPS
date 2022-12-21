'''
PyAPS module to compute InSAR phase delay maps from weather models.

Written by Romain Jolivet <rjolivet@gps.caltech.edu> and Piyush Agram <piyush@gps.caltech.edu>.
The original Fortran package was written by Romain Jolivet and the Python version including
support for different models was written by Piyush Agram.

.. note::
    Details of the python module can be obtained `here. <http://code.google.com/p/pyaps>`
'''

# get version info [requires python >=3.8]
from importlib.metadata import PackageNotFoundError, version
try:
    __version__ = version(__name__)
except PackageNotFoundError:
    print('package is not installed!\n'
          'Please follow the installation instructions in the README.md.\n'
          'Or, to just get the version number, use:\n'
          '   python -m setuptools_scm')

# top-level functions
from pyaps3.objects import PyAPS
from pyaps3.autoget import ECMWFdload, MERRAdload, NARRdload

__all__ = ['__version__', 'PyAPS', 'ECMWFdload', 'MERRAdload', 'NARRdload']

############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
