'''
PyAPS module to compute InSAR phase delay maps from weather models.

Written by Romain Jolivet <rjolivet@gps.caltech.edu> and Piyush Agram <piyush@gps.caltech.edu>. The original Fortran package was written by Romain Jolivet and the Python version including support for different models was written by Piyush Agram.

.. note::
    Details of the python module can be obtained `here. <http://code.google.com/p/pyaps>`_
'''

__all__ = ['autoget','objects']

from pyaps3.objects import PyAPS
from pyaps3.autoget import ECMWFdload, MERRAdload, NARRdload

# get version info
from pyaps3.version import release_version as __version__

############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
