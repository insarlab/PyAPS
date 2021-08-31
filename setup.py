############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################

# always prefer setuptools over distutils
from setuptools import setup, find_packages

setup(
    name='pyaps3',
    version='0.2.0',
    description="Python based Atmospheric Phase Screen Estimation",
    url="https://github.com/AngeliqueBenoit/pyaps3",
    author="Angelique Benoit, Romain Jolivet",
    author_email="insar@geologie.ens.fr",

    # package discovery
    packages=find_packages(),

    # dependencies
    install_requires=[
        'cdsapi',
        'ecCodes',
        'matplotlib',
        'netcdf4',
        'numpy',
        'pygrib',
        #'pyhdf',   #for MERRA
        'scipy',
    ],

    # data files
    include_package_data=True,
    package_data={'': ['*.cfg']},

) 
