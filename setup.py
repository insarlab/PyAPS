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
    version='0.2.1',
    description="Python based Atmospheric Phase Screen Estimation",
    url="https://github.com/insarlab/pyaps3",
    author="Romain Jolivet, Angelique Benoit",
    author_email="insar@geologie.ens.fr",

    # package discovery
    packages=find_packages("src"),  # include all packages under src
    package_dir={"": "src"},        # tell distutils packages are under src

    # dependencies
    install_requires=[
        'cdsapi',
        'ecCodes',
        'matplotlib',
        'netcdf4',
        'numpy',
        'pygrib',
        # 'pyhdf',   #for MERRA
        'scipy',
    ],

    # data files
    include_package_data=True,
    package_data={"pysolid": ["*.cfg"]},
) 
