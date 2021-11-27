# Author: Zhang Yunjun, Panji Brotoisworo, Jan 2021

# always prefer setuptools over distutils
from setuptools import setup, find_packages


# Grab from version.py file: version
with open("src/pyaps3/version.py", "r") as f:
    lines = f.readlines()
    line = [line for line in lines if line.strip().startswith("Tag(")][0].strip()
    version = line.replace("'",'"').split('"')[1]


# Grab from README file: long_description
with open("README.md", "r") as f:
    long_description = f.read()


setup(
    name='pyaps3',
    version=version,
    description="Python based Atmospheric Phase Screen Estimation",
    long_description=long_description,
    long_description_content_type="text/markdown",

    author="Romain Jolivet, Angelique Benoit",
    author_email="insar@geologie.ens.fr",

    license='GPL-3.0-or-later',
    llicense_files=('LICENSE',),

    url="https://github.com/insarlab/PyAPS",
    project_urls={
        "Bug Reports": "https://github.com/insarlab/PyAPS/issues",
        "Documentation": "https://github.com/insarlab/PyAPS/tree/main/docs",
        "Source": "https://github.com/insarlab/PyAPS",
    },

    # package discovery
    packages=find_packages("src"),  # include all packages under src
    package_dir={"": "src"},        # tell distutils packages are under src

    # dependencies
    install_requires=[
        'cdsapi',
        'matplotlib',
        'numpy',
        'pygrib',
        'scipy',
        # 'netcdf4', #for MERRA, which is currently not supported
        # 'pyhdf',   #for MERRA, which is currently not supported
    ],

    # data files
    include_package_data=True,
    package_data={"pyaps3": ["*.cfg"]},
)
