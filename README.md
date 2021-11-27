[![Language](https://img.shields.io/badge/python-3.6%2B-blue.svg)](https://www.python.org/)
[![CircleCI](https://img.shields.io/circleci/build/github/insarlab/PyAPS.svg?logo=circleci&label=test)](https://circleci.com/gh/insarlab/PyAPS)
[![Version](https://img.shields.io/badge/version-v0.3.1-green.svg)](https://github.com/insarlab/PyAPS/releases)
[![License](https://img.shields.io/badge/license-GPLv3-yellow.svg)](https://github.com/insarlab/PyAPS/blob/main/LICENSE)
[![Citation](https://img.shields.io/badge/doi-10.1029%2F2011GL048757-blue)](https://doi.org/10.1029/2011GL048757)

## PyAPS - Python based Atmospheric Phase Screen estimation

This python 3 module estimates differential phase delay maps due to the stratified atmosphere for correcting radar interferograms. It is rewritten in Python 3 language from PYAPS source code and adapted for ECMWF's ERA-5 corrections. 

WARNING: The current version does not work with NARR and MERRA datasets. Contributions are welcomed.

This is research code provided to you "as is" with NO WARRANTIES OF CORRECTNESS. Use at your own risk.

### 1. Installation

`pyaps3` is available on the `conda-forge` channel and PyPI. The released version can be installed via `conda` as:

```bash
conda install -c conda-forge pyaps3
```

or via `pip` as:

```bash
pip install pyaps3
```

#### Build from source

The development version can be installed via `pip` as:

```bash
pip install git+https://github.com/insarlab/PyAPS.git
```

Or build from source manually as:

```bash
git clone https://github.com/insarlab/PyAPS.git
conda install -c conda-forge --file PyAPS/requirements.txt
python -m pip install -e PyAPS
```

Test the installation by running:

```bash
python PyAPS/tests/test_calc.py
```

### 2. Account setup for [ERA5](https://retostauffer.org/code/Download-ERA5/)

ERA5 data set is redistributed over the Copernicus Climate Data Store (CDS). Registration is required for the data access and downloading.

+ [Create a new account](https://cds.climate.copernicus.eu/user/register) on the CDS website if you don't own a user account yet. 
+ Create local key file. Create a file named `.cdsapirc` in your home directory and add the following two lines:

```shell
url: https://cds.climate.copernicus.eu/api/v2
key: 12345:abcdefghij-134-abcdefgadf-82391b9d3f
```

where 12345 is your personal user ID (UID), the part behind the colon is your personal API key. More details can be found [here](https://cds.climate.copernicus.eu/api-how-to).

+ **Make sure** that you accept the data license in the Terms of use on ECMWF website.

+ Test the account setup by running:

```bash
git clone https://github.com/insarlab/PyAPS.git --depth 1
python PyAPS/tests/test_dload.py
```

### 3. Citing this work

The methodology and validation can be found in:

+ Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311, doi:10.1029/2011GL048757.
