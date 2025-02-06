[![Language](https://img.shields.io/badge/python-3.8%2B-blue.svg?style=flat-square)](https://www.python.org/)
[![CircleCI](https://img.shields.io/circleci/build/github/insarlab/PyAPS.svg?logo=circleci&label=tests&style=flat-square)](https://circleci.com/gh/insarlab/PyAPS)
[![Conda Download](https://img.shields.io/conda/dn/conda-forge/pyaps3?color=green&label=conda%20downloads&style=flat-square)](https://anaconda.org/conda-forge/pyaps3)
[![Version](https://img.shields.io/github/v/release/insarlab/PyAPS?color=yellow&label=version&style=flat-square)](https://github.com/insarlab/PyAPS/releases)
[![License](https://img.shields.io/badge/license-GPLv3+-blue.svg?style=flat-square)](https://github.com/insarlab/PyAPS/blob/main/LICENSE)
[![Citation](https://img.shields.io/badge/doi-10.1029%2F2011GL048757-blue?style=flat-square)](https://doi.org/10.1029/2011GL048757)

## PyAPS - Python based Atmospheric Phase Screen estimation

This Python 3 module estimates differential phase delay maps due to the stratified atmosphere for correcting radar interferograms. It is rewritten in Python 3 language from PyAPS source code and adapted for ECMWF's ERA-5 corrections.

WARNING: The current version does not work with NARR and MERRA datasets. Contributions are welcomed.

This is research code provided to you "as is" with NO WARRANTIES OF CORRECTNESS. Use at your own risk.

### 1. Installation

#### a. Install the released version [recommended]

`pyaps3` is available on the [conda-forge](https://anaconda.org/conda-forge/pyaps3) channel, [PyPI](https://pypi.org/project/pyaps3/) and the main archive of the [Debian](https://tracker.debian.org/pkg/pyaps3) GNU/Linux OS. The released version can be installed via `conda` as:

```bash
conda install -c conda-forge pyaps3
```

or via `pip` as:

```bash
pip install pyaps3
```

or via `apt` (or other package managers) for [Debian-derivative OS](https://wiki.debian.org/Derivatives/Census) users, including [Ubuntu](https://ubuntu.com), as:

```bash
apt install python3-pyaps3
```

#### b. Install the development version

<p>
<details>
<p><summary>Click to expand for more details</summary></p>

The development version can be installed via `pip` as:

```bash
pip install git+https://github.com/insarlab/PyAPS.git
```

or build from source manually as:

```bash
git clone https://github.com/insarlab/PyAPS.git
conda install -c conda-forge --file PyAPS/requirements.txt --file PyAPS/tests/requirements.txt
python -m pip install -e PyAPS
```

Test the installation by running:

```bash
python PyAPS/tests/test_calc.py
```
</details>
</p>

### 2. Account setup for [ERA5](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5)

ERA5 data set is redistributed over the Copernicus Climate Data Store (CDS). Registration is required for the data access and downloading.

+ [Create a new account](https://cds.climate.copernicus.eu/) on the CDS website if you don't own a user account yet. Note: the old CDS account won't work ([Goodbye legacy CDS, Hellow New CDS!](https://forum.ecmwf.int/t/goodbye-legacy-climate-data-store-hello-new-climate-data-store-cds/6380)).
+ [CDS API setup](https://cds.climate.copernicus.eu/how-to-api): Create the local file `$HOME/.cdsapirc` (in your Unix/Linux environment) and add the following two lines:

```shell
url: https://cds.climate.copernicus.eu/api
key: your-personal-access-token
```

Your Personal Access Token can be found under [Your profile > Personal Access Token](https://cds.climate.copernicus.eu/profile) section or on the [setup guide](https://cds.climate.copernicus.eu/how-to-api) page. Alternatively, you could add the token to the `[CDS]` section in `model.cfg` file in the package directory, `site-packages/pyaps3` if installed via conda. Note: using your legacy CDS API key will lead to a 401 Client Error and Authentication failed.

+ **Make sure** that you accept the data license in the Terms of use on ECMWF website: Login, under [Datasets > ERA5 hourly data on pressure levels from 1940 to present > Download > Terms of use](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=download), click **Accept** to accespt the license to use Copernicus Products.

+ Test the account setup by running:

```bash
git clone https://github.com/insarlab/PyAPS.git --depth 1
python PyAPS/tests/test_dload.py
```

### 3. Citing this work

The methodology and validation can be found in:

+ Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric phase delay corrections from global meteorological reanalysis data, _Geophys. Res. Lett., 38,_ L17311, doi:[10.1029/2011GL048757](https://doi.org/10.1029/2011GL048757).
