## PyAPS - Python based Atmospheric Phase Screen estimation

This python 3 module estimates differential phase delay maps due to the stratified atmosphere for correcting radar interferograms. It is rewritten in Python 3 language from PYAPS source code and adapted for ECMWF's ERA-5 corrections. 

WARNING: PyAPS3 could not work with NARR and MERRA weather models. Contributions are welcomed.

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

One could also installed the development version via `pip` as:

```bash
pip install git+https://github.com/insarlab/PyAPS.git
```

Run the following to 

### 2. Account setup for [ERA5](https://retostauffer.org/code/Download-ERA5/)

ERA5 data set is redistributed over the Copernicus Climate Data Store (CDS). Registration is required for the data downloading.

+ [Create a new account](https://cds.climate.copernicus.eu/user/register) on the CDS website if you don't own a user account yet. 
+ Create local key file. Create a file named `.cdsapirc` in your home directory and add the following two lines:

```shell
url: https://cds.climate.copernicus.eu/api/v2
key: 12345:abcdefghij-134-abcdefgadf-82391b9d3f
```

where 12345 is your personal user ID (UID), the part behind the colon is your personal API key. More details can be found [here](https://cds.climate.copernicus.eu/api-how-to).

+ **Make sure** that you accept the data license in the Terms of use on ECMWF website.

### 3. Testing

Run the following to test the installation and account setup:

```bash
git clone https://github.com/insarlab/PyAPS.git --depth 1
python PyAPS/tests/test_dload.py
python PyAPS/tests/test_calc.py
```

### 4. Citing this work

The metholody and validation of DelayPackage can be found in:

+ Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311, doi:10.1029/2011GL048757.
