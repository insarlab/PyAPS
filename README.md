## PyAPS - Python based Atmospheric Phase Screen estimation

This python 3 module estimates differential phase delay maps due to the stratified atmosphere for correcting radar interferograms. It is rewritten in Python 3 language from PYAPS source code and adapted for ECMWF's ERA-5 corrections. 

WARNING: PyAPS3 could not work with NARR and MERRA weather models. Feel free to modify these parts of the code.

### 1. Installation

#### via pip

```bash
pip install git+https://github.com/yunjunz/PyAPS.git
```

#### via conda

```bash
# download source code
cd ~/tools
git clone https://github.com/yunjunz/PyAPS.git

# install dependencies
conda install --file PyAPS/requirements.txt
```

Add the following variables in your source file (**_~/.bash_profile_** for _bash_ users or **_~/.cshrc_** for _csh/tcsh_ users).

```bash
export PYTHONPATH=${PYTHONPATH}:~/tools/PyAPS
```

### 2. Account setup for global atmospheric models

#### [ERA-5](https://retostauffer.org/code/Download-ERA5/)
ERA-5 data set is redistributed over the Copernicus Climate Data Store (CDS), [create a new account](https://cds.climate.copernicus.eu/user/register) on the CDS website if you don't own a user account yet. On the profile, you will find your user id (**UID**) and your personal **API Key**. Add the key to the `model.cfg` file as below:

```
#####key for ECMWF in Climate Data Store Application Program Interface
#Get it from https://retostauffer.org/code/Download-ERA5/
[CDS]
key = 12345:abcdefghij-134-abcdefgadf-82391b9d3f
```

where 12345 is your personal user ID (UID), the part behind the colon your personal API key. More details on CDSAPI can be found [here](https://cds.climate.copernicus.eu/api-how-to).

Run `examples/TestECMWF.ipynb` in Jupyter Notebook in your local machine to check if everything works.

### 3. Citing this work
The metholody and validation of DelayPackage can be found in:
+ Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311, doi:10.1029/2011GL048757.

Examples in the example directory.
