# pyaps3
Python 3 Atmospheric Phase Screen

This python 3 module estimates differential phase delay maps due to the stratified atmosphere for correcting radar interferograms. It is rewritten in Python 3 language from PYAPS source code and adapted for ECMWF's ERA-5 corrections. 

Warning: PyAPS3 could not work with NARR and MERRA weather models. Feel free to modify these parts of the code.

## Installation

git clone 

wget 
https://pypi.org/project/cdsapi/


## Account setup for global atmospheric models
### [ERA-5](https://retostauffer.org/code/Download-ERA5/)
ERA-5 data set is redistributed over the Copernicus Climate Data Store (CDS), [create a new account](https://cds.climate.copernicus.eu/user/register) on the CDS website if you don't own a user account yet. On the profile, you will find your user id (**UID**) and your personal **API Key**.

For batch script data downloads you’ll have to create a local ASCII file with your user information (UID, API key) which is used by the python package (cdsapi). To do so (linux) simply create a file called .cdsapirc in your home directory and add the following two lines:

```
url: https://cds.climate.copernicus.eu/api/v2
key: 1234:abcdefghij-134-abcdefgadf-82391b9d3f
```

where 1234 is your personal user ID (UID), the part behind the colon your personal API key. Line one simply contains the URL to the web API. More details can be found [here](https://cds.climate.copernicus.eu/api-how-to).

## Citing this work
The metholody and validation of DelayPackage can be found in:
+ Jolivet, R., R. Grandin, C. Lasserre, M.-P. Doin and G. Peltzer (2011), Systematic InSAR tropospheric phase delay corrections from global meteorological reanalysis data, Geophys. Res. Lett., 38, L17311, doi:10.1029/2011GL048757.

Examples in the example directory.
