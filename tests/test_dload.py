#!/usr/bin/env python
# Author: Zhang Yunjun, Nov 15, 2021

import os
import pyaps3 as pa

print('------------------------------------------------')
print('import pyaps3 from {}'.format(pa.__file__))
print('------------------------------------------------')
print('test ERA5 data download')
print('NOTE: Account setup is required on the Copernicus Climate Data Store (CDS).')
print('      More detailed info can be found on: https://cds.climate.copernicus.eu/how-to-api')
print('      Add your account info to ~/.cdsapirc file.')
filedir = os.path.join(os.path.dirname(__file__), 'data', 'ERA5')
pa.ECMWFdload(['20200601','20200901'], hr='14', filedir=filedir, model='ERA5', snwe=(30,40,120,140))

print('------------------------------------------------')
print('Downloads OK')
print('------------------------------------------------')
