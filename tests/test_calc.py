#!/usr/bin/env python
# Author: Zhang Yunjun, Nov 15, 2021


import os
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})

import pyaps3 as pa


print('------------------------------------------------')
print('import pyaps3 from {}'.format(pa.__file__))
print('------------------------------------------------')
print('test tropospheric delay calculation from ERA5.')
print('------------------------------------------------')
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# read geometry files
print('read ISCE geometry files: hgt/los/lat/lon.rdr')
dem = pa.utils.read_data(os.path.join(data_dir, 'hgt.rdr'))
inc = pa.utils.read_data(os.path.join(data_dir, 'los.rdr'), dname='inc')
lat = pa.utils.read_data(os.path.join(data_dir, 'lat.rdr'))
lon = pa.utils.read_data(os.path.join(data_dir, 'lon.rdr'))

# calculate
print('calculate tropospheric delay from GRB files...')
print('------------------------------------------------')
grb_file1 = os.path.join(data_dir, 'ERA5/ERA5_N30_N40_E120_E140_20101017_14.grb')
grb_file2 = os.path.join(data_dir, 'ERA5/ERA5_N30_N40_E120_E140_20110117_14.grb')
obj1 = pa.PyAPS(grb_file1, dem=dem, inc=inc, lat=lat, lon=lon, grib='ERA5', verb=True)
obj2 = pa.PyAPS(grb_file2, dem=dem, inc=inc, lat=lat, lon=lon, grib='ERA5', verb=True)
phs = obj2.getdelay() - obj1.getdelay()

# plot
date12 = '{}_{}'.format(grb_file1.split('_')[-2], grb_file2.split('_')[-2])
titles = [f'Tropospheric delay\n{date12}', 'Height']
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[10, 6])
for ax, data, title in zip(axs, [phs, dem], titles):
    im = ax.imshow(data, interpolation='nearest')
    # axis format
    ax.invert_yaxis()
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('meter')
fig.tight_layout()

print('------------------------------------------------')
print('Passed tropospheric delay calculation from ERA5.')
print('------------------------------------------------')
plt.show()
