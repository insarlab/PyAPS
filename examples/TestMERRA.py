import pyaps as pa
import numpy as np
import matplotlib.pyplot as plt
import sys

print '------------------------------------------------'
print 'You are using PyAPS from %s'%pa.__file__
print '------------------------------------------------'

print 'Testing Download Methods'
pa.ECMWFdload(['20040526','20030426'],'12','./ECMWF/')

print 'Testing MERRA Downloads'
pa.MERRAdload(['20040526','20030426'],'12','./MERRA/')

print 'Downloads OK'
print '------------------------------------------------'
print '------------------------------------------------'

print 'Testing ECMWF in Radar geometry, with a RMG dem'
aps1 = pa.PyAPS_rdr('ECMWF/ERA-Int_20030426_12.grb','dem_16rlks.hgt',grib='ECMWF',verb=True)
aps2 = pa.PyAPS_rdr('ECMWF/ERA-Int_20040526_12.grb','dem_16rlks.hgt',grib='ECMWF',verb=True)

print 'With Lat Lon files'
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getgeodelay(phs1,inc=23.0,wvl=0.056,lat='lat.flt',lon='lon.flt')
aps2.getgeodelay(phs2,inc=23.0,wvl=0.056,lat='lat.flt',lon='lon.flt')

LLphs = phs2-phs1
plt.subplot(121)
plt.imshow(LLphs[10:600,10:350])
plt.title('ECMWF')
plt.colorbar()

print 'Testing MERRA in Radar geometry, with a RMG dem'
aps1 = pa.PyAPS_rdr('MERRA/merra-20030426-12.hdf','dem_16rlks.hgt',grib='MERRA',verb=True)
aps2 = pa.PyAPS_rdr('MERRA/merra-20040526-12.hdf','dem_16rlks.hgt',grib='MERRA',verb=True)

print 'With Lat Lon files'
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getgeodelay(phs1,inc=23.0,wvl=0.056,lat='lat.flt',lon='lon.flt')
aps2.getgeodelay(phs2,inc=23.0,wvl=0.056,lat='lat.flt',lon='lon.flt')

phs = phs2-phs1
plt.subplot(122)
plt.imshow(phs[10:600,10:350])
plt.title('MERRA')
plt.colorbar()
plt.show()

print 'Radar Coordinates OK'
print '------------------------------------------------'
print '------------------------------------------------'

print 'Testing ECMWF in GEO geometry, with a FLT dem'
aps1 = pa.PyAPS_geo('ECMWF/ERA-Int_20030426_12.grb','Dem_crop.dem',grib='ECMWF',verb=True)
aps2 = pa.PyAPS_geo('ECMWF/ERA-Int_20040526_12.grb','Dem_crop.dem',grib='ECMWF',verb=True)

print 'With Lat Lon files'
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getdelay(phs1,inc=23.0,wvl=0.056)
aps2.getdelay(phs2,inc=23.0,wvl=0.056)

LLphs = phs2-phs1
plt.subplot(121)
plt.imshow(LLphs)
plt.title('ECMWF')
plt.colorbar()

print 'Testing MERRA in Geo geometry, with a FLT dem'
aps1 = pa.PyAPS_geo('MERRA/merra-20030426-12.hdf','Dem_crop.dem',grib='MERRA',verb=True)
aps2 = pa.PyAPS_geo('MERRA/merra-20040526-12.hdf','Dem_crop.dem',grib='MERRA',verb=True)

print 'With Lat Lon files'
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getdelay(phs1,inc=23.0,wvl=0.056)
aps2.getdelay(phs2,inc=23.0,wvl=0.056)

phs = phs2-phs1
plt.subplot(122)
plt.imshow(phs)
plt.title('MERRA')
plt.colorbar()
plt.show()




############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
