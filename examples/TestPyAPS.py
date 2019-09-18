#!/usr/bin/env python
import pyaps3 as pa
import numpy as np
import matplotlib.pyplot as plt
import sys

print('------------------------------------------------')
print('You are using PyAPS from %s'%pa.__file__)
print('------------------------------------------------')

print('Testing Download Methods')
print('Testing ECMWF Downloads')
pa.ECMWFdload(['20150503'],'00','./ECMWF','interim','fc')

print 'Testing NARR Downloads'
pa.NARRdload(['20040526','20030426'],'12','./NARR/')

print 'Testing ERA@UCAR Downloads'
pa.ERAdload(['20050708','20041201'],'12','./ERA/')

print('Downloads OK')
print('------------------------------------------------')
print('------------------------------------------------')

print('Testing ECMWF in Radar geometry, with a RMG dem')
aps1 = pa.PyAPS_rdr('ECMWF/ERA-Int_20150503_12.grb','/data/angel/TestISCE/stack/INTERFEROGRAMS/hgt.9alks_27rlks.full.rdr',grib='ECMWF',verb=True)
aps2 = pa.PyAPS_rdr('ECMWF/ERA-Int_20040526_12.grb','dem_16rlks.hgt',grib='ECMWF',verb=True)

print('With Lat Lon files')
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getgeodelay(phs1,inc='/data/angel/TestISCE/stack/INTERFEROGRAMS/losFull.9alks_27rlks.rdr',wvl=0.05546576,lat='/data/angel/TestISCE/stack/INTERFEROGRAMS/latFull.9alks_27rlks.rdr',lon='/data/angel/TestISCE/stack/INTERFEROGRAMS/lonFull.9alks_27rlks.rdr')
aps2.getgeodelay(phs2,inc=23.0,wvl=0.056,lat='lat.flt',lon='lon.flt')

LLphs = phs2-phs1
plt.subplot(121)
plt.imshow(LLphs[10:600,10:350])
plt.title('With latlon file')
plt.colorbar()
plt.show()

print('Testing ECMWF in Radar geometry, with a FLT dem')
aps1 = pa.PyAPS_rdr('ECMWF/ERA-Int_20030426_12.grb','DEM.flt',grib='ECMWF',demfmt='HGT')
aps2 = pa.PyAPS_rdr('ECMWF/ERA-Int_20040526_12.grb','DEM.flt',grib='ECMWF',demfmt='HGT')

print('With Lat Lon files')
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getgeodelay(phs1,inc=23.0,wvl=0.056,lat='lat.flt',lon='lon.flt')
aps2.getgeodelay(phs2,inc=23.0,wvl=0.056,lat='lat.flt',lon='lon.flt')

LLphs = phs2-phs1
plt.subplot(121)
plt.imshow(LLphs[10:600,10:350])
plt.title('With latlon file')
plt.colorbar()
plt.show()

print('Testing ECMWF in Geographic coordinates, with a FLT dem')
aps1 = pa.PyAPS_geo('ECMWF/ERA-Int_20030426_12.grb','Dem_crop.dem',grib='ECMWF')
aps2 = pa.PyAPS_geo('ECMWF/ERA-Int_20040526_12.grb','Dem_crop.dem',grib='ECMWF')

print('With Lat Lon files')
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getdelay(phs1,inc=23.0,wvl=0.056)
aps2.getdelay(phs2,inc=23.0,wvl=0.056)

del phs1
del phs2

print('ECMWF OK')
print('------------------------------------------------')
print('------------------------------------------------')

sys.exit(1)

print('Testing NARR in Radar geometry, with a RMG dem')
aps1 = pa.PyAPS_rdr('NARR/narr-a_221_20030426_1200_000.grb','hgt.flt',grib='NARR',demfmt='HGT')
aps2 = pa.PyAPS_rdr('NARR/narr-a_221_20040526_1200_000.grb','hgt.flt',grib='NARR',demfmt='HGT')
print('With Lat Lon files')
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getgeodelay(phs1,inc=23.0,wvl=0.056,lat='lat2.flt',lon='lon2.flt')
aps2.getgeodelay(phs2,inc=23.0,wvl=0.056,lat='lat2.flt',lon='lon2.flt')

phs = phs2-phs1
plt.imshow(phs)
plt.colorbar()
plt.show()

print('------------------------------------------------')
print('------------------------------------------------')

print('All Test OK')



############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
