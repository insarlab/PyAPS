from pyhdf.HDF import *
from pyhdf.SD import *
from pyhdf.V import *
from pyhdf.VS import *
import numpy as np
import netCDF4


def get_merra(fname,minlat,maxlat,minlon,maxlon,cdic,verbose=False):
    '''Read data from MERRA hdf file. Note that the Lon values
       should be between [0-360]. Hdf file with weather model
       data can be downloaded from
       http://disc.sci.gsfc.nasa.gov/daac-bin/FTPSubset.pl

    Args:
        * fname   (str)       : Path to the grib file
        * minlat  (float)     : Minimum latitude
        * maxlat  (float)     : Maximum latitude
        * minlon  (float)     : Minimum longitude
        * maxlon  (float)     : Maximum longitude
        * cdic    (float)     : Dictionary of constants

    Returns:
        * lvls    (np.ndarray): Pressure levels
        * latlist (np.ndarray): Latitudes of the stations
        * lonlist (np.ndarray): Longitudes of the stations
        * gph     (np.ndarray): Geopotential height
        * tmp     (np.ndarray): Temperature
        * vpr     (np.ndarray): Vapor pressure

        '''

    if fname[-3:]=='hdf':
        print(minlat)
        print(maxlat)
        print(minlon)
        print(maxlon)

        # Read the hdf file
        file = SD(fname)
        if verbose:
            print('PROGRESS: READING HDF FILE')
        lvl = file.select('levels')
        # Pressure levels are from lowest to highest
        rlvls = lvl.get()
        lvls = []
        # Reverse the pressure levels to be consistent with other GAMss
        for i in range(len(rlvls)):
            index = len(rlvls) - i - 1
            lvls.append(rlvls[index])
        nlvls = len(lvls)
        lvls = np.array(lvls)

        alpha = cdic['Rv']/cdic['Rd']

        # Select latitutde and longitude
        lat = file.select('latitude')
        lon = file.select('longitude')
        lats = lat.get()
        lons = lon.get()
        mask1 = (lats > minlat) & (lats < maxlat)
        mask2 = (lons > minlon) & (lons < maxlon)
        [ii] = np.where(mask1 == True)
        [jj] = np.where(mask2 == True)
        del mask1
        del mask2
        iimemo = []
        for m in range(len(ii)):
            for i in range(len(jj)):
                iimemo.append(ii[m])
        jjmemo = []
        for i in range(len(ii)):
            jjmemo.append(jj)
        jjmemo = np.array(jjmemo)
        jjmemo = jjmemo.flatten()
        iimemo = np.array(iimemo)
        ii = iimemo
        jj = jjmemo
        latlist = lats[ii]
        lonlist = lons[jj]
        nstn = len(latlist)

        # Create arrays for 3D storage
        gph = np.zeros((nlvls, nstn))     #Potential height
        tmp = gph.copy()                  #Temperature
        vpr = gph.copy()                  #Vapor pressure
        if verbose:
            print('Number of stations:', nstn)

        # Get data from files
        h = file.select('h')
        height = h.get()[0]
        qv = file.select('qv')
        humidity = qv.get()[0]
        sp = file.select('ps')
        spressure = sp.get()[0]
        t = file.select('t')
        temp = t.get()[0]

        # Lvls are in hecto pascals, convert to pascals
        lvls = 100.0*lvls

        # Reverse altitude
        for i in range(nlvls):
            index = nlvls - i - 1
            gph[i,:] = height[index][ii,jj]

        idx = np.zeros(nstn)
        for i in range(nstn):
            for m in range(nlvls):
                if spressure[ii[i]][jj[i]] > lvls[m]:
                    idx[i] = m

        # extrapolation of temperature and humidity data at pressure levels under surface at each grid point
        tk = np.zeros(nstn)
        tb = np.zeros(nstn)
        for i in range(nstn):
            t = temp[:,ii[i],jj[i]]
            x = [lvls[idx[i]],lvls[idx[i] -1 ]]
            y = [t[nlvls - idx[i] - 1],t[nlvls - idx[i] ]]
            coef = np.polyfit(x,y,1)
            tk[i] = coef[0]
            tb[i] = coef[1]

        hk = np.zeros(nstn)
        hb = np.zeros(nstn)
        for i in range(nstn):
            hum = humidity[:,ii[i],jj[i]]
            x = [lvls[idx[i]],lvls[idx[i] -1 ]]
            y = [hum[nlvls - idx[i] - 1],hum[nlvls - idx[i]]]
            coef = np.polyfit(x,y,1)
            hk[i] = coef[0]
            hb[i] = coef[1]


        #fill out the tmp and vpr array
        for i in range(nstn):
            ind = int(idx[i]+1)
            for n in range(ind):
                tmp[n,i] = temp[nlvls - 1 - n,ii[i],jj[i]]
            exl = nlvls - ind
            for m in range(exl):
                tmp[ind+m,i] =  tk[i]*lvls[ind+m] + tb[i]

        for i in range(nstn):
            ind = int(idx[i]+1)
            for n in range(ind):
                vpr[n,i] = humidity[nlvls - 1- n,ii[i],jj[i]]
            exl = nlvls - ind
            for m in range(exl):
                vpr[ind+m,i] =  hk[i]*lvls[ind+m] + hb[i]

        memo = list(vpr)
        memo = np.array(memo)
        for i in range(nlvls):
            vpr[i,:] = memo[i,:]*lvls[i]*alpha/(1+(alpha - 1)*memo[i,:])

        # Close the hdf file
        file.end()

    if fname[-3:]=='nc4':
        # Read the netcdf file
        file = netCDF4.Dataset(fname)
        if verbose:
            print('PROGRESS: READING netcdf FILE')
        ncv = file.variables
        # Pressure levels are from lowest to highest
        rlvls = ncv['lev'][:]
        lvls = []
        # Reverse the pressure levels to be consistent with other GAMs
        for i in range(len(rlvls)):
            index = len(rlvls) - i - 1
            lvls.append(rlvls[index])
        nlvls = len(lvls)
        lvls = np.array(lvls)

        alpha = cdic['Rv']/cdic['Rd']

        # Select latitutde and longitude
        lats = ncv['lat'][:]
        lons = ncv['lon'][:]
        mask1 = (lats > minlat) & (lats < maxlat)
        mask2 = (lons > minlon) & (lons < maxlon)
        [ii] = np.where(mask1 == True)
        [jj] = np.where(mask2 == True)
        del mask1
        del mask2
        iimemo = []
        for m in range(len(ii)):
            for i in range(len(jj)):
                iimemo.append(ii[m])
        jjmemo = []
        for i in range(len(ii)):
            jjmemo.append(jj)
        jjmemo = np.array(jjmemo)
        jjmemo = jjmemo.flatten()
        iimemo = np.array(iimemo)
        ii = iimemo
        jj = jjmemo
        latlist = lats[ii]
        lonlist = lons[jj]
        nstn = len(latlist)

        # Create arrays for 3D storage
        gph = np.zeros((nlvls, nstn))     #Potential height
        tmp = gph.copy()                  #Temperature
        vpr = gph.copy()                  #Vapor pressure
        if verbose:
            print('Number of stations:', nstn)

        # Get data from files
        height = ncv['H'][:]
        height = height[0]
        humidity = ncv['QV'][:]
        humidity = humidity[0]
        spressure = ncv['PS'][:]
        spressure = spressure[0]
        temp = ncv['T'][:]
        temp = temp[0]

        # Lvls are in hecto pascals, convert to pascals
        lvls = 100.0*lvls

        # Reverse altitude
        for i in range(nlvls):
            index = nlvls - i - 1
            gph[i,:] = height[index][ii,jj]

        idx = np.zeros(nstn)
        for i in range(nstn):
            for m in range(nlvls):
                if spressure[ii[i]][jj[i]] > lvls[m]:
                    idx[i] = m

        # extrapolation of temperature and humidity data at pressure levels under surface at each grid point
        tk = np.zeros(nstn)
        tb = np.zeros(nstn)
        for i in range(nstn):
            t = temp[:,ii[i],jj[i]]
            x = [lvls[idx[i]],lvls[idx[i] -1 ]]
            y = [t[nlvls - idx[i] - 1],t[nlvls - idx[i] ]]
            coef = np.polyfit(x,y,1)
            tk[i] = coef[0]
            tb[i] = coef[1]

        hk = np.zeros(nstn)
        hb = np.zeros(nstn)
        for i in range(nstn):
            hum = humidity[:,ii[i],jj[i]]
            x = [lvls[idx[i]],lvls[idx[i] -1 ]]
            y = [hum[nlvls - idx[i] - 1],hum[nlvls - idx[i]]]
            coef = np.polyfit(x,y,1)
            hk[i] = coef[0]
            hb[i] = coef[1]


        #fill out the tmp and vpr array
        for i in range(nstn):
            ind = int(idx[i]+1)
            for n in range(ind):
                tmp[n,i] = temp[nlvls - 1 - n,ii[i],jj[i]]
            exl = nlvls - ind
            for m in range(exl):
                tmp[ind+m,i] =  tk[i]*lvls[ind+m] + tb[i]

        for i in range(nstn):
            ind = int(idx[i]+1)
            for n in range(ind):
                vpr[n,i] = humidity[nlvls - 1- n,ii[i],jj[i]]
            exl = nlvls - ind
            for m in range(exl):
                vpr[ind+m,i] =  hk[i]*lvls[ind+m] + hb[i]
        memo = list(vpr)
        memo = np.array(memo)
        for i in range(nlvls):
            vpr[i,:] = memo[i,:]*lvls[i]*alpha/(1+(alpha - 1)*memo[i,:])

        # Close the netcdf file
        file.close()

    # Send data
    return lvls,latlist,lonlist,gph,tmp,vpr


############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
