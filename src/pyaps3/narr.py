import numpy as np
import pygrib
import scipy.interpolate as si
import sys

#############Clausius-Clapeyron for NARR ###########
def cc_narr(tmp,cdic):
    '''Clausius Clayperon law used in the NARR model.

    Args:
        * tmp  (np.ndarray): Temperature.
        * cdic (dict)      : Dictionnary of constants.

    Returns:
        * esat (np.ndarray): Saturation water vapor partial pressure.'''

    a1w=cdic['a1w']
    a3w=cdic['a3w']
    a4w=cdic['a4w']
    T3=cdic['T3']
    Rv=cdic['Rv']
    esat = a1w*np.exp(a3w*(tmp-T3)/(tmp-a4w))

    return esat
###############Completed CC_NARR#####################################


########Read in ERA data from a given ERA Interim file##################
def get_narr(fname,minlat,maxlat,minlon,maxlon,cdic,verbose=False):
    '''Read data from NARR grib file. Note that the Lon values should be between [0-360].
    GRB file with weather model data can be downloaded from http://nomads.ncdc.noaa.gov/data/narr .

    Args:
        * fname  (str)  :  Path to the grib file
        * minlat (float):  Minimum latitude
        * maxlat (float):  Maximum latitude
        * minlon (float):  Minimum longitude
        * maxlon (float):  Maximum longitude
        * cdic   (float):  Dictionary of constants
    
    Kwargs:
        * humidity (str): Specific ('Q') or relative humidity ('R').

    Returns:
        * lvls    (np.ndarray): Pressure levels
        * latlist (np.ndarray): Latitudes of the stations
        * lonlist (np.ndarray): Longitudes of the stations
        * gph     (np.ndarray): Geopotential height
        * tmp     (np.ndarray): Temperature
        * vpr     (np.ndarray): Vapor pressure

    .. note::
        Uses cc_narr by default.
        '''

    if verbose:
        print('PROGRESS: READING GRIB FILE')
    lvls = np.array([100, 125, 150, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550,
                     600, 650, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000])
    nlvls = len(lvls)

    alpha = cdic['Rv']/cdic['Rd']
    gphind = np.array([16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,137,
                       146,155,164,173,182,191,200,210,219,228,237,246,255])
    
    grbs = pygrib.open(fname)
    grbs.seek(gphind[0])
    grb=grbs.read(1)[0]
    lats,lons = grb.latlons()
    lons[lons<0] += 360.
    g = cdic['g']    
    mask = (lats > minlat) & (lats < maxlat) & (lons > minlon) & (lons < maxlon)
    [ii,jj] = np.where(mask == True)
    del mask
    latlist = lats[ii,jj]
    lonlist = lons[ii,jj]
    nstn = len(ii)
    
    ####Create arrays for 3D storage
    gph = np.zeros((nlvls, nstn))     #Potential height
    tmp = gph.copy()                  #Temperature
    vpr = gph.copy()                  #Vapor pressure
    if verbose:
        print('Number of stations:', nstn)

    lvls = 100.0*lvls                 #Conversion to absolute pressure
    for i in range(nlvls):
        grbs.seek(gphind[i])          #Reading potential height.
        grb = grbs.read(3)
        val = grb[0].values
        gph[i,:] = val[ii,jj]

        val = grb[1].values           #Reading temperature
        temp = val[ii,jj]
        tmp[i,:] = temp

                
        val = grb[2].values  #Specific humidity
        temp = val[ii,jj]
        vpr[i,:] = temp*lvls[i]*alpha/(1+(alpha - 1)*temp)
        
    return lvls,latlist,lonlist,gph,tmp,vpr
###############Completed GET_ERA########################################


########Interpolates the NARR delay to a regular grid####################
def intdel(hgt,latlin,lonlin,delcin,spacing=0.3):
    '''Interpolates the NARR data to a regular grid with a grid spacing being the average of previous grid spaces

    Args:
            * hgt     (np.ndarray): Altitude levels
            * latlin  (np.ndarray): Latitudes of the stations
            * lonlin  (np.ndarray): Longitudes of the stations
            * delcin  (np.ndarray): Delay cube len(latlist)xlen(hgt)
    
    Returns:
            * latlist (np.ndarray): Latitudes of the stations
            * lonlist (np.ndarray): Longitudes of the stations
            * delc    (np.ndarray): Delay cube len(latlist)xlen(hgt)

    .. note::
    '''

    # Points array
    nstn = len(latlin)
    Points = np.zeros((nstn,2))
    Points[:,0] = lonlin
    Points[:,1] = latlin

    # Output arrays
    minlat = latlin.min()
    maxlat = latlin.max()
    minlon = lonlin.min()
    maxlon = lonlin.max()
    latlist = np.arange(minlat,maxlat,spacing)
    lonlist = np.arange(minlon,maxlon,spacing)
    [lonlist,latlist] = np.meshgrid(lonlist,latlist)
    latlist = latlist.flatten()
    lonlist = lonlist.flatten()

    # Xi arrays
    nstno = len(latlist)
    Xi = np.zeros((nstno,2))
    Xi[:,0] = lonlist
    Xi[:,1] = latlist

    # Create the arrays
    nlvls = len(hgt)
    delc = np.zeros((nstno,nlvls))

    # Loop on the pressure levels
    for n in range(nlvls):
        delc[:,n] = si.griddata(
            Points,
            delcin[:,n],
            Xi,
            method='cubic',
            #fill_value=10*float(n),
        )

    return delc,latlist,lonlist


############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
############################################################
