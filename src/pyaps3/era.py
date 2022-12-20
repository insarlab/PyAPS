############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
import numpy as np
import pygrib


#############Clausis-Clayperon for ECMWF###########################
def cc_era(tmp,cdic):
    '''Clausius Clayperon law used by ERA Interim.

    Args:
        * tmp  (np.ndarray) : Temperature.
        * cdic (dict)       : Dictionary of constants

    Returns:
        * esat (np.ndarray) : Water vapor saturation partial pressure.'''


    a1w = cdic['a1w']
    a3w = cdic['a3w']
    a4w = cdic['a4w']
    a1i = cdic['a1i']
    a3i = cdic['a3i']
    a4i = cdic['a4i']
    T3  = cdic['T3']
    Ti  = cdic['Ti'] 

    esatw = a1w*np.exp(a3w*(tmp-T3)/(tmp-a4w))
    esati = a1i*np.exp(a3i*(tmp-T3)/(tmp-a4i))
    esat = esati.copy()
    for k in range(len(tmp)):
        if (tmp[k] >= T3):
            esat[k] = esatw[k]
        elif (tmp[k] <= Ti):
            esat[k] = esati[k]
        else:
            wgt = (tmp[k]-Ti)/(T3-Ti)
            esat[k] = esati[k] + (esatw[k]-esati[k])*wgt*wgt

    return esat


########Read in ERA data from a given ERA Interim file##################
def get_era(fname,minlat,maxlat,minlon,maxlon,cdic, humidity='Q',verbose=False):
    '''Read data from ERA interim grib file. 
    Note that the Lon values should be between [0-360]. 
    GRB file with weather model data can be downloaded from http://rda.ucar.edu/datasets/ds627.0/

    Args:
        * fname   (str)       : Path to the grib file
        * minlat  (float)     : Minimum latitude
        * maxlat  (float)     : Maximum latitude
        * minlon  (float)     : Minimum longitude
        * maxlon  (float)     : Maximum longitude
        * cdic    (float)     : Dictionary of constants
    
    Kwargs:
        * humidity(str)       : Specific ('Q') or relative humidity ('R').

    Returns:
        * lvls    (np.ndarray): Pressure levels
        * latlist (np.ndarray): Latitudes of the stations
        * lonlist (np.ndarray): Longitudes of the stations
        * gph     (np.ndarray): Geopotential height
        * tmp     (np.ndarray): Temperature
        * vpr     (np.ndarray): Vapor pressure

    .. note::
        Uses cc_era by default.
        '''
    
    assert humidity in ('Q','R'), 'Undefined humidity field in get_era.'
    if verbose:
        print('PROGRESS: READING GRIB FILE {}'.format(fname))
    lvls = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 
                     200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 
                     800, 825, 850, 875, 900, 925, 950, 975, 1000])
    nlvls = len(lvls)

    alpha = cdic['Rv']/cdic['Rd']
    gphind = np.arange(nlvls)*12+1
    
    grbs = pygrib.open(fname)
    grbs.seek(gphind[0])
    grb=grbs.read(1)[0]
    lats, lons = grb.latlons()
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
        gph[i,:] = val[ii,jj]/g

        val = grb[1].values           #Reading temperature
        temp = val[ii,jj]
        tmp[i,:] = temp

        if humidity in ('R'):
            esat = cc_era(temp,cdic)       
            grbs.seek(gphind[i]+6)
            grb = grbs.read(1)
            val = grb[0].values
            temp = val[ii,jj]/100.0
            vpr[i,:] = temp*esat
                
        elif humidity in ('Q'):
            val = grb[2].values       #Specific humidity
            temp = val[ii,jj]
            vpr[i,:] = temp*lvls[i]*alpha/(1+(alpha - 1)*temp)
        
        else:
             assert 1==0, 'Undefined Humidity in get_era().' 

    return lvls,latlist,lonlist,gph,tmp,vpr
###############Completed GET_ERA########################################


########Read in ERA data from a given ERA Interim file##################
def get_ecmwf(model,fname,minlat,maxlat,minlon,maxlon,cdic, humidity='Q',verbose=False):
    '''Read data from ERA Interim, ERA-5 or HRES grib file. Note that Lon values should be between [-180, 180].
    Modified by A. Benoit, January 2019.

    Args:
        * model    (str)       : Model used (ERA5, ERAINT or HRES)
        * fname    (str)       : Path to the grib file
        * minlat   (float)     : Minimum latitude
        * maxlat   (float)     : Maximum latitude
        * minlon   (float)     : Minimum longitude
        * maxlon   (float)     : Maximum longitude
        * cdic     (float)     : Dictionary of constants
    
    Kwargs:
        * humidity (str)       : Specific ('Q') or relative humidity ('R').

    Returns:
        * lvls     (np.ndarray): Pressure levels
        * latlist  (np.ndarray): Latitudes of the stations
        * lonlist  (np.ndarray): Longitudes of the stations
        * gph      (np.ndarray): Geopotential height
        * tmp      (np.ndarray): Temperature
        * vpr      (np.ndarray): Vapor pressure

    .. note::
        Uses cc_era by default.
    '''

    assert humidity in ('Q','R'), 'Undefined humidity field in get_era.'
    assert model in ('ERA5', 'ERAINT','HRES'), 'Model not recognized.'
    if verbose:
        print('PROGRESS: READING GRIB FILE')
    if model in 'HRES':
        if verbose:
            print('INFO: USING PRESSURE LEVELS OF HRES DATA')
        lvls = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 150, 
                         200, 250, 300, 400, 500, 600, 700,
                         800, 850, 900, 925, 950, 1000])
    else:
        if verbose:
            print('INFO: USING PRESSURE LEVELS OF ERA-INT OR ERA-5 DATA')
        lvls = np.array([1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 
                         200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775,
                         800, 825, 850, 875, 900, 925, 950, 975, 1000])
    nlvls = len(lvls)

    alpha = cdic['Rv']/cdic['Rd']
    gphind = np.arange(nlvls)*3

    grbs = pygrib.open(fname)
    grbs.seek(gphind[0])
    grb = grbs.read(1)[0]
    lats, lons = grb.latlons()
    #if model == 'ERA5':
    #    lons[lons < 0.] += 360.
    g = cdic['g']

    #extract indices 
    mask = ((lats > minlat) & (lats < maxlat)) & ((lons > minlon) & (lons < maxlon))
    uu = [i for i in list(range(np.shape(mask)[0])) if any(mask[i,:])]
    vv = [j for j in list(range(np.shape(mask)[1])) if any(mask[:,j])]
    
    latlist = lats[uu,:][:,vv]
    lonlist = lons[uu,:][:,vv]
    nlat, nlon = latlist.shape

    ####Create arrays for 3D storage
    gph = np.zeros((nlvls, nlat, nlon))   #Potential height
    tmp = gph.copy()                  #Temperature
    vpr = gph.copy()                  #Vapor pressure
    if verbose:
        print('INFO: IMAGE DIMENSIONS: {} LATITUDES AND {} LONGITUDES'.format(nlat, nlon))

    lvls = 100.0*lvls              #Conversion to absolute pressure
    for i in range(nlvls):
        grbs.seek(gphind[i])   #Reading potential height.
        grb = grbs.read(3)
        gph[i,:,:] = grb[0].values[uu,:][:,vv]/g

        #Reading temperature
        temp = grb[1].values[uu,:][:,vv]
        tmp[i,:,:] = temp

        if humidity in ('R'):   # Relative humidity
            esat = cc_era(temp,cdic)       
            temp = grb[2].values[uu,:][:,vv]/100.0
            vpr[i,:,:] = temp*esat
            
        elif humidity in ('Q'):
            val = grb[2].values  #Specific humidity
            temp = grb[2].values[uu,:][:,vv]
            vpr[i,:,:] = temp*lvls[i]*alpha/(1+(alpha - 1)*temp)
            
        else:
            assert 1==0, 'Undefined Humidity in get_ecmwf().'     

    return lvls,latlist,lonlist,gph,tmp,vpr
###############Completed GET_ECMWF########################################
