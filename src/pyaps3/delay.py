############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
import numpy as np

#############Clausius-Clapeyron for ECMWF as used in Jolivet et al 2011#
def cc_eraorig(tmp,cdic):
    '''This routine takes temperature profiles and returns Saturation water vapor
    partial pressure using the Clausius-Clapeyron law applied in Jolivet et al. 2011,
    GRL, doi:10.1029/2011GL048757. It can be used in case you are using Relative.
    Humidity from ECMWF models.

    Args:
        * tmp  (np.ndarray) : Temperature
        * cdic (dict)       : Dictionnary of constants

    Returns:
        * esat (np.ndarray) : Saturation water vapor partial pressure.
    '''

    a1w=cdic['a1w']
    T3=cdic['T3']
    Rv=cdic['Rv']

    esat=a1w*np.exp( (2.5e6/Rv) * ( (1/T3) - (1/tmp) ) )

    return esat
###############Completed CC_ERAORIG#####################################


##############Read Input text file #####################################
def read_eratxt(fname,cdic):
    '''Read ECMWF files from 0.75 degree grid similar to Romain Jolivet's Delay Package.

    Args:
        * fname    (str)       :  Path to the delay file
        * cdic     (float)     :  Dictionary of constants

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
        Uses cc_eraorig by default.
    '''

    lvls=[]
    latlist=[]
    lonlist=[]
    gpht=[]
    tmpt=[]
    reht=[]

    g=cdic['g']

    f=open(fname,'r')
    tmp=f.readlines()
    i=0
    nstn=0
    maxloop=int(len(tmp))
    while i<maxloop:
        if (tmp[i][0]=='-'):
            nstn=nstn+1
            lonlat=tmp[i+3].rsplit(' ')
            lonlist.append(float(lonlat[3]))
            latlist.append(float(lonlat[9]))
            i=i+5
            new='y'
        else:
            if new in ('y'):
                n=1
                val=tmp[i].rsplit(' ')
                gpht.append(float(val[0]))
                lvls.append(float(val[1]))
                tmpt.append(float(val[2]))
                reht.append(float(val[3]))
                i=i+1
                new='n'
            else:
                n=n+1
                val=tmp[i].rsplit(' ')
                gpht.append(float(val[0]))
                lvls.append(float(val[1]))
                tmpt.append(float(val[2]))
                reht.append(float(val[3]))
                i=i+1

    gpht=np.array(gpht)/g
    gph=np.flipud(gpht.reshape((n,nstn),order='F'))
    del gpht

    tmpt=np.array(tmpt)
    esat=cc_eraorig(tmpt,cdic)
    tmp=np.flipud(tmpt.reshape((n,nstn),order='F'))
    del tmpt

    vprt=(np.array(reht)/100.)*esat
    vpr=np.flipud(vprt.reshape((n,nstn),order='F'))
    del vprt
    del esat

    lvls=np.flipud(np.array(lvls))
    lvls=lvls[0:n]

    lonlist=np.array(lonlist)
    latlist=np.array(latlist)

    return lvls,latlist,lonlist,gph,tmp,vpr
