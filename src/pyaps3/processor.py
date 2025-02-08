############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################


import numpy as np
import scipy.interpolate as intp
import scipy.integrate as intg


def initconst():
    '''Initialization of various constants needed for computing delay.

    Args:
        * None

    Returns:
        * constdict (dict): Dictionary of constants
    '''

    constdict = {}
    constdict['k1'] = 0.776   #(K/Pa)
    constdict['k2'] = 0.716   #(K/Pa)
    constdict['k3'] = 3750    #(K^2.Pa)
    constdict['g'] = 9.81     #(m/s^2)
    constdict['Rd'] = 287.05  #(J/Kg/K)
    constdict['Rv'] = 461.495 #(J/Kg/K)
    constdict['mma'] = 29.97  #(g/mol)
    constdict['mmH'] = 2.0158 #(g/mol)
    constdict['mmO'] = 16.0   #(g/mol)
    constdict['Rho'] = 1000.0 #(kg/m^3)

    constdict['a1w'] = 611.21 # hPa
    constdict['a3w'] = 17.502 #
    constdict['a4w'] = 32.19  # K
    constdict['a1i'] = 611.21 # hPa
    constdict['a3i'] = 22.587 #
    constdict['a4i'] = -0.7   # K
    constdict['T3'] = 273.16  # K
    constdict['Ti'] = 250.16  # K
    constdict['nhgt'] = 300   # Number of levels for interpolation (OLD 151)
    constdict['minAlt'] = -200.0
    constdict['maxAlt'] = 50000.0
    constdict['minAltP'] = -200.0
    return constdict    

###############Completed the list of constants################


##########Interpolating to heights from Pressure levels###########
def intP2H(lvls,hgt,gph,tmp,vpr,cdic,verbose=False):
    '''Interpolates the pressure level data to altitude.

    Args:
        * lvls (np.ndarray) : Pressure levels.
        * hgt  (np.ndarray) : Height values for interpolation.
        * gph  (np.ndarray) : Geopotential height.
        * tmp  (np.ndarray) : Temperature.
        * vpr  (np.ndarray) : Vapor pressure.
        * cdic (dict)       : Dictionary of constants.

    .. note::
        gph,tmp,vpr are of size (nstn,nlvls).

    Returns:
        * Presi (np.ndarray): Interpolated pressure.
        * Tempi (np.ndarray): Interpolated temperature.
        * Vpri  (np.ndarray): Interpolated vapor pressure.

    .. note::
        Cubic splines are used to convert pressure level data to height level data.'''

    minAlt = cdic['minAlt']      #Hardcoded parameter.
    maxAlt = gph.max().round()

    if verbose:
        print('PROGRESS: INTERPOLATING FROM PRESSURE TO HEIGHT LEVELS')
    nlat = gph.shape[1]           #Number of stations
    nlon = gph.shape[2]
    nhgt = len(hgt)               #Number of height points
    Presi = np.zeros((nlat,nlon,nhgt))
    Tempi = np.zeros((nlat,nlon,nhgt))
    Vpri  = np.zeros((nlat,nlon,nhgt))

    for i in range(nlat):
        for j in range(nlon):
            temp = gph[:,i,j]       #Obtaining height values
            hx = temp.copy()
            sFlag = False
            eFlag = False
            #Add point at start
            if (hx.min() > minAlt):
                sFlag = True
                #changed from 1 to 0 (-1 should also work), CL
                hx = np.concatenate((hx,[minAlt-1]),axis = 0)

            #Add point at end
            if (hx.max() < maxAlt):
                eFlag = True
                #changed from 1 to 0 (-1 should also work), CL
                hx = np.concatenate(([maxAlt+1],hx),axis=0)

            #Splines needs monotonically increasing.
            hx = -hx

            #Interpolating pressure values
            hy = lvls.copy()
            if (sFlag == True):
                val = hy[-1] + (hx[-1] - hx[-2]) * (hy[-1] - hy[-2]) / (hx[-2] - hx[-3])
                #changed from 1 to 0 (-1 should also work), CL
                hy = np.concatenate((hy,[val]),axis=0)

            if (eFlag == True):
                val = hy[0] - (hx[0] - hx[1]) * (hy[0] - hy[1]) / (hx[1] - hx[2])
                #changed from 1 to 0 (-1 should also work), CL
                hy = np.concatenate(([val],hy),axis=0)

            tck = intp.interp1d(hx,hy,kind='cubic')
            #Again negative for consistency with hx
            temp = tck(-hgt)
            Presi[i,j,:] = temp.copy()
            del temp

            #Interpolating temperature
            temp = tmp[:,i,j]
            hy = temp.copy()
            if (sFlag == True):
                val = hy[-1] + (hx[-1] - hx[-2]) * (hy[-1] - hy[-2]) / (hx[-2] - hx[-3])
                #changed from 1 to 0 (-1 should also work), CL
                hy = np.concatenate((hy,[val]),axis=0)

            if (eFlag == True):
                val = hy[0] - (hx[0] - hx[1]) * (hy[0] - hy[1]) / (hx[1] - hx[2])
                #changed from 1 to 0 (-1 should also work), CL
                hy = np.concatenate(([val],hy),axis=0)

            tck = intp.interp1d(hx,hy,kind='cubic')
            temp = tck(-hgt)
            Tempi[i,j,:] = temp.copy()
            del temp

            #Interpolating vapor pressure
            temp = vpr[:,i,j]
            hy = temp.copy()
            if (sFlag == True):
                val = hy[-1] + (hx[-1] - hx[-2]) * (hy[-1] - hy[-2]) / (hx[-2] - hx[-3])
                #changed from 1 to 0 (-1 should also work), CL
                hy = np.concatenate((hy,[val]),axis=0)

            if (eFlag == True):
                val = hy[0] - (hx[0] - hx[1]) * (hy[0] - hy[1]) / (hx[1] - hx[2])
                #changed from 1 to 0 (-1 should also work), CL
                hy = np.concatenate(([val],hy),axis=0)

            tck = intp.interp1d(hx,hy,kind='cubic')
            temp = tck(-hgt)
            Vpri[i,j,:] = temp.copy()
            del temp

    ## Save into files for plotting (test for era5 interpolation)
    #out_dir = '/home/angel/Tests/test_era5interpolation/Outputs'
    #suffix = '20141023_20141210_era5'
    #np.savetxt(f'{out_dir}/alt_{suffix}.txt',   gph,     fmt=' '.join(['%s']*380)) # altitude (of lvls)
    #np.savetxt(f'{out_dir}/lvls_{suffix}.txt',  lvls.T,  fmt='%s')                 # pressure
    #np.savetxt(f'{out_dir}/tmp_{suffix}.txt',   tmp,     fmt=' '.join(['%s']*380)) # temperature
    #np.savetxt(f'{out_dir}/vpr_{suffix}.txt',   vpr,     fmt=' '.join(['%s']*380)) # vapor
    #np.savetxt(f'{out_dir}/Alti_{suffix}.txt',  hgt.T,   fmt='%s')                 # interpo altitude
    #np.savetxt(f'{out_dir}/Presi_{suffix}.txt', Presi.T, fmt=' '.join(['%s']*380)) # interpo pressure
    #np.savetxt(f'{out_dir}/Tempi_{suffix}.txt', Tempi.T, fmt=' '.join(['%s']*380)) # interpo temperature
    #np.savetxt(f'{out_dir}/Vpri_{suffix}.txt',  Vpri.T,  fmt=' '.join(['%s']*380)) # interpo vapor

    return Presi,Tempi,Vpri
###########Completed interpolation to height levels #####################


###########Computing the delay function ###############################
def PTV2del(Presi,Tempi,Vpri,hgt,cdict,verbose=False):
    '''Computes the delay function given Pressure, Temperature and Vapor pressure.

    Args:
        * Presi (np.ndarray) : Pressure at height levels.
        * Tempi (np.ndarray) : Temperature at height levels.
        * Vpri  (np.ndarray) : Vapor pressure at height levels.
        * hgt   (np.ndarray) : Height levels.
        * cdict (np.ndarray) : Dictionary of constants.

    Returns:
        * DDry2 (np.ndarray) : Dry component of atmospheric delay.
        * DWet2 (np.ndarray) : Wet component of atmospheric delay.

    .. note::
        Computes refractive index at each altitude and integrates the delay using cumtrapz.'''

    if verbose:
        print('PROGRESS: COMPUTING DELAY FUNCTIONS')
    nhgt = len(hgt)              #Number of height points
    nlat = Presi.shape[0]        #Number of stations
    nlon = Presi.shape[1]
    WonT = Vpri/Tempi
    WonT2 = WonT/Tempi

    k1 = cdict['k1']
    Rd = cdict['Rd']
    Rv = cdict['Rv']
    k2 = cdict['k2']
    k3 = cdict['k3']
    g = cdict['g']

    #Dry delay
    DDry2 = np.zeros((nlat,nlon,nhgt))
    DDry2[:,:,:] = k1*Rd*(Presi[:,:,:] - Presi[:,:,-1][:,:,np.newaxis])*1.0e-6/g

    #Wet delay
    S1 = intg.cumulative_trapezoid(WonT,x=hgt,axis=-1)
    val = 2*S1[:,:,-1]-S1[:,:,-2]
    val = val[:,:,None]
    S1 = np.concatenate((S1,val),axis=-1)
    del WonT

    S2 = intg.cumulative_trapezoid(WonT2,x=hgt,axis=-1)
    val = 2*S2[:,:,-1]-S2[:,:,-2]
    val = val[:,:,None]
    S2 = np.concatenate((S2,val),axis=-1)
    DWet2 = -1.0e-6*((k2-k1*Rd/Rv)*S1+k3*S2)
    
    for i in range(nlat):
        for j in range(nlon):
            DWet2[i,j,:]  = DWet2[i,j,:] - DWet2[i,j,-1]

    return DDry2,DWet2
####################Completed delay function#################


####Setting up 3D interpolation function in geo/xy coordinates#######
def make3dintp(Delfn,lonlist,latlist,hgt,hgtscale):
    '''Returns a 3D interpolation function that can be used to interpolate using llh coordinates.

    Args:
        * Delfn    (np.ndarray) : Array of delay values.
        * lonlist  (np.ndarray) : Array of station longitudes.
        * latlist  (np.ndarray) : Array of station latitudes.
        * hgt      (np.ndarray) : Array of height levels.
        * hgtscale (float)      : Height scale factor for interpolator.

    Returns:
        * fnc  (function) : 3D interpolation function.

    .. note::
        We currently use the LinearNDInterpolator from scipy.
    '''

    ##Delfn   = Ddry + Dwet. Delay function.
    ##lonlist = list of lons for stations. / x
    ##latlist = list of lats for stations. / y
    nstn = Delfn.shape[0]
    nhgt = Delfn.shape[1]
    xyz = np.zeros((nstn*nhgt,3))
    Delfn = np.reshape(Delfn,(nstn*nhgt,1))
    print('3D interpolation')
    count = 0
    for m in range(nstn):
        for n in range(nhgt):
            xyz[count,0] = lonlist[m]
            xyz[count,1] = latlist[m]
            xyz[count,2] = hgt[n]/hgtscale     #For same grid spacing as lat/lon
            count += 1

    #xyz[:,2] = xyz[:,2] #+ 1e-30*np.random.rand((nstn*nhgt))/hgtscale #For unique Delaunay    
    del latlist
    del lonlist
    del hgt
    if verbose:
        print('PROGRESS: BUILDING INTERPOLATION FUNCTION')
    fnc = intp.LinearNDInterpolator(xyz,Delfn)

    return fnc
###########Completed 3D interpolation in geo coordinates############


########### Class for bilinear interpolation at a certain level in a 3d cube
class Bilinear2DInterpolator:
    '''Bilinear interpolation in 2D. The code is modified from mpl_toolkits.basemap.interp and scipy.interpolate'''
    def __init__(self, xin, yin, datain,cube=False):
        '''Setting up the interpolator.

        .. Args:

            * xin     -> Monotonic array of x coordinates
            * yin     -> Monotonic array of y coordinates
            * datain  -> 2D array corresponding to (y,x) if cube=False, else 3D array corresponding to (nz,y,x)
        
        .. Kwargs:
            
            * cube    -> 2D array or 3D array cube.'''

        if cube:
            if xin.shape[0] != datain.shape[2]:
                raise ValueError('Shapes of datain and x do not match')

            if yin.shape[0] != datain.shape[1]:
                raise ValueError('Shapes of datain and y do not match')

        else:
            if xin.shape[0] != datain.shape[1]:
                raise ValueError('Shapes of datain and x do not match')

            if yin.shape[0] != datain.shape[0]:
                raise ValueError('Shapes of datain and y do not match')

        if xin[-1] < xin[0]:
            raise ValueError('Array x not sorted')

        if yin[-1] < yin[0]:
            raise ValueError('Array y not sorted')

        self.xin = xin.copy()
        self.yin = yin.copy()
        
        delx = xin[1:] - xin[0:-1]
        dely = yin[1:] - yin[0:-1]

        if max(delx)-min(delx) < 1.e-4 and max(dely)-min(dely) < 1.e-4:
            self.regular = True
        else:
            self.regular = False

        self.xinlist = self.xin.tolist()
        self.yinlist = self.yin.tolist()
        self.nx = len(self.xinlist)
        self.ny = len(self.yinlist)
        self.cube = cube
        self.zin = datain.copy()

        # All done
        return

    def __call__(self,xi,yi,iz=0):
        '''Function call to actually interpolate.'''
        if xi.shape != yi.shape:
            raise ValueError('xi and yi must have same shape.')

        if self.regular:
            xcoords = (self.nx-1)*(xi-self.xin[0])/(self.xin[-1]-self.xin[0])
            ycoords = (self.ny-1)*(yi-self.yin[0])/(self.yin[-1]-self.yin[0])
        else:
            xiflat = xi.flatten()
            yiflat = yi.flatten()
            ix = (np.searchsorted(self.xin,xiflat)-1).tolist()
            iy = (np.searchsorted(self.yin,yiflat)-1).tolist()
            xiflat = xiflat.tolist()
            yiflat = yiflat.tolist()

            xin = self.xinlist
            yin = self.yinlist
                
            xcoords = []
            ycoords = []
            for n,i in enumerate(ix):
                if i < 0:
                    xcoords.append(-1)
                elif i >= self.nx-1:
                    xcoords.append(self.nx)
                else:
                    xcoords.append(float(i)+(xiflat[n]-xin[i])/(xin[i+1]-xin[i]))
            for m,j in enumerate(iy):
                if j < 0:
                    ycoords.append(-1)
                elif j >= self.ny-1:
                    ycoords.append(self.ny)
                else:
                    ycoords.append(float(j)+(yiflat[m]-yin[j])/(yin[j+1]-yin[j]))

            xcoords = np.reshape(xcoords, xi.shape)
            ycoords = np.reshape(ycoords, yi.shape)

        xcoords = np.clip(xcoords,0,self.nx-1)
        ycoords = np.clip(ycoords,0,self.ny-1)

        xint = xcoords.astype(np.int32)
        yint = ycoords.astype(np.int32)
        xip1 = np.clip(xint+1,0,self.nx-1)
        yip1 = np.clip(yint+1,0,self.ny-1)

        delx = xcoords - xint.astype(np.float32)
        dely = ycoords - yint.astype(np.float32)

        zin = self.zin
        if self.cube:
            dataout = (1.-delx)*(1.-dely)*zin[iz,yint,xint] + \
                      delx*dely*zin[iz,yip1,xip1] + \
                      (1.-delx)*dely*zin[iz,yip1,xint] + \
                      delx*(1.-dely)*zin[iz,yint,xip1]
        else:
            dataout = (1.-delx)*(1.-dely)*zin[yint,xint] + \
                      delx*dely*zin[yip1,xip1] + \
                      (1.-delx)*dely*zin[yip1,xint] + \
                      delx*(1.-dely)*zin[yint,xip1]

        return dataout
