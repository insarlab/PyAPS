############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################


import os.path
import sys
import numpy as np
import scipy.integrate as intg
import scipy.interpolate as si
import scipy.spatial as ss

import pyaps3.processor as processor
import pyaps3.utils as utils


def sanity_check(model):
    if model in ('ECMWF','ERA','NARR'):
        from . import era, narr
    if model=='MERRA':
        from . import merra
    return


##############Creating a class object for PyAPS use.
class PyAPS_geo:
    '''
    Class for dealing with Atmospheric phase corrections in geo-coded space.
    Operates on one weather model file and one Geo.rsc file at a time.

    This use the relative humidity information from the specified ERA file over the area specified by the dem rsc file.
    '''

    def __init__(self,gribfile,demfile,grib='ECMWF',humidity='Q',demtype=np.int16,demfmt='HGT',verb=False, Del='comb'):
        '''
        Initiates the data structure for atmos corrections in geocoded domain.
        Args:
            * gribfile (str)   : Path to the weather model file.
            * demfile  (str)   : Path to the SIM_nRLKS.hgt file.

        Kwargs:
            * grib     (str)   : Can be ECMWF, ERA, NARR or MERRA depending on the model.
            * humidity (str)   : Can be 'Q' or 'R' depdending on specific or relative humidity.
            * demtype  (dtype) : Data type of the DEM File.
            * demfmt   (str)   : Can be 'HGT' or 'RMG' depending on type of file.
                                 If RMG, altitude is assumed to be second channel.
            * Del      (str)   : Can be 'comb', 'dry' or 'wet'. Specifies which delay you want (both, dry or wet)

        .. note :: 
        The DEMFile is only used to set up the geographic boundaries of the problem and the size of the output matrix
        during init. The GRIBfile is read and the 3D interpolation function is prepared.
        '''

        sanity_check(grib)
        grib = grib.upper()
        humidity = humidity.upper()
        demfmt = demfmt.upper()

        assert grib in ('ERA','NARR','ECMWF','MERRA'), 'PyAPS: Undefined grib file source.'
        self.grib = grib
        '''String with the name of the weather model.'''

        assert humidity in ('Q','R'), 'PyAPS: Undefined humidity.'
        self.humidity = humidity

        if self.grib in ('NARR','MERRA'):
            assert self.humidity in ('Q'), 'PyAPS: Relative humidity not provided by NARR/MERRA.'

        assert os.path.isfile(gribfile), 'PyAPS: GRIB File does not exist.'
        self.gfile = gribfile

        assert os.path.isfile(demfile), 'PyAPS: DEM file does not exist.'
        self.hfile = demfile


        assert demfmt in ('RMG','HGT','VAR'), 'PyAPS: DEM Format can be RMG or HGT, VAR'
        self.fmt = demfmt

        self.demtype = demtype
        self.dict = processor.initconst()
        self.bufspc = 1.2           ####Hardcoded.

        if grib in ('ERA','ECMWF'):
            self.hgtscale = ((self.dict['maxAlt']-self.dict['minAlt'])/self.dict['nhgt'])/0.703
        elif grib in ('NARR'):
            self.hgtscale = ((self.dict['maxAlt']-self.dict['minAlt'])/self.dict['nhgt'])/0.3
        elif grib in ('MERRA'):
            self.hgtscale = ((self.dict['maxAlt']-self.dict['minAlt'])/self.dict['nhgt'])/0.5

        [lon,lat,nx,ny] = utils.geo_rsc(self.hfile,verbose=verb)


        ######Problems near the international dateline
        if grib in ('MERRA'):
            #MERRA's longtidue range -180~180.
            #This is to correct the value from utils.geo_rsc() function. Cunren, OCT-2015
            lon[lon > 180] -= 360.0
            #if(lon[0] > 180):
            #    lon[0] = lon[0] - 360.0
            #if(lon[1] > 180):
            #    lon[1] = lon[1] - 360.0
    
            self.minlon = lon.min()-self.bufspc
            self.maxlon = lon.max()+self.bufspc
            self.minlat = lat.min()-self.bufspc
            self.maxlat = lat.max()+self.bufspc
            self.nx = nx
            self.ny = ny
        else:
            lon[lon < 0] += 360.0
            self.minlon = lon.min()-self.bufspc
            self.maxlon = lon.max()+self.bufspc
            self.minlat = lat.min()-self.bufspc
            self.maxlat = lat.max()+self.bufspc
            self.nx = nx
            self.ny = ny

        if self.grib in ('ERA'):
            [lvls,latlist,lonlist,gph,tmp,vpr] = era.get_era(self.gfile,
                                                             self.minlat,
                                                             self.maxlat,
                                                             self.minlon,
                                                             self.maxlon,
                                                             self.dict,
                                                             humidity=self.humidity,
                                                             verbose=verb)

        elif self.grib in ('ECMWF'):
            [lvls,latlist,lonlist,gph,tmp,vpr] = era.get_ecmwf(self.gfile,
                                                               self.minlat,
                                                               self.maxlat,
                                                               self.minlon,
                                                               self.maxlon,
                                                               self.dict,
                                                               humidity=self.humidity,
                                                               verbose=verb)

        elif self.grib in ('NARR'):
            [lvls,latlist,lonlist,gph,tmp,vpr] = narr.get_narr(self.gfile,
                                                               self.minlat,
                                                               self.maxlat,
                                                               self.minlon,
                                                               self.maxlon,
                                                               self.dict,
                                                               verbose=verb)
        elif self.grib in ('MERRA'):
            [lvls,latlist,lonlist,gph,tmp,vpr] = merra.get_merra(self.gfile,
                                                                 self.minlat,
                                                                 self.maxlat,
                                                                 self.minlon,
                                                                 self.maxlon,
                                                                 self.dict,
                                                                 verbose=verb)
            lonlist[lonlist < 0.] += 360.0

        hgt = np.linspace(self.dict['minAlt'],self.dict['maxAlt'],self.dict['nhgt'])

        [Pi,Ti,Vi] = processor.intP2H(lvls,hgt,gph,tmp,vpr,self.dict,verbose=verb)

        [DDry,DWet] = processor.PTV2del(Pi,Ti,Vi,hgt,self.dict,verbose=verb)

        if Del in ('comb','Comb'):
            Delfn=DDry+DWet
        elif Del in ('dry','Dry'):
            Delfn=DDry
        elif Del in ('wet','Wet'):
            Delfn=DWet
        else:
            print('Unrecognized delay type')
            sys.exit(1)

        if self.grib in ('NARR'):
            [Delfn,latlist,lonlist] = narr.intdel(hgt,latlist,lonlist,Delfn)

        self.Delfn = Delfn
        self.lonlist = lonlist
        self.latlist = latlist
        self.Pi = Pi
        self.Vi = Vi
        self.Ti = Ti
        self.hgt = hgt
        self.verb = verb


    def merisfactor(self,dataobj,inc=0.0,wvl=4*np.pi):
            '''
            Write pi-factor from Li et al 2012 to a matrix / HDF5 object or a file directly.

             Args:
                     * dataobj  (str or HDF5 or np.array): Final output. If str, output is written to file.

             Kwargs:
                     * inc  (np.float): Incidence angle in degrees. Default is vertical.
                     * wvl  (np.float): Wavelength in meters. Default output results in delay in meters.

             .. note::
                     If dataobj is string, output is written to the file.
                     If np.array or HDF5 object, it should be of size (ny,nx).
            '''

            minAltp = self.dict['minAltP']

            # Incidence
            cinc = np.cos(inc*np.pi/180.0)

            # Compute the two integrals
            WonT = self.Vi/self.Ti
            WonT2 = WonT/self.Ti

            S1 = intg.cumtrapz(WonT,x=self.hgt,axis=-1)
            val = 2*S1[:,-1]-S1[:,-2]
            val = val[:,None]
            S1 = np.concatenate((S1,val),axis=-1)
            del WonT
            S2 = intg.cumtrapz(WonT2,x=self.hgt,axis=-1)
            val = 2*S2[:,-1]-S2[:,-2]
            val = val[:,None]
            S2 = np.concatenate((S2,val),axis=-1)
            del WonT2
            Tm = S1/S2
            self.Tm = Tm

            # Reading in the DEM
            if self.verb:
                print('PROGRESS: READING DEM')
            fin = open(self.hfile,'rb');
            if self.fmt in ('HGT'):
                dem = np.fromfile(file=fin,dtype=self.demtype,count=self.nx*self.ny).reshape(self.ny,self.nx)
            elif self.fmt in ('RMG'):
                dem = np.fromfile(file=fin,dtype=self.demtype,count=2*self.nx*self.ny).reshape(self.ny,2*self.nx)
                dem = dem[:,self.nx:]
            dem = np.round(dem).astype(np.int)
            fin.close()

            # check output, and open file if necessary
            outFile = isinstance(dataobj,str)
            if outFile:
                fout = open(dataobj,'wb')
                dout = np.zeros((self.ny,self.nx))
            else:
                assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), 'PyAPS: Not a valid data object.'
                dout = dataobj

            # Create the lon/lat arrays
            laty = np.linspace(self.maxlat-self.bufspc,self.minlat+self.bufspc,self.ny)
            lonx = np.linspace(self.minlon+self.bufspc,self.maxlon-self.bufspc,self.nx)

            # Create the 1d interpolator
            if self.verb:
                print('PROGRESS: FINE INTERPOLATION OF HEIGHT LEVELS')
            intp_1d = si.interp1d(self.hgt,Tm,kind='cubic',axis=1)

            # Interpolate the Tm variable every meter
            dem[dem < minAltp] = minAltp
            minH = dem.min()
            maxH = dem.max()+1
            kh = np.arange(minH,maxH)
            Tm_1m = intp_1d(kh)
            self.alti = kh

            # Reshape Tm
            Lonu = np.unique(self.lonlist)
            Latu = np.unique(self.latlist)
            nLon = len(Lonu)
            nLat = len(Latu)
            Tm_1m = np.reshape(Tm_1m.T,(len(kh),nLat,nLon))
            self.Tm_1m = Tm_1m

            # Create the cube interpolator for the bilinear method
            if self.verb:
                print('PROGRESS: CREATE THE BILINEAR INTERPOLATION FUNCTION')
            bilicube = processor.Bilinear2DInterpolator(Lonu,Latu,Tm_1m,cube=True)

            # Get the values from the dictionnary
            k1 = self.dict['k1']
            k2 = self.dict['k2']
            k3 = self.dict['k3']
            mmO = self.dict['mmO']
            mmH = self.dict['mmH']
            mma = self.dict['mma']
            w = (2*mmH + mmO)/mma
            Rv = self.dict['Rv']
            Rho = self.dict['Rho']

            # Loop on the lines
            if self.verb:
                toto = utils.ProgressBar(maxValue=self.ny)
                print('PROGRESS: MAPPING THE DELAY')

            for m in range(self.ny):
                if self.verb:
                    toto.update(m,every=5)

                # Get Lon/Lat
                loni = lonx
                lati = laty[m]*np.ones((loni.shape))

                # Make the bilinear interpolation
                D = dem[m,:] - minH
                val = bilicube(loni,lati,D)
                val = 0.000001 * Rho * Rv * ( k3/val + k2 - w*k1) * np.pi*4.0/(cinc*wvl)

                if outFile:
                    resy = val.astype(np.float32)
                    resy.tofile(fout)
                else:
                    dataobj[m,:] = val
                if self.verb:
                    toto.close()

            # Close if outfile              
            if outFile:
                fout.close()

###########End of PyAPS_geo class###############################
