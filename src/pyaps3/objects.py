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
import scipy.interpolate as si
from pyaps3 import utils, processor


##############Creating a class object for PyAPS use.
class PyAPS:
    '''Class for dealing with Atmospheric phase corrections in radar/geo coordinates.
    Operates on one weather model file and one Geo.rsc file at a time.'''

    def __init__(self, gribfile, dem, lat, lon, inc=0.0, mask=None,
                 grib='era5', humidity='Q', Del='comb', model='ERA5', verb=False):
        '''Initiates the data structure for atmos corrections in geocoded domain.
        Args:
            * gribfile (str)        : path to downloaded grib file
            * dem      (np.ndarray) : height    in size of (length, width)
            * lat      (np.ndarray) : latitude  in size of (length, width)
            * lon      (np.ndarray) : longitude in size of (length, width)

        Kwargs:
            * inc      (np.ndarray) : incidence angle (in size of (length, width) for np.ndarray)
            * mask     (np.ndarray) : mask of valid pixels in size of (length, width)
            * grib     (str)        : grib name in ['ERA5', 'ERAINT', 'HRES', 'NARR', 'MERRA']
            * humidity (str)        : ['Q', 'R']
            * Del      (str)        : ['comb', 'wet', 'dry']
            * model    (str)        : ECMWF dataset name in ['era5', 'eraint', 'hres']
            * verb     (bool)       : True or False

        .. note::
            For ISCE products, lat/lon can be read from lat/lon.rdr file
            For ROIPAC products, lat, lon = utils.get_lat_lon('radar_16rlks.hgt.rsc')
        '''

        #--------- Check files exist and we have what's needed for computation
        # check grib type and import module
        grib = grib.upper()
        if grib in ['ERA5','ERAINT','HRES']:
            from pyaps3 import era
        elif grib == 'NARR':
            from pyaps3 import narr
        elif grib == 'MERRA':
            from pyaps3 import merra
        else:
            raise ValueError('PyAPS: Undefined grib file source: {}'.format(grib))
        self.grib = grib

        # Check Humidity variable
        humidity = humidity.upper()
        assert humidity in ('Q','R'), 'PyAPS: Undefined humidity.'
        self.humidity = humidity
        if self.grib in ('NARR','MERRA'):
            assert self.humidity in ('Q'), 'PyAPS: Relative humidity not provided by NARR/MERRA.'

        # Check the model for ECMWF
        self.model = model

        # Check grib file exists
        assert os.path.isfile(gribfile), 'PyAPS: GRIB File does not exist.'
        self.gfile = gribfile

        # Get altitude, lon, lat, etc
        self.dem = dem
        self.lon = lon
        self.lat = lat
        self.inc = inc
        self.mask = np.ones(self.dem.shape) if mask is None else mask

        # Get size
        self.ny, self.nx = self.dem.shape
        assert self.lon.shape  == (self.ny, self.nx), 'PyAPS: Longitude array size mismatch'
        assert self.lat.shape  == (self.ny, self.nx), 'PyAPS: Latitude array size mismatch'
        assert self.mask.shape == (self.ny, self.nx), 'PyAPS: Mask array size mismatch'

        # check incidence angle size and type
        if isinstance(self.inc, np.ndarray):
            assert self.inc.shape == (self.ny, self.nx), 'PyAPS: Incidence array size mismatch'
            if verb:
                print('INFO: INCIDENCE ANGLE AS AN ARRAY')
        elif isinstance(self.inc, (int, float, np.float32, np.float64)):
            if verb:
                print('INFO: INCIDENCE ANGLE AS A NUMBER: {} DEG'.format(self.inc))
        else:
            raise ValueError('PyAPS: unrecognized incidence data type: {}'.format(type(self.inc)))

        #--------- Initialize variables
        self.dict = processor.initconst()

        # Get some scales
        if self.grib in ('ERA5','ERAINT','HRES'):
            self.hgtscale = ((self.dict['maxAlt'] - self.dict['minAlt']) / self.dict['nhgt']) / 0.703
            self.bufspc = 1.2
        elif self.grib in ('NARR'):
            self.hgtscale = ((self.dict['maxAlt'] - self.dict['minAlt']) / self.dict['nhgt']) / 0.3
            self.bufspc = 1.2
        elif self.grib in ('MERRA'):
            self.hgtscale = ((self.dict['maxAlt'] - self.dict['minAlt']) / self.dict['nhgt']) / 0.5
            self.bufspc = 1.0 

        # Problems in isce when lon and lat arrays have weird numbers
        if self.grib in ('ERA5','ERAINT','HRES'):
            self.lon[self.lon > 180.] -= 360.
        else:
            self.lon[self.lon < 0.] += 360.0
        self.minlon = np.nanmin(self.lon[np.nonzero(self.mask)]) - self.bufspc
        self.maxlon = np.nanmax(self.lon[np.nonzero(self.mask)]) + self.bufspc
        self.minlat = np.nanmin(self.lat[np.nonzero(self.mask)]) - self.bufspc
        self.maxlat = np.nanmax(self.lat[np.nonzero(self.mask)]) + self.bufspc
        if verb:
            print('INFO: AREA COVERAGE IN SNWE: ({:.2f}, {:.2f}, {:.2f}, {:.2f})'.format(
                self.maxlat, self.minlat, self.minlon, self.maxlon))

        #--------- Extract infos from gribfiles
        if self.grib in ('ERA'):
            assert False, 'Need to modify get_era to fit with the new standards'
            [lvls,latlist,lonlist,gph,tmp,vpr] = era.get_era(
                self.gfile,
                self.minlat,
                self.maxlat,
                self.minlon,
                self.maxlon,
                self.dict,
                humidity=self.humidity,
                verbose=verb,
            )

        elif self.grib in ('ERA5','ERAINT','HRES'):
            [lvls,latlist,lonlist,gph,tmp,vpr] = era.get_ecmwf(
                self.model,
                self.gfile,
                self.minlat,
                self.maxlat,
                self.minlon,
                self.maxlon,
                self.dict,
                humidity=self.humidity,
                verbose=verb,
            )

        elif self.grib in ('NARR'):
            assert False, 'Need to modify get_narr to fit with the new standards'
            [lvls,latlist,lonlist,gph,tmp,vpr] = narr.get_narr(
                self.gfile,
                self.minlat,
                self.maxlat,
                self.minlon,
                self.maxlon,
                self.dict,
                verbose=verb,
            )

        elif self.grib in ('MERRA'):
            assert False, 'Need to modify get_merra to fit with the new standards'
            [lvls,latlist,lonlist,gph,tmp,vpr] = merra.get_merra(
                self.gfile,
                self.minlat,
                self.maxlat,
                self.minlon,
                self.maxlon,
                self.dict,
                verbose=verb,
            )
            lonlist[lonlist < 0.] += 360.0

        # Make a height scale
        hgt = np.linspace(self.dict['minAltP'], gph.max().round(), self.dict['nhgt'])

        # Interpolate pressure, temperature and Humidity of hgt
        [Pi,Ti,Vi] = processor.intP2H(lvls, hgt, gph, tmp, vpr, self.dict, verbose=verb)

        # Calculate the delays
        [DDry,DWet] = processor.PTV2del(Pi,Ti,Vi,hgt,self.dict,verbose=verb)
        if Del.lower() == 'comb':
            Delfn = DDry+DWet
        elif Del.lower() == 'dry':
            Delfn = DDry
        elif Del.lower() == 'wet':
            Delfn = DWet
        else:
            raise ValueError('Unrecognized delay type: {}'.format(Del))

        if self.grib in ('NARR'):
            assert False, 'Need to check narr.intdel'
            [Delfn, latlist, lonlist] = narr.intdel(hgt, latlist, lonlist, Delfn)

        #--------- Save things
        self.Delfn = Delfn
        self.latlist = latlist
        self.lonlist = lonlist
        self.lat = lat
        self.lon = lon
        self.hgt = hgt
        self.Pi = Pi
        self.Ti = Ti
        self.Vi = Vi
        self.verb = verb

        # All done
        return


    def getdelay(self, dout=None, outFile=None, wvl=np.pi*4., writeStations=True):
        '''Get the 2D matrix of tropospheric delay with bilinear interpolation.
        Kwargs:
            * dout    : 2D np.ndarray, output delay matrix
            * outFile : str, file path of output delay matrix
            * wvl     : Wavelength in meters.
                        4*pi (by default) --> output results in delay in meters.
                        0.056 --> output results in delay in radians for C-band SAR.
        Returns:
            * dout    : 2D np.ndarray in size of (ny, nx) in float32.
        '''

        # To know the type of incidence (float, array)
        if isinstance(self.inc, (int, float, np.float32, np.float64)):
            cinc = np.cos(self.inc*np.pi/180.)
            incFileFlag = 'number'
        else:
            incFileFlag = 'array'

        # Get some info from the dictionary
        minAltp = self.dict['minAltP']

        # initiate output
        dout = np.zeros((self.ny, self.nx), dtype=np.float32) if dout is None else dout
        fout = open(outFile, 'wb') if outFile else None

        #######################################################################################
        # BILINEAR INTERPOLATION

        # Create the 1d interpolator to interpolate delays in altitude direction
        if self.verb:
            print('PROGRESS: FINE INTERPOLATION OF HEIGHT LEVELS')
        intp_1d = si.interp1d(self.hgt, self.Delfn, kind='cubic', axis=-1)

        # Interpolate the delay function every meter, for each station
        self.dem[np.isnan(self.dem)] = minAltp
        self.dem[self.dem < minAltp] = minAltp
        minH = np.max([np.nanmin(self.dem*self.mask), self.hgt.min()])
        maxH = int(np.nanmax(self.dem*self.mask)) + 100.
        kh = np.arange(minH,maxH)
        self.Delfn_1m = intp_1d(kh)
        self.alti = kh

        # no reshape
        Lonu = self.lonlist[0,:]
        Latu = self.latlist[:,0]

        # Create the cube interpolator for the bilinear method, to interpolate delays into a grid (x,y,z)
        if self.verb:
            print('PROGRESS: CREATE THE BILINEAR INTERPOLATION FUNCTION')

        # Define a linear interpolating function on the 3D grid: ((x, y, z), data)
        # We do the weird trick of [::-1,:,:] because Latu has to be in increasing order 
        # for the RegularGridInterpolator method of scipy.interpolate
        linearint = si.RegularGridInterpolator(
            (Latu[::-1], Lonu,kh),
            self.Delfn_1m[::-1,:,:],
            method='linear',
            bounds_error=False,
            fill_value=0.0,
        )

        # Show progress bar
        if self.verb:
            toto = utils.ProgressBar(maxValue=self.ny)
            print('PROGRESS: MAPPING THE DELAY')

        # Loop on the lines
        for m in range(self.ny):

            # Update progress bar
            if self.verb:
                toto.update(m+1, every=5)

            # Get latitude and longitude arrays
            lati = self.lat[m,:]*self.mask[m,:]
            loni = self.lon[m,:]*self.mask[m,:]

            # Remove negative values
            if self.grib in ('ERA5','ERAINT','HRES'):
                loni[loni > 180.] -= 360.
            else:
                loni[loni < 0.] += 360.

            # Remove NaN values
            ii = np.where(np.isnan(lati))
            jj = np.where(np.isnan(loni))
            xx = np.union1d(ii,jj)
            lati[xx]=0.0
            loni[xx]=0.0

            # Get incidence if file provided
            if incFileFlag == 'array':
                cinc = np.cos(self.mask[m,:]*self.inc[m,:]*np.pi/180.)

            # Make the bilinear interpolation
            D = self.dem[m,:]
            val = linearint(np.vstack((lati, loni, D)).T)*np.pi*4.0/(cinc*wvl)
            val[xx] = np.nan

            # save output
            dout[m,:] = val
            if outFile:
                val.astype(np.float32).tofile(fout)

        if self.verb:
            toto.close()

        if outFile:
            fout.close()

        return dout
