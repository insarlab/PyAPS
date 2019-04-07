############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
import os.path
import numpy as np
import pyaps3.utils as utils
import pyaps3.processor as processor
import scipy.integrate as intg

GRIBflag = True
MERRAflag = True

try:
    import pyaps3.era as era
except:
    GRIBflag = False

try:
    import pyaps3.narr as narr
except:
    GRIBflag = False

try:
    import pyaps3.merra as merra
except:
    MERRAflag = False

import scipy.interpolate as si
import matplotlib.pyplot as plt
import sys

def sanity_check(model):
    if model in ('ECMWF','ERA','NARR') and GRIBflag==False:
        print('No module name pygrib found.')
        print('To use ECMWF and NARR, please install pygrib.')

    if model=='MERRA' and MERRAflag==False:
        print('No module name pyhdf found')
        print('To use MERRA, please install pyhdf')
        print('Visit: http://pysclint.sourceforge.net/pyhdf')


##############Creating a class object for PyAPS use.
class PyAPS_rdr:
    '''Class for dealing with Atmospheric phase corrections in radar coordinates.
    Operates on one weather model file and one Geo.rsc file at a time.'''
    
    def __init__(self,gribfile,
                      dem, lon, lat, inc,
                      grib='ECMWF',
                      humidity='Q',
                      verb=False, 
                      mask=None,
                      box=None,
                      Del='comb',
                      model=None):
        '''Initiates the data structure for atmos corrections in geocoded domain.'''

        # Check files exist and we have what's needed for computation
        sanity_check(grib)
        grib = grib.upper()
        humidity = humidity.upper()
        
        # Check grib type name
        assert grib in ('ERA','NARR','ECMWF','MERRA'), \
                'PyAPS: Undefined grib file source.'
        self.grib = grib

        # Check the model for ECMWF
        self.model = model
        
        # Check Humidity variable
        assert humidity in ('Q','R'), 'PyAPS: Undefined humidity.'
        self.humidity = humidity
        if self.grib in ('NARR','MERRA'):
            assert self.humidity in ('Q'), \
                    'PyAPS: Relative humidity not provided by NARR/MERRA.'

        # Check grib file exists
        assert os.path.isfile(gribfile), 'PyAPS: GRIB File does not exist.'
        self.gfile = gribfile
        
        # Get altitude, lon, lat, etc
        self.dem = dem
        self.lon = lon 
        self.lat = lat
        self.inc = inc
        if mask is None:
            self.mask = np.ones(self.dem.shape)
        else:
            self.mask = mask

        # Get size
        self.ny, self.nx = self.dem.shape
        #self.inc = np.ones(self.dem.shape) * self.inc
        assert self.lon.shape==(self.ny, self.nx), 'Longitude array size mismatch'
        assert self.lat.shape==(self.ny, self.nx), 'Latitude array size mismatch'
        assert self.mask.shape==(self.ny, self.nx), 'Mask array mismatch'

        # Initialize variables
        self.dict = processor.initconst()

        # Get some scales
        if grib in ('ERA','ECMWF'):
            if GRIBflag:
                self.bufspc = 1.2
            else:
                print('================================')
                print('********************************')
                print('  pyGrib needs to be installed  ')
                print('    No ECMWF or NARR possible   ')
                print('********************************')
                print('================================')
                sys.exit(1)
        elif grib in ('NARR'):
            if GRIBflag:
                self.bufspc = 1.2
            else:
                print('================================')
                print('********************************')
                print('  pyGrib needs to be installed  ')
                print('    No ECMWF or NARR possible   ')
                print('********************************')
                print('================================')
                sys.exit(1)
            self.bufspc = 1.2
        elif grib in ('MERRA'):
            if MERRAflag:
                self.bufspc = 1.0 
            else:
                print('================================')
                print('********************************')
                print('  pyHdf needs to be installed   ')
                print('        No MERRA possible       ')
                print('********************************')
                print('================================')
                sys.exit(1)
            
        ######Problems in isce when lon and lat arrays have weird numbers
        self.lon[self.lon < 0.] += 360.0
        if box is None:
            self.minlon = np.nanmin(self.lon*self.mask)-self.bufspc
            self.maxlon = np.nanmax(self.lon*self.mask)+self.bufspc
            self.minlat = np.nanmin(self.lat*self.mask)-self.bufspc
            self.maxlat = np.nanmax(self.lat*self.mask)+self.bufspc
        else:
            print('Box specified')
            self.minlon, self.maxlon, self.minlat, self.maxlat = box
            self.minlon -= self.bufspc
            self.maxlon += self.bufspc
            self.minlat -= self.bufspc
            self.maxlat += self.bufspc

        # Extract infos from gribfiles
        if self.grib in ('ERA'):
            assert False, 'Need to modify get_era to fit with the new standards'
            [lvls,latlist,lonlist,gph,tmp,vpr] = era.get_era(self.gfile,
                                                             self.minlat,
                                                             self.maxlat,
                                                             self.minlon,
                                                             self.maxlon,
                                                             self.dict,
                                                             humidity=self.humidity,
                                                             verbose=verb)
        elif self.grib in ('ECMWF'):
            [lvls,latlist,lonlist,gph,tmp,vpr] = era.get_ecmwf(self.model,
                                                               self.gfile,
                                                               self.minlat,
                                                               self.maxlat,
                                                               self.minlon,
                                                               self.maxlon,
                                                               self.dict,
                                                               humidity=self.humidity,
                                                               verbose=verb)
        elif self.grib in ('NARR'):
            assert False, 'Need to modify get_narr to fit with the new standards'
            [lvls,latlist,lonlist,gph,tmp,vpr] = narr.get_narr(self.gfile,
                                                               self.minlat,
                                                               self.maxlat,
                                                               self.minlon,
                                                               self.maxlon,
                                                               self.dict,
                                                               verbose=verb)
        elif self.grib in ('MERRA'):
            assert False, 'Need to modify get_merra to fit with the new standards'
            [lvls,latlist,lonlist,gph,tmp,vpr] = merra.get_merra(self.gfile,
                                                                 self.minlat,
                                                                 self.maxlat,
                                                                 self.minlon,
                                                                 self.maxlon,
                                                                 self.dict,
                                                                 verbose=verb)
            lonlist[lonlist < 0.] += 360.0 

        # Make a height scale
        hgt = np.linspace(self.dict['minAltP'],gph.max().round(),self.dict['nhgt'])
        
        # Interpolate pressure, temperature and Humidity of hgt
        [Pi,Ti,Vi] = processor.intP2H(lvls,hgt,gph,tmp,vpr,self.dict,verbose=verb)
        
        # Calculate the delays
        [DDry,DWet] = processor.PTV2del(Pi,Ti,Vi,hgt,self.dict,verbose=verb)
        if Del in ('comb','Comb'):
            Delfn = DDry+DWet
        elif Del in ('dry','Dry'):
            Delfn = DDry
        elif Del in ('wet','Wet'):
            Delfn = DWet
        else:
            print('Unrecognized delay type')
            sys.exit(1)

        if self.grib in ('NARR'):
            assert False, 'Need to check narr.intdel'
            [Delfn,latlist,lonlist] = narr.intdel(hgt,latlist,lonlist,Delfn)
        
        # Save things
        self.Delfn = Delfn
        self.latlist = latlist
        self.lonlist = lonlist
        self.hgt = hgt
        self.Pi = Pi
        self.Ti = Ti
        self.Vi = Vi
        self.verb = verb

        # All done
        return

    def getdelay(self, dataobj, wvl=np.pi*4., writeStations=True):
        '''Write delay to a matrix / HDF5 object or a file directly. 
           Bilinear Interpolation is used.
                
                Args:   
                        * dataobj  : Final output. (str or HDF5 or np.array)
                                     If str, output is written to file.
                Kwargs:         
                        * wvl      : Wavelength in meters. 
                                     Default output results in delay in meters.
                
                .. note::
                        If dataobj is string, output is written to the file.
                        If np.array or HDF5 object, it should be of size (ny,nx).
        '''
        
        # To know the type of incidence (float, array)
        if isinstance(self.inc,float) \
                or isinstance(self.inc,np.float64) \
                or isinstance(self.inc,np.float32):
            if self.verb:
                print('Incidence is a number: {}'.format(self.inc))
            cinc = np.cos(self.inc*np.pi/180.)
            incFileFlag = 'number'
        else:
            if self.verb:
                print('Assuming incidence can be directly accessed like an array')
            incFileFlag = 'array'
            assert self.inc.shape==(self.ny, self.nx), \
                    'Input table for incidence is not of the right shape'
        
        # Get some info from the dictionary
        minAltp = self.dict['minAltP']

        # Check output and open file if necessary
        outFile = isinstance(dataobj,str)
        if outFile:
            fout = open(dataobj,'wb')
            dout = np.zeros((self.ny,self.nx))
        else:
            assert ((dataobj.shape[0]==self.ny) & (dataobj.shape[1]==self.nx)), \
                    'PyAPS: Not a valid data object.'
            dout = dataobj

        #######################################################################################
        # BILINEAR INTERPOLATION

        # Create the 1d interpolator to interpolate delays in altitude direction
        if self.verb:
            print('PROGRESS: FINE INTERPOLATION OF HEIGHT LEVELS')
        intp_1d = si.interp1d(self.hgt,self.Delfn,kind='cubic',axis=2)

        # Interpolate the delay function every meter, for each station
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
        linearint = si.RegularGridInterpolator((Latu[::-1], Lonu,kh), 
                                               self.Delfn_1m[::-1,:,:], 
                                               method='linear', 
                                               bounds_error=False, 
                                               fill_value = 0.0)

        # Show progress bar
        if self.verb:
            toto = utils.ProgressBar(maxValue=self.ny)
            print('PROGRESS: MAPPING THE DELAY')

        # Loop on the lines
        for m in range(self.ny):

            # Update progress bar
            if self.verb:
                toto.update(m,every=5)
            
            # Get latitude and longitude arrays
            lati = self.lat[m,:]*self.mask[m,:]
            loni = self.lon[m,:]*self.mask[m,:]

            # Remove negative values
            loni[loni<0.] += 360.
            
            # Remove NaN values
            ii = np.where(np.isnan(lati))
            jj = np.where(np.isnan(loni))
            xx = np.union1d(ii,jj)
            lati[xx]=0.0
            loni[xx]=0.0

            # Get incidence if file provided
            if incFileFlag=='array':
                cinc = np.cos(self.mask[m,:]*self.inc[m,:]*np.pi/180.)

            # Make the bilinear interpolation
            D = self.dem[m,:]
            val = linearint(np.vstack((lati, loni, D)).T)*np.pi*4.0/(cinc*wvl)
            val[xx] = np.nan

            # Write outfile
            if outFile:
               resy = val.astype(np.float32)
               resy.tofile(fout)
            else:
                dataobj[m,:] = val

        if self.verb:
            toto.close()

        if outFile:
            fout.close()

        # All done
        return

        #######################################################################################
        # WRITE STATIONS INFOS TO FILE 

        #if writeStations:
        #    stations_latfile = os.path.join(os.path.dirname(self.gfile),'latStations.txt')
        #    stations_lonfile = os.path.join(os.path.dirname(self.gfile),'lonStations.txt')
        #    if self.verb:
        #        print('SAVING STATIONS LATITUDES in: {}'.format(stations_latfile))
        #        print('SAVING STATIONS LONGITUDES in: {}'.format(stations_lonfile))
        #    np.savetxt(stations_latfile, self.latlist, fmt=['%1.2f'])
        #    np.savetxt(stations_lonfile, self.lonlist, fmt=['%1.2f'])

        #    for station in range(0,len(self.lonlist)):
        #        # Open file output
        #        sfile = 'station{}_{}.txt'.format(station, 
        #                    os.path.splitext(os.path.basename(self.gfile))[0])
        #        stationFile = os.path.join(os.path.dirname(self.gfile), sfile)
        #        myfile = open(stationFile, 'w')
        #        # Iterate over altitude level
        #        for i in range(self.Delfn_1m.shape[0]):
        #            altitudeValue = kh[i]
        #            phaseValue = self.Delfn_1m[i,:,:].flatten()[station]
        #            # Write file
        #            myfile.write("{} {}\n".format(altitudeValue, phaseValue))
        #        myfile.close()
        #
        #    # Save kh into file
        #    khfile = os.path.join(os.path.dirname(self.gfile), 'kh.txt')
        #    np.savetxt(khfile, kh, newline='\n', fmt="%s")

        ########################################################################################
        ## CUBIC INTERPOLATION 

        ## Create bicubic interpolator
        #if self.verb:
        #    print('PROGRESS: CREATE THE CUBIC INTERPOLATION FUNCTION')

        ## Resize
        #lonn, hgtn = np.meshgrid(self.lonlist, self.hgt)
        #latn, hgtn = np.meshgrid(self.latlist, self.hgt)

        ## Define a cubic interpolating function on the 3D grid: ((x, y, z), data)
        #cubicint = si.Rbf(lonn.flatten(), latn.flatten(), hgtn.flatten(), self.Delfn.flatten(),kind='cubic',fill_value = 0.0)

        ## Show progress bar
        #if self.verb:
        #    toto = utils.ProgressBar(maxValue=self.ny)
        #    print('PROGRESS: MAPPING THE DELAY')

        ## Loop on the lines
        #for m in range(self.ny):

        #    # Update progress bar
        #    if self.verb:
        #        toto.update(m,every=5)

        #    ###############################################
        #    ## Get values of the m line

        #    ## Longitude
        #    #print(self.lonlist.shape)
        #    #Lonu = np.unique(self.lonlist)
        #    #nLon = len(Lonu)
        #    #lonlisti = Lonu

        #    ## Latitude by iterating over self.latlist
        #    #if m == 0:
        #    #    pos1 = 0
        #    #    pos2 = nLon
        #    #else:
        #    #    pos1 = pos1 + nLon
        #    #    pos2 = pos2 + nLon
        #    #lonlisti = Lonu
        #    #latlisti = self.latlist[pos1:pos2]

        #    ## Height
        #    #hgtlisti = self.hgt
        #    #
        #    ## Delay by iterating over self.Delfn
        #    #print(self.Delfn.shape)
        #    #Delfni = (self.Delfn[pos1:pos2,:]).T

        #    ###############################################
        #    ## Define a cubic interpolating function on the 3D grid: ((x, y, z), data)
        #    #
        #    #lonn, hgtn = np.meshgrid(lonlisti, hgtlisti)
        #    #latn, hgtn = np.meshgrid(latlisti, hgtlisti)
        #    #cubicint = si.Rbf(lonn.flatten(), latn.flatten(), hgtn.flatten(), Delfni.flatten(), kind='cubic',fill_value = 0.0)
        #    #
        #    ###############################################

        #    # Get latitude and longitude arrays
        #    lati = self.lat[m,:]*self.mask[m,:]
        #    loni = self.lon[m,:]*self.mask[m,:]

        #    # Remove negative values
        #    loni[loni<0.] += 360.
        #    
        #    # Remove NaN values
        #    ii = np.where(np.isnan(lati))
        #    jj = np.where(np.isnan(loni))
        #    xx = np.union1d(ii,jj)
        #    lati[xx]=0.0
        #    loni[xx]=0.0

        #    # Get incidence if file provided
        #    if incFileFlag=='array':
        #        cinc = np.cos(self.mask[m,:]*self.inc[m,:]*np.pi/180.)

        #    # Make the interpolation
        #    hgti = self.dem[m,:]
        #    val = cubicint(loni.flatten(), lati.flatten(), hgti.flatten())
        #    val[xx] = np.nan

        #    # Write outfile
        #    if outFile:
        #       resy = val.astype(np.float32)
        #       resy.tofile(fout)
        #    else:
        #        dataobj[m,:] = val

        #if self.verb:
        #    toto.close()

        #if outFile:
        #    fout.close()

        ## All done
        #return
