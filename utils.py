############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
###############Reading input RSC file for radar###############
def rd_rsc(inname,full=False,verbose=False):
    '''Reading a ROI-PAC style RSC file.

    Args:
        * inname (str): Path to the RSC file.

    Returns:
        * lat (np.array) : Array of lat of the 4 corners.
        * lon (np.array) : Array of lon of the 4 corners.
        * nx  (np.int)   : Number of range bins.
        * ny  (np.int)   : Number of azimuth lines.
        * dpix (np.float): Average pixel spacing.

    .. note:: 
        Currently set up to work with SIM_nrlks.hgt from ROI-PAC.'''

    if verbose:
        print("PROGRESS: READING %s RSC FILE" %inname)

    rpacdict = {}
    infile = open(inname+'.rsc','r')
    line = infile.readline()
    rpacdict['LAT_REF1']=0.0
    rpacdict['LON_REF1']=0.0
    rpacdict['LAT_REF2']=0.0
    rpacdict['LON_REF2']=0.0
    rpacdict['LAT_REF3']=0.0
    rpacdict['LON_REF3']=0.0
    rpacdict['LAT_REF4']=0.0
    rpacdict['LON_REF4']=0.0
    while line:
       llist = line.split()
       if len(llist)>0 :
          rpacdict[llist[0]] = llist[1]
       line = infile.readline()
    infile.close()

    nx = np.int(rpacdict['WIDTH'])
    ny = np.int(rpacdict['FILE_LENGTH'])
    lat=np.zeros((4,1))
    lon=np.zeros((4,1))
    lat[0] = np.float(rpacdict['LAT_REF1'])
    lon[0] = np.float(rpacdict['LON_REF1'])
    lat[1] = np.float(rpacdict['LAT_REF2'])
    lon[1] = np.float(rpacdict['LON_REF2'])
    lat[2] = np.float(rpacdict['LAT_REF3'])
    lon[2] = np.float(rpacdict['LON_REF3'])
    lat[3] = np.float(rpacdict['LAT_REF4'])
    lon[3] = np.float(rpacdict['LON_REF4'])
#    rgpix = np.float(rpacdict['RANGE_PIXEL_SIZE'])
#    azpix = np.float(rpacdict['AZIMUTH_PIXEL_SIZE'])
    rgpix = 90.0
    azpix = 90.0
    dpix = np.sqrt(6.25*rgpix*rgpix+azpix*azpix)
    if full:
        return lon,lat,nx,ny,dpix,rpacdict
    else:
        return lon,lat,nx,ny,dpix


#######################Finished rd_rsc###########################

###############Reading input RSC file for geo###############
def geo_rsc(inname,full=False,verbose=False):
    '''Reading a ROI-PAC style geocoded rsc file.

    Args:
        * inname (str): Path to the RSC file.

    Returns:
        * lon (np.array) : Array of min and max lon values.
        * lat (np.array) : Array of min and max lat values.
        * nx  (np.int)   : Number of lon bins.
        * ny  (np.int)   : Number of lat bins.

        .. note:: 
            Currently set up to work with dem.rsc file from ROI-PAC.'''

    if verbose:
        print("PROGRESS: READING %s RSC FILE" %inname)

    rpacdict = {}
    infile = open(inname+'.rsc','r')
    line = infile.readline()
    while line:
        llist = line.split()
        if len(llist)>0 :
            rpacdict[llist[0]] = llist[1]
        line = infile.readline()
    infile.close()

    nx = np.int(rpacdict['WIDTH'])
    ny = np.int(rpacdict['FILE_LENGTH'])
    lat=np.zeros((2,1))
    lon=np.zeros((2,1))
    lat[1] = np.float(rpacdict['Y_FIRST'])
    lon[0] = np.float(rpacdict['X_FIRST'])
    if(lon[0] < 0):
        lon[0] = lon[0] + 360.0

    dx = np.float(rpacdict['X_STEP'])
    dy = np.float(rpacdict['Y_STEP'])
    
    lat[0] = lat[1] + dy*ny
    lon[1] = lon[0] + dx*nx
    
    if full:
        return lon,lat,nx,ny,rpacdict
    else:
        return lon,lat,nx,ny

#######################Finished geo_rsc###########################


##########Conversion from Geo coordinate to Radar coordinates######
def lonlat2rdr(lon, lat, lonlist, latlist, plotflag=False):
    '''
    Transfer the lat lon coordinates of the weather stations into the index
    range and azimuth coordinates of the radar scene.

    Args:
        * lon               : Longitude array of the radar scene size=(ny,nx)
        * lat               : Latitude array of the radar scene size=(ny,nx)
        * lonlist           : Longitude list of weeather stations
        * latlist           : Latitude list of the weather stations.

    Kwargs:
        * plotflag          : Plot something to check. (default is False)
    
    note:
        Mapping function is :
        * range = a1*lat+b1*lon+c1
        * azimu = a2*lat+b2*lon+c2
        * First  point is (1,1)   i.e. Near Range, First Lane <==> Lat[0,0],Lon[0,0]
        * Second point is (nx,1)  i.e. Far Range, First Lane <==> Lat[0,-1],Lon[0,-1]
        * Third  point is (1,ny)  i.e. Near Range, Last Lane <==> Lat[-1,0],Lon[-1,0]
        * Fourth point is (nx,ny) i.e. Far Range, Last Lane <==> Lat[-1,-1],Lon[-1,-1] 
    '''

    # Size
    ny, nx = lon.shape

    # Mapping function
    A=np.array([ [lat[0,0], lon[0,0], 1.],
                 [lat[0,-1], lon[0,-1], 1.],
                 [lat[-1,0], lon[-1,0], 1.],
                 [lat[-1,-1], lon[-1,-1], 1.]])
    b=np.array([[1., 1.],[nx, 1.],[1., ny],[nx,ny]])
    mfcn=np.linalg.lstsq(A,b)[0]

    #Get grid points xi yi coordinates from this mapping function
    A=np.vstack([latlist, lonlist, np.ones(len(lonlist),)]).T
    xi=np.dot(A,mfcn[:,0])
    yi=np.dot(A,mfcn[:,1])
    
    # Plot y/n
    if plotflag:
        import matplotlib.patches as ptch

        plt.figure()
        plt.subplot(211)
        plt.scatter(lonlist,latlist,s=8,c=np.cumsum(np.ones(len(lonlist),)))
        xline=[lon[0,0],lon[0,-1],lon[-1,-1],lon[-1,0],lon[0,0]]
        yline=[lat[0,0],lat[0,-1],lat[-1,-1],lat[-1,0],lat[0,0]]
        plt.plot(xline,yline,'-r')
        plt.title('Area of interest and {} stations used'.format(len(lonlist)))
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.subplot(212)
        plt.scatter(xi,yi,s=8,c=np.cumsum(np.ones(len(lonlist),)))
        p = ptch.Rectangle((1,1),lon.shape[1],lon.shape[0],
                                 edgecolor="Red",fill=False)
        plt.gca().add_patch(p)
        plt.title('Area of interest in Radar Geometry')
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.show()

    # All done
    return np.array(xi).astype(float), np.array(yi).astype(float)

def glob2rdr(nx,ny,lat,lon,latl,lonl,plotflag='n'):
    '''Transfert these latlon coordinates into radar geometry (xi,yi) with a 
       simple linear transformation given the first pixel and the pixel 
       spacing of the simulation.

    Args:
        * nx   (np.int)   : Number of range bins.
        * ny   (np.int)   : Number of azimuth lines.
        * lat  (np.array) : Array of latitudes of the corners
        * lon  (np.array) : Array of longitudes of the corners
        * latl (np.array) : Latitudes of the stations.
        * lonl (np.array) : Longitudes of the stations.

    Kwargs:
        * plotflag (bool) : Plot the stations distribution.

    Returns:
        * xi (np.array)   : Position of stations in range.
        * yi (np.array)   : Position of stations in azimuth.

    .. note::
        Mapping function is :
        * range = a1*lat+b1*lon+c1
        * azimu = a2*lat+b2*lon+c2
        * First  point is (1,1)   i.e. Near Range, First Lane <==> Lat[1],Lon[1]
        * Second point is (nx,1)  i.e. Far Range, First Lane <==> Lat[2],Lon[2]
        * Third  point is (1,ny)  i.e. Near Range, Last Lane <==> Lat[3],Lon[3]
        * Fourth point is (nx,ny) i.e. Far Range, Last Lane <==> Lat[4],Lon[4] '''

    # Mapping function
    A=np.hstack([lat,lon,np.ones((4,1))])
    b=np.array([[1, 1],[nx, 1],[1, ny],[nx,ny]])
    mfcn=np.linalg.lstsq(A,b)[0]

    #Get grid points xi yi coordinates from this mapping function
    nstn=len(latl)
    A=np.array([latl, lonl, np.ones(nstn)]).T
    xi=np.dot(A,mfcn[:,0])
    yi=np.dot(A,mfcn[:,1])

    if plotflag in ('y','Y'):
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(lonl,latl,s=8,c='k');
        xline=[lon[0],lon[1],lon[3],lon[2],lon[0]]
        yline=[lat[0],lat[1],lat[3],lat[2],lat[0]]
        plt.plot(xline,yline,'-r')
        plt.title('Area of interest and %d stations used'%(nstn))
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.subplot(212)
        plt.scatter(xi,yi,s=8,c='k')
        import matplotlib.patches as ptch
        p = ptch.Rectangle((1,1),nx,ny,edgecolor="Red",fill=False)
        plt.gca().add_patch(p)
        plt.title('Area of interest in Radar Geometry')
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.show()

    return xi,yi
#################Completed transforming geo 2 radar##################

##########Conversion from Geo coordinate to Radar coordinates######
def rdr2glob(wid,lgt,lat,lon,x,y,plotflag='n'):   #nx,ny,lat,lon,latl,lonl,plotflag='n'):
    '''Transfert these radar geometry (x,y) coordinates in lat/lon coordinates with a simple linear transformation given the image width/length and the lat/lon coordinates of the 4 corners

    Args:
            * wid  (np.int)   : Width of the image (i.e. number of range bins)
            * lgt  (np.int)   : Length of the image (i.e. number of azimuth lines)
            * lat  (np.array) : Array of latitudes of the corners
            * lon  (np.array) : Array of longitudes of the corners
            * rang (np.array) : Range of the points to transfert
            * azim (np.array) : Azimuth of the points to transfert

    Kwargs:
    * plotflag (bool) : Plot the stations distribution.

    Returns:
            * loni (np.array)   : Longitude of the points.
            * lati (np.array)   : Latitude of the points.

    .. note::
            Mapping function is :
    * lat = a1*rang + b1*azim + c1
    * lon = a2*rang + b2*azim + c2
            * First  point is (1,1)   i.e. Near Range, First Lane <==> Lat[1],Lon[1]
            * Second point is (nx,1)  i.e. Far Range, First Lane <==> Lat[2],Lon[2]
            * Third  point is (1,ny)  i.e. Near Range, Last Lane <==> Lat[3],Lon[3]
            * Fourth point is (nx,ny) i.e. Far Range, Last Lane <==> Lat[4],Lon[4] '''

    A=np.array([[1, 1, 1.],[wid, 1, 1.],[1, lgt, 1.],[wid, lgt, 1.]])
    #b=np.array([[lat[0],lon[0]],[lat[1],lon[1]],[lat[2],lon[2]],[lat[3],lon[3]]])
    b = np.hstack((lat,lon))
    mfcn=np.linalg.lstsq(A,b)[0]

    #Get grid points xi yi coordinates from this mapping function
    nstn=len(x)
    A=np.array([x, y, np.ones((nstn,1))]).T
    lati=np.dot(A,mfcn[:,0])
    loni=np.dot(A,mfcn[:,1])

    if plotflag in ('y','Y'):
        plt.figure(1)
        plt.subplot(211)
        plt.scatter(lon,lat,s=8,c='k');
        xline=[lon[0],lon[1],lon[3],lon[2],lon[0]]
        yline=[lat[0],lat[1],lat[3],lat[2],lat[0]]
        plt.scatter(loni,lati,s=8,c='r')
        plt.plot(xline,yline,'-r')
        plt.title('Area of interest in geographic coordinates')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        plt.subplot(212)
        plt.scatter(x,y,s=8,c='k')
        import matplotlib.patches as ptch
        p = ptch.Rectangle((1,1),wid,lgt,edgecolor="Red",fill=False)
        plt.gca().add_patch(p)
        plt.title('Area of interest in Radar Geometry')
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.show()

    return loni,lati
#################Completed transforming geo 2 radar##################


###########################Simple progress bar######################
class ProgressBar:
    """ Creates a text-based progress bar. Call the object with 
    the simple `print'command to see the progress bar, which looks 
    something like this:
    [=======> 22% ]
    You may specify the progress bar's width, min and max values on init.
    
    .. note::
        Code originally from http://code.activestate.com/recipes/168639/"""

    def __init__(self, minValue = 0, maxValue = 100, totalWidth=80):
        self.progBar = "[]" # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.reset()

    def reset(self):
        self.start_time = time.time()
        self.amount = 0 # When amount == max, we are 100% done
        self.updateAmount(0) # Build progress bar string

    def updateAmount(self, newAmount = 0):
        """ Update the progress bar with the new amount (with min and max
        values set at initialization; if it is over or under, it takes the
        min or max value as a default. """
        if newAmount < self.min:
            newAmount = self.min
        if newAmount > self.max:
            newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = np.float(self.amount - self.min)
        percentDone = (diffFromMin / np.float(self.span)) * 100.0
        percentDone = np.int(np.round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2 - 18
        numHashes = (percentDone / 100.0) * allFull
        numHashes = np.int(np.round(numHashes))

        # Build a progress bar with an arrow of equal signs; special cases for
        # empty and full
        if numHashes == 0:
            self.progBar = '[>%s]' % (' '*(allFull-1))
        elif numHashes == allFull:
            self.progBar = '[%s]' % ('='*allFull)
        else:
            self.progBar = '[%s>%s]' % ('='*(numHashes-1),
                            ' '*(allFull-numHashes))
            # figure out where to put the percentage, roughly centered
            percentPlace = (len(self.progBar) // 2) - len(str(percentDone))
            percentString = ' ' + str(percentDone) + '% '
            elapsed_time = time.time() - self.start_time
            # slice the percentage into the bar
            self.progBar = ''.join([self.progBar[0:percentPlace],percentString, self.progBar[percentPlace+len(percentString):], ])
            if percentDone > 0:
                self.progBar += ' %6ds / %6ds' % (int(elapsed_time),int(elapsed_time*(100.//percentDone-1)))

    def update(self, value, every=1):
        """ Updates the amount, and writes to stdout. Prints a
         carriage return first, so it will overwrite the current
          line in stdout."""
        if value % every == 0 or value >= self.max:
            self.updateAmount(newAmount=value)
            sys.stdout.write('\r' + self.progBar)
            sys.stdout.flush()

    def close(self):
        """Prints a blank space at the end to ensure proper printing
        of future statements."""
        print(' ')

################################End of progress bar class####################################
