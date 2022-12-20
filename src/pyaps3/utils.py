############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
# Add more utils functions, Zhang Yunjun, April 2019
#     get_lat_lon()
#     read_(meta)data()
#     read_isce_xml/data/lalo_ref()
#     read_roipac_data/rsc()


import os
import sys
import time
import numpy as np



def snwe2str(snwe):
    """Get area extent in string"""
    if not snwe:
        return None
    s, n, w, e = snwe

    area = ''
    area += '_S{}'.format(abs(s)) if s < 0 else '_N{}'.format(abs(s))
    area += '_S{}'.format(abs(n)) if n < 0 else '_N{}'.format(abs(n))
    area += '_W{}'.format(abs(w)) if w < 0 else '_E{}'.format(abs(w))
    area += '_W{}'.format(abs(e)) if e < 0 else '_E{}'.format(abs(e))

    return area



############### Read ISCE / ROIPAC file ###############

def get_lat_lon(metafile):
    """Get 2D lat and lon from rsc/xml file
    Args:
        * metafile (str or dict) : path to metadata file, or dict of metadata
    """
    if isinstance(metafile, str):
        fext = os.path.splitext(metafile)[1]
        if fext == '.rsc':
            meta = read_roipac_rsc(metafile)
        elif fext == '.xml':
            meta = read_isce_xml(metafile)
        else:
            raise ValueError('unrecognized file extension: {}'.format(fext))
    elif isinstance(metafile, dict):
        meta = {}
        for key, value in metafile.items():
            meta[key] = value
    length, width = int(meta['FILE_LENGTH']), int(meta['WIDTH'])

    # get bbox
    if 'Y_FIRST' in meta.keys():
        # geo coordinates
        lat0 = float(meta['Y_FIRST'])
        lon0 = float(meta['X_FIRST'])
        lat_step = float(meta['Y_STEP'])
        lon_step = float(meta['X_STEP'])
        lat1 = lat0 + lat_step * (length - 1)
        lon1 = lon0 + lon_step * (width - 1)
    else:
        # radar coordinates
        lats = [float(meta['LAT_REF{}'.format(i)]) for i in [1,2,3,4]]
        lons = [float(meta['LON_REF{}'.format(i)]) for i in [1,2,3,4]]
        lat0 = np.mean(lats[0:2])
        lat1 = np.mean(lats[2:4])
        lon0 = np.mean(lons[0:3:2])
        lon1 = np.mean(lons[1:4:2])

    # bbox --> 2D array
    lat, lon = np.mgrid[lat0:lat1:length*1j,
                        lon0:lon1:width*1j]
    return lat, lon


def read_data(fname, dname='inc'):
    """Read 2D data from ISCE / ROIPAC file"""
    # get meta file extension
    meta_exts = [i for i in ['.xml', '.rsc'] if os.path.isfile(fname+i)]
    if len(meta_exts) == 0:
        raise FileNotFoundError('No metadata file found for data file: {}'.format(fname))
    meta_ext = meta_exts[0]

    if meta_ext == '.xml':
        data = read_isce_data(fname, dname=dname)
    elif meta_ext == '.rsc':
        data = read_roipac_data(fname)
    return data


def read_metadata(fname, latfile=None, lonfile=None, full=False, verbose=False):
    '''Reading metadata for ISCE or ROIPAC style data file

    Args:
        * fname (str): Path to the ROIPAC or ISCE data file.

    Returns:
        * lat  (np.ndarray): Array of lat of the 4 corners.
        * lon  (np.ndarray): Array of lon of the 4 corners.
        * nx   (int)       : Number of range bins.
        * ny   (int)       : Number of azimuth lines.
        * dpix (float)     : Average pixel spacing.
        * meta (dict)      : Dictionary of metadata [return if full==True]
    '''
    # get meta file extension
    meta_exts = [i for i in ['.rsc','.xml'] if os.path.isfile(fname+i)]
    if len(meta_exts) == 0:
        raise FileNotFoundError('No metadata file found for data file: {}'.format(fname))
    meta_ext = meta_exts[0]

    if verbose:
        print("PROGRESS: READING {} FILE".format(fname+meta_ext))

    # read metadata
    meta = {}
    if meta_ext == '.rsc':
        meta = read_roipac_rsc(fname+meta_ext)

    elif meta_ext == '.xml':
        meta = read_isce_xml(fname+meta_ext)

        # default latfile
        if not latfile:
            latfile = os.path.join(os.path.dirname(fname),'lat.rdr')
            if not os.path.isfile(latfile):
                raise FileNotFoundError("No latitude file found in ISCE style!")

        # default lonfile
        if not lonfile:
            lonfile = os.path.join(os.path.dirname(fname),'lon.rdr')
            if not os.path.isfile(lonfile):
                raise FileNotFoundError("No longitude file found in ISCE style!")

        meta.update(get_isce_lalo_ref(latfile, lonfile))

    # prepare output
    nx = int(meta['WIDTH'])
    ny = int(meta['FILE_LENGTH'])
    lat = np.zeros((4,1))
    lon = np.zeros((4,1))
    lat[0] = float(meta['LAT_REF1'])
    lon[0] = float(meta['LON_REF1'])
    lat[1] = float(meta['LAT_REF2'])
    lon[1] = float(meta['LON_REF2'])
    lat[2] = float(meta['LAT_REF3'])
    lon[2] = float(meta['LON_REF3'])
    lat[3] = float(meta['LAT_REF4'])
    lon[3] = float(meta['LON_REF4'])
    rgpix = 90.0
    azpix = 90.0
    dpix = np.sqrt(6.25*rgpix*rgpix + azpix*azpix)
    if full:
        return lon,lat,nx,ny,dpix,meta
    else:
        return lon,lat,nx,ny,dpix


def read_isce_xml(xmlfile):
    """Read ISCE XML file into dict.
    Add from mintpy.utils.readfile.py by Zhang Yunjun
    """
    import xml.etree.ElementTree as ET
    meta = {}
    root = ET.parse(xmlfile).getroot()
    if root.tag.startswith('image'):
        for child in root.findall('property'):
            key = child.get('name')
            value = child.find('value').text
            meta[key] = value

        # Read lat/lon info for geocoded file
        # in form: root/component coordinate*/property name/value
        for coord_name, prefix in zip(['coordinate1', 'coordinate2'], ['X', 'Y']):
            child = root.find("./component[@name='{}']".format(coord_name))
            if ET.iselement(child):
                v_step  = float(child.find("./property[@name='delta']").find('value').text)
                v_first = float(child.find("./property[@name='startingvalue']").find('value').text)
                if abs(v_step) < 1. and abs(v_step) > 1e-7:
                    meta['{}_STEP'.format(prefix)] = v_step
                    meta['{}_FIRST'.format(prefix)] = v_first - v_step / 2.

    # convert key name from isce to roipac
    isce2roipacKeyDict = {
        'width':'WIDTH',
        'length':'FILE_LENGTH',
    }
    for key,value in isce2roipacKeyDict.items():
        meta[value] = meta[key]
    return meta


def read_isce_data(fname, dname=None):
    """Read ISCE data file"""
    # read xml file
    xml_file = fname+'.xml'
    meta = read_isce_xml(xml_file)

    # get data_type
    dataTypeDict = {
        'byte': 'bool_',
        'float': 'float32',
        'double': 'float64',
        'cfloat': 'complex64',
    }
    data_type = dataTypeDict[meta['data_type'].lower()]

    width = int(meta['width'])
    length = int(meta['length'])
    num_band = int(meta['number_bands'])
    band = 1
    if fname.startswith('los') and dname and dname.startswith('az'):
        band = 2

    # read
    data = np.fromfile(fname, dtype=data_type, count=length*width*num_band).reshape(-1, width*num_band)
    data = data[:, width*(band-1):width*band]
    return data


def get_isce_lalo_ref(lat_file, lon_file):
    """Get LAT/LON_REF1/2/3/4 value from ISCE lat/lon.rdr file
    Add from mintpy.prep_isce.py by Zhang Yunjun
    """

    def get_nonzero_row_number(data, buffer=2):
        """Find the first and last row number of rows without zero value
        for multiple swaths data
        """
        if np.all(data):
            r0, r1 = 0 + buffer, -1 - buffer
        else:
            row_flag = np.sum(data != 0., axis=1) == data.shape[1]
            row_idx = np.where(row_flag)[0]
            r0, r1 = row_idx[0] + buffer, row_idx[-1] - buffer
        return r0, r1

    meta = {}
    # read LAT/LON_REF1/2/3/4 from lat/lonfile
    lat = read_isce_data(lat_file)
    r0, r1 = get_nonzero_row_number(lat)
    meta['LAT_REF1'] = str(lat[r0, 0])
    meta['LAT_REF2'] = str(lat[r0, -1])
    meta['LAT_REF3'] = str(lat[r1, 0])
    meta['LAT_REF4'] = str(lat[r1, -1])

    lon = read_isce_data(lon_file)
    r0, r1 = get_nonzero_row_number(lon)
    meta['LON_REF1'] = str(lon[r0, 0])
    meta['LON_REF2'] = str(lon[r0, -1])
    meta['LON_REF3'] = str(lon[r1, 0])
    meta['LON_REF4'] = str(lon[r1, -1])
    return meta


def read_roipac_rsc(rscfile):
    """Read ROIPAC style RSC file into dict"""
    meta = {}
    f = open(rscfile,'r')
    line = f.readline()
    while line:
        llist = line.split()
        if len(llist)>0 :
            meta[llist[0]] = llist[1]
        line = f.readline()
    f.close()
    return meta


def read_roipac_data(fname):
    """Read ROIPAC data file"""
    # get length / width
    rsc_file = fname+'.rsc'
    meta = read_roipac_rsc(rsc_file)
    length = int(meta['FILE_LENGTH'])
    width = int(meta['WIDTH'])

    # get datatype / band
    data_type = 'float32'
    num_band = 1
    band = 1
    fext = os.path.splitext(fname)[1]
    if fext == '.dem':
        data_type = 'int16'
    elif fext == '.hgt':
        num_band = 2
        band = 2

    # read
    data = np.fromfile(fname, dtype=data_type, count=length*width*num_band).reshape(-1, width*num_band)
    data = data[:, width*(band-1):width*band]
    return data

############### Read ISCE / ROIPAC file - end #########



############### obsolete: RSC for radar ###############

def rd_rsc(inname,full=False,verbose=False):
    '''Reading a ROI-PAC style RSC file.

    Args:
        * inname (str): Path to the RSC file.

    Returns:
        * lat  (np.ndarray): Array of lat of the 4 corners.
        * lon  (np.ndarray): Array of lon of the 4 corners.
        * nx   (int)       : Number of range bins.
        * ny   (int)       : Number of azimuth lines.
        * dpix (float)     : Average pixel spacing.

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

    # prepare output
    nx = int(rpacdict['WIDTH'])
    ny = int(rpacdict['FILE_LENGTH'])
    lat=np.zeros((4,1))
    lon=np.zeros((4,1))
    lat[0] = float(rpacdict['LAT_REF1'])
    lon[0] = float(rpacdict['LON_REF1'])
    lat[1] = float(rpacdict['LAT_REF2'])
    lon[1] = float(rpacdict['LON_REF2'])
    lat[2] = float(rpacdict['LAT_REF3'])
    lon[2] = float(rpacdict['LON_REF3'])
    lat[3] = float(rpacdict['LAT_REF4'])
    lon[3] = float(rpacdict['LON_REF4'])
    rgpix = 90.0   # float(rpacdict['RANGE_PIXEL_SIZE'])
    azpix = 90.0   # float(rpacdict['AZIMUTH_PIXEL_SIZE'])
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
        * lon (np.ndarray): Array of min and max lon values.
        * lat (np.ndarray): Array of min and max lat values.
        * nx  (int)       : Number of lon bins.
        * ny  (int)       : Number of lat bins.

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

    nx = int(rpacdict['WIDTH'])
    ny = int(rpacdict['FILE_LENGTH'])
    lat=np.zeros((2,1))
    lon=np.zeros((2,1))
    lat[1] = float(rpacdict['Y_FIRST'])
    lon[0] = float(rpacdict['X_FIRST'])
    if(lon[0] < 0):
        lon[0] = lon[0] + 360.0

    dx = float(rpacdict['X_STEP'])
    dy = float(rpacdict['Y_STEP'])

    lat[0] = lat[1] + dy*ny
    lon[1] = lon[0] + dx*nx
    
    if full:
        return lon,lat,nx,ny,rpacdict
    else:
        return lon,lat,nx,ny

############### obsolete: RSC for radar - end #########



############### obsolete: simple geo2radar ############

def lonlat2rdr(lon, lat, lonlist, latlist, plotflag=False):
    '''
    Transfer the lat lon coordinates of the weather stations into the index
    range and azimuth coordinates of the radar scene.

    Args:
        * lon               : Longitude array of the radar scene size=(ny,nx)
        * lat               : Latitude array of the radar scene size=(ny,nx)
        * lonlist           : Longitude list of the weeather stations
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
        from matplotlib import pyplot as plt, patches

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
        p = patches.Rectangle((1,1),lon.shape[1],lon.shape[0],
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
        * nx   (int)        : Number of range bins.
        * ny   (int)        : Number of azimuth lines.
        * lat  (np.ndarray) : Array of latitudes of the corners
        * lon  (np.ndarray) : Array of longitudes of the corners
        * latl (np.ndarray) : Latitudes of the stations.
        * lonl (np.ndarray) : Longitudes of the stations.

    Kwargs:
        * plotflag (bool)   : Plot the stations distribution.

    Returns:
        * xi (np.ndarray)   : Position of stations in range.
        * yi (np.ndarray)   : Position of stations in azimuth.

    .. note::
        Mapping function is :
        * range = a1*lat+b1*lon+c1
        * azimu = a2*lat+b2*lon+c2
        * First  point is (1,1)   i.e. Near Range, First Lane <==> Lat[1],Lon[1]
        * Second point is (nx,1)  i.e. Far Range, First Lane <==> Lat[2],Lon[2]
        * Third  point is (1,ny)  i.e. Near Range, Last Lane <==> Lat[3],Lon[3]
        * Fourth point is (nx,ny) i.e. Far Range, Last Lane <==> Lat[4],Lon[4] '''

    # Mapping function
    A = np.hstack([lat,lon,np.ones((4,1))])
    b = np.array([[1, 1],[nx, 1],[1, ny],[nx,ny]])
    mfcn = np.linalg.lstsq(A,b,rcond=None)[0]

    #Get grid points xi yi coordinates from this mapping function
    nstn = latl.size
    A = np.array([latl.reshape(-1,1), lonl.reshape(-1,1), np.ones((nstn,1))]).T
    xi = np.dot(A,mfcn[:,0]).reshape(latl.shape)
    yi = np.dot(A,mfcn[:,1]).reshape(latl.shape)

    if plotflag in ('y','Y'):
        from matplotlib import pyplot as plt, patches

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
        p = patches.Rectangle((1,1),nx,ny,edgecolor="Red",fill=False)
        plt.gca().add_patch(p)
        plt.title('Area of interest in Radar Geometry')
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.show()

    return xi,yi

############### obsolete: simple geo2radar - end ######



############### obsolete: simple radar2geo ############

def rdr2glob(wid,lgt,lat,lon,x,y,plotflag='n'):   #nx,ny,lat,lon,latl,lonl,plotflag='n'):
    '''Transfert these radar geometry (x,y) coordinates in lat/lon coordinates with a 
    simple linear transformation given the image width/length and the lat/lon coordinates
    of the 4 corners

    Args:
            * wid  (int)        : Width of the image (i.e. number of range bins)
            * lgt  (int)        : Length of the image (i.e. number of azimuth lines)
            * lat  (np.ndarray) : Array of latitudes of the corners
            * lon  (np.ndarray) : Array of longitudes of the corners
            * rang (np.ndarray) : Range of the points to transfert
            * azim (np.ndarray) : Azimuth of the points to transfert

    Kwargs:
            * plotflag (bool)   : Plot the stations distribution.

    Returns:
            * loni (np.ndarray) : Longitude of the points.
            * lati (np.ndarray) : Latitude of the points.

    .. note::
            Mapping function is :
            * lat = a1*rang + b1*azim + c1
            * lon = a2*rang + b2*azim + c2
            * First  point is (1,1)   i.e. Near Range, First Lane <==> Lat[1],Lon[1]
            * Second point is (nx,1)  i.e. Far Range, First Lane <==> Lat[2],Lon[2]
            * Third  point is (1,ny)  i.e. Near Range, Last Lane <==> Lat[3],Lon[3]
            * Fourth point is (nx,ny) i.e. Far Range, Last Lane <==> Lat[4],Lon[4] '''

    A = np.array([[1, 1, 1.],[wid, 1, 1.],[1, lgt, 1.],[wid, lgt, 1.]])
    b = np.hstack((lat,lon))
    mfcn = np.linalg.lstsq(A,b)[0]

    #Get grid points xi yi coordinates from this mapping function
    nstn = len(x)
    A = np.array([x, y, np.ones((nstn,1))]).T
    lati = np.dot(A,mfcn[:,0])
    loni = np.dot(A,mfcn[:,1])

    if plotflag in ('y','Y'):
        from matplotlib import pyplot as plt, patches

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
        p = patches.Rectangle((1,1),wid,lgt,edgecolor="Red",fill=False)
        plt.gca().add_patch(p)
        plt.title('Area of interest in Radar Geometry')
        plt.xlabel('Range')
        plt.ylabel('Azimuth')
        plt.show()

    return loni,lati

############### obsolete: simple radar2geo - end ######



############### Simple progress bar ###################

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
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = int(np.round(percentDone))

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2 - 18
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(np.round(numHashes))

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
            self.progBar = ''.join([self.progBar[0:percentPlace],
                                    percentString,
                                    self.progBar[percentPlace+len(percentString):]])
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

############### Simple progress bar - end #############
