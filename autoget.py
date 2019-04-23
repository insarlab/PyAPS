############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################

import configparser as ConfigParser
import sys
import os.path
import cdsapi
# disable InsecureRequestWarning message from cdsapi
import urllib3
urllib3.disable_warnings()


######Set up variables in model.cfg before using
dpath = os.path.dirname(__file__)
config = ConfigParser.RawConfigParser(delimiters='=')
config.read('%s/model.cfg'%(dpath))


def ECMWFdload(bdate,hr,filedir,model='ERA5',datatype='fc',humidity='Q',snwe=None,flist=None):
    '''
    ECMWF data downloading.

    Args:
        * bdate     : date to download (str)
        * hr        : hour to download (str)
        * filedir   : files directory (str)
        * model     : weather model ('ERA5', 'ERAINT', 'HRES')
        * datatype  : reanalysis data type (an)
        * snwe      : area extent (tuple of int)
        * humidity  : humidity
    '''
    
    #-------------------------------------------
    # Initialize

    # Check data
    assert model in ('ERA5', 'ERAINT', 'HRES'), 'Unknown model for ECMWF: {}'.format(model)

    # Infos for downloading
    if model in 'ERAINT':
        print('WARNING: you are downloading from the old ECMWF platform. '
              'ERA-Interim is deprecated, use ERA-5 instead.')
    if model in 'ERA5':
        print('INFO: You are using the latest ECMWF platform for downloading datasets: '
              'https://cds.climate.copernicus.eu/api/v2')

    #-------------------------------------------
    # Define parameters

    # Humidity
    assert humidity in ('Q','R'), 'Unknown humidity field for ECMWF'
    if humidity in 'Q':
        if model in 'ERA5':
            humidparam = 'specific_humidity'
        else:
            humidparam = 133
    elif humidity in 'R':
        if model in 'ERA5':
            humidparam = 'relative_humidity'
        else:
            humidparam = 157

    # Grid size (only for HRES and ERA-Interim)
    if model in 'HRES':
        gridsize = '0.10/0.10'
    elif model in 'ERAINT':
        gridsize = '0.75/0.75'

    #-------------------------------------------
    # file name
    if not flist:
        flist = []
        for k in range(len(bdate)):
            day = bdate[k]

            if model == 'ERA5':
                fname = os.path.join(filedir, 'ERA-5_{}_{}.grb'.format(day, hr))
            elif model == 'ERAINT':
                fname = os.path.join(filedir, 'ERA-Int_{}_{}.grb'.format(day, hr))
            elif model in 'HRES':
                fname = os.path.join(filedir, 'HRES_{}_{}.grb'.format(day, hr))
            else:
                raise ValueError('unrecognized model input: {}'.format(model))
            flist.append(fname)

    # Iterate over dates
    for k in range(len(bdate)):
        day = bdate[k]
        fname = flist[k]

        #-------------------------------------------
        # CASE 1: request for CDS API client (new ECMWF platform, for ERA-5)    
        if model in 'ERA5':
            url = 'https://cds.climate.copernicus.eu/api/v2'
            key = config.get('CDS', 'key')

            # Contact the server
            c = cdsapi.Client(url=url, key=key)

            # Pressure levels
            pressure_lvls = ['1','2','3','5','7','10','20','30','50', 
                            '70','100','125','150','175','200','225',
                            '250','300','350','400','450','500','550',
                            '600','650','700','750','775','800','825',
                            '850','875','900','925','950','975','1000']

            # Dictionary
            indict = {'product_type'   :'reanalysis',
                      'format'         :'grib',
                      'variable'       :['geopotential','temperature','{}'.format(humidparam)],
                      'pressure_level' : pressure_lvls,
                      'year'           :'{}'.format(day[0:4]),
                      'month'          :'{}'.format(day[4:6]),
                      'day'            :'{}'.format(day[6:8]),
                      'time'           :'{}:00'.format(hr)}

            # download a geographical area subset
            if snwe is not None:
                s, n, w, e = snwe
                indict['area'] = '/'.join(['{:.2f}'.format(x) for x in [n, w, s, e]])

            # Assert grib file not yet downloaded
            if not os.path.exists(fname):
                print('Downloading %d of %d: %s '%(k+1,len(bdate), fname))
                print(indict)

                # Make the request
                c.retrieve('reanalysis-{}-pressure-levels'.format(model.lower()),indict,target=fname)

        #-------------------------------------------
        # CASE 2: request for WEB API client (old ECMWF platform, deprecated, for ERA-Int and HRES)
        else:
            # Contact the server
            from pyaps3.ecmwfapi import ECMWFDataServer
            url = "https://api.ecmwf.int/v1"
            emid = config.get('ECMWF', 'email')
            key = config.get('ECMWF', 'key')
            server = ECMWFDataServer(url=url, key=key, email=emid)

            # Dictionary
            indict = {'dataset'  : '{}'.format(model),
                      'type'     : '{}'.format(datatype),
                      'date'     : '{}'.format(day),
                      'time'     : '{}'.format(hr),
                      'step'     : '0',
                      'levtype'  : "pl",
                      'levelist' : "all",
                      'grid'     : '{}'.format(gridsize),
                      'param'    : '129/130/{}'.format(humidparam),
                      'target'   : '{}'.format(fname)}
            
            # Assert grib file not yet downloaded
            if not os.path.exists(fname):
                print('Downloading %d of %d: %s '%(k+1,len(bdate),fname))
                print(indict)

                # Make the request
                server.retrieve(indict)

    return flist


def MERRAdload(bdate,hr,filedir, hdf=False):
    print('By default, downloads MERRA-2 data.')
    user = config.get('MERRA', 'user')
    pw = config.get('MERRA', 'password')

    if not os.path.isdir(filedir):
        os.makedirs(filedir)
        print('create foler: {}'.format(filedir))

    flist = []
    for i in range(len(bdate)):
        date = bdate[i]
        filename = '%s/merra-%s-%s.nc4' %(filedir,date,hr)
        flist.append(filename)
        print('Downloading %d of %d: %s'%((i+1),len(bdate),filename))
        yr = date[0:4]
        mon = date[4:6]
        hr = hr
        year = int(yr)
        url1 = 'http://goldsmr5.gesdisc.eosdis.nasa.gov/daac-bin/OTF/'
        url1 += 'HTTP_services.cgi?FILENAME=%2Fdata%2Fs4pa%2FMERRA2%2FM2I6NPANA.5.12.4%2F'
        url2 = '%2F'
        url3 = '%2FMERRA2_300.inst6_3d_ana_Np.'
        url3n = '%2FMERRA2_400.inst6_3d_ana_Np.'
        url3o = '%2FMERRA2_100.inst6_3d_ana_Np.'

        url4 = '.nc4&FORMAT=bmM0Yy8&BBOX=-90%2C-180%2C90%2C180&TIME=1979-01-01T'

        url5 = '%3A00%3A00%2F1979-01-01T'
        url6 = '%3A00%3A00&LABEL=svc_MERRA2_300.inst6_3d_ana_Np.'
        url6n = '%3A00%3A00&LABEL=svc_MERRA2_400.inst6_3d_ana_Np.'
        url6o = '%3A00%3A00&LABEL=svc_MERRA2_100.inst6_3d_ana_Np.'

        url7 = '.nc4&FLAGS=&SHORTNAME=M2I6NPANA&SERVICE=SUBSET_MERRA2&LAYERS=&VERSION=1.02&VARIABLES='

        if year < 1992:
            weburl = '%s%s%s%s%s%s%s%s%s%s%s%s%s' %(url1,yr,url2,mon,url3o,date,url4,hr,url5,hr,url6o,date,url7)
        elif year < 2001:
            weburl = '%s%s%s%s%s%s%s%s%s%s%s%s%s' %(url1,yr,url2,mon,url3,date,url4,hr,url5,hr,url6,date,url7)
        else:
            weburl = '%s%s%s%s%s%s%s%s%s%s%s%s%s' %(url1,yr,url2,mon,url3n,date,url4,hr,url5,hr,url6n,date,url7)
        dir = '%s' %(filename)
                
        if not os.path.exists(dir):
            #urllib.urlretrieve(weburl,dir)
            dloadCmd = 'wget "{}" --user {} --password {} -O {}'.format(weburl, user, pw, filename)
            dloadCmd += ' --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies'
            dloadCmd += ' --auth-no-challenge=on --keep-session-cookies'
            os.system(dloadCmd)

    return flist


def NARRdload(bdate,hr,filedir):
    if not os.path.isdir(filedir):
        os.makedirs(filedir)
        print('create foler: {}'.format(filedir))

    flist = []      
    for k in range(len(bdate)):
            day = bdate[k]
            webdir = day[0:6]
            fname = 'narr-a_221_%s_%s00_000.grb'%(day,hr)
            flist.append('%s/%s'%(filedir,fname))
            weburl='http://nomads.ncdc.noaa.gov/data/narr/%s/%s/%s'%(webdir,day,fname)
            dname = '%s/%s'%(filedir,fname)
            print('Downloading %d of %d: %s'%(k+1,len(bdate),fname))
            if not os.path.exists(dname):
                    urlretrieve(weburl,dname) #,reporthook)

    return flist

