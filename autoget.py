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

def ECMWFdload(bdate,hr,filedir,model,datatype,humidity='Q'):
    '''
    ECMWF data downloading.

    Args:
        * bdate     : date to download (str)
        * hr        : hour to download (str)
        * filedir   : files directory (str)
        * model     : weather model ('interim' or 'era5' or 'hres')
        * datatype  : reanalysis data type (an)
        * humidity  : humidity
    '''
    
    #-------------------------------------------
    # Initialize

    # Check data
    assert model in ('era5','interim','hres'), 'Unknown model for ECMWF'

    # Infos for downloading
    if model in 'interim':
        print('WARNING: you are downloading from the old ECMWF platform. ERA-Interim is deprecated, use ERA-5 instead.')
    if model in 'era5':
        print('INFO: You are using the latest ECMWF platform for downloading datasets: https://cds.climate.copernicus.eu/api/v2')

    #-------------------------------------------
    # Define parameters

    # Humidity
    assert humidity in ('Q','R'), 'Unknown humidity field for ECMWF'
    if humidity in 'Q':
        if model in 'era5':
            humidparam = 'specific_humidity'
        else:
            humidparam = 133
    elif humidity in 'R':
        if model in 'era5':
            humidparam = 'relative_humidity'
        else:
            humidparam = 157

    # Grid size (only for HRES and ERA-Interim)
    if model in 'hres':
        gridsize = '0.10/0.10'
    elif model in 'interim':
        gridsize = '0.75/0.75'

    #-------------------------------------------
    # Iterate over dates

    flist = []
    for k in range(len(bdate)):
        day = bdate[k]

        #-------------------------------------------
        # CASE 1: request for CDS API client (new ECMWF platform, for ERA-5)    
        if model in 'era5':
            
            # Contact the server
            c = cdsapi.Client()

            # Pressure levels
            pressure_lvls = ['1','2','3','5','7','10','20','30','50', 
                            '70','100','125','150','175','200','225',
                            '250','300','350','400','450','500','550',
                            '600','650','700','750','775','800','825',
                            '850','875','900','925','950','975','1000']

            # Dictionary
            indict = {'product_type':'reanalysis',
                      'format':'grib',
                      'variable':['geopotential','temperature','{}'.format(humidparam)],
                      'pressure_level': pressure_lvls,
                      'year':'{}'.format(day[0:4]),
                      'month':'{}'.format(day[4:6]),
                      'day':'{}'.format(day[6:8]),
                      'time':'{}:00'.format(hr)}

            # Output
            fname = '%s/ERA-5_%s_%s_%s.grb'%(filedir,day,hr,datatype)
            flist.append(fname)

            # Assert grib file not yet downloaded
            if not os.path.exists(fname):
                print('Downloading %d of %d: %s '%(k+1,len(bdate),fname))
                print(indict)
                
                # Make the request
                c.retrieve(
                        'reanalysis-{}-pressure-levels'.format(model),
                        indict,
                        fname)
        
        #-------------------------------------------
        # CASE 2: request for WEB API client (old ECMWF platform, deprecated, for ERA-Int and HRES)
        else:

            # Contact the server
            from pyaps3.ecmwfapi import ECMWFDataServer
            dpath = os.path.dirname(__file__)
            config = ConfigParser.RawConfigParser()
            config.read(os.path.expanduser('~/.ecmwfcfg'))
            url = "https://api.ecmwf.int/v1"
            emid = config.get('ECMWF', 'email')
            key = config.get('ECMWF', 'key')
            server = ECMWFDataServer(url=url, key=key, email=emid)

            # Output
            if model in 'interim':
                fname = '%s/ERA-Int_%s_%s_%s.grb'%(filedir,day,hr,datatype)
            elif model in 'hres':
                fname = '%s/HRES_%s_%s_%s.grb'%(filedir,day,hr,datatype)
            flist.append(fname)

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
