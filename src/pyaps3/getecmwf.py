############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
#!/usr/bin/python
from ecmwf import ECMWFDataServer

def getfiles(bdate,hr,filedir,model,humidity='Q'):

        server = ECMWFDataServer('https://api.ecmwf.int/v1',
                'yourkey',
                'youremailadress')

    # Define dataset type
    assert datatype in ('fc','an'), 'Unknown dataset type field for ECMWF'
    if datatype in 'fc':
        dstep = '6'
    elif datatype in 'an':
        dstep = '0'

    # Define grid size according to weather model
    assert model in ('era5','interim', 'hres'), 'Unknown model for ECMWF'
    if model in 'era5':
        gridsize = '0.3/0.3'
    elif model in 'hres':
        gridsize = '0.1/0.1'
    elif model in 'interim':
        gridsize = '0.75/0.75'

    # Define humidity param
    assert humidity in ('Q','R'), 'Unknown humidity field for ECMWF'
    if humidity in 'Q':
        humidparam = 133
    elif humidity in 'R':
        humidparam = 157

    for day in bdate:
        server.retrieve({
          'dataset'  : "%s"%(model),
          'type'     : "%s"%(datatype),
          'date'     : "%s"%(bdate),
          'time'     : "%s"%(hr),
          'step'     : "%s"%(dstep),
          'levtype'  : "pl",
          'levelist' : "all",
          'grid'     : "%s"%(gridsize),
          'param'    : "129/130/%d"%(humidparam),
          'target'   : "%s/%s_%s_%s_%s.grb"%(fileloc,model,bdate,hr,datatype),
        })
