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

#####Replace with your ECMWF key and email address.
####server = ECMWFDataServer(
####       'http://data-portal.ecmwf.int/data/d/dataserver/',
####       'a9a3c89533c035e5a92c5758fbf27875',
####       'abc@def.com'
####    )


def getfiles(bdate,hr,filedir,humidity='Q'):

        server = ECMWFDataServer('https://api.ecmwf.int/v1',
        'yourkey',
        'youremailadress')

    # Define humidity param
    assert humidity in ('Q','R'), 'Unknown humidity field for ECMWF'
    if humidity in 'Q':
        humidparam = 133
    elif humidity in 'R':
        humidparam = 157

    for day in bdate:
        server.retrieve({
          'class'    : 'od',
          'stream'   : 'oper',
          'expver'   : '1',
          'date'     : "%s"%(bdate),
          'time'     : "%s"%(hr),
          'step'     : "0",
          'type'     : "an",
          'levtype'  : "pl",
          'levelist' : "all",
          'param'    : "129/130/%d"%(humidparam),
          'target'   : "%s/%s_%s_%s.grb"%(fileloc,bdate,hr),
        })
