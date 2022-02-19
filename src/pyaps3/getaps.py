############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
'''Script to automatically download weather model data corresponding to our interferograms.'''

import PyAPS
import tsinsar as ts
import datetime as dt

dates,sat,bperp = ts.load_list('ifg.list')
days,usat,Jmat = ts.ConnectMatrix(dates,sat)


daylist = []
for k in days:
    dobj = dt.date.fromordinal(int(k))
    strobj = '%4d%02d%02d'%(dobj.year,dobj.month,dobj.day)
    daylist.append(strobj)

PyAPS.ECMWFdload(daylist,'18','./Atmos/ECMWF/')

#PyAPS.NARRdload(daylist,'18','./Atmos/NARR/')
#PyAPS.ERAdload(daylist,'18','./Atmos/ERA')
