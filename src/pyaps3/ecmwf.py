############################################################
# Program is part of PyAPS                                 #
# Copyright 2012, by the California Institute of Technology#
# Contact: earthdef@gps.caltech.edu                        #
# Modified by A. Benoit and R. Jolivet 2019                #
# Ecole Normale Superieure, Paris                          #
# Contact: insar@geologie.ens.fr                           #
############################################################
import urllib.request
import urllib.parse
import time
import datetime

class ECMWFDataServer:
    def __init__(self,portal,token,email):
        self.version = '0.3'
        self.portal  = portal
        self.token   = token
        self.email   = email

    def _call(self,action,args):

        params = {'_token'   : self.token,
                  '_email'   : self.email,
                  '_action'  : action,
                  '_version' : self.version}
        params.update(args)
            
        data = urllib.parse.urlencode(params)
        req = urllib.request.Request(self.portal, data)
        response = urllib.request.urlopen(req)

        json = response.read();
        json = eval(json)

        if json != None:
            if 'error' in json:
                raise RuntimeError(json['error'])
            if 'message' in json:
                self.put(json['message'])

        return json

    def put(self,*args):
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')),
        for a in args:
            print(a),
        print()

        
    def retrieve(self,args):
        self.put("ECMWF data server batch tool version",self.version);
        user = self._call("user_info",{});
        self.put("Welcome to",user['name'], "from" , user['organisation']);

        r = self._call('retrieve',args)
        rid = r['request']

        last  = ''
        sleep = 0
        while r['status'] != 'complete' and r['status'] != 'aborted':
            text = r['status'] + '.'
            if 'info' in r and r['info'] != None:
                text = text + ' ' + r['info'] 

            if text != last:
                self.put("Request",text)
                last = text

            time.sleep(sleep)
            r = self._call('status',{'request':rid})
            if sleep < 60:
                sleep = sleep + 1

        if r['status'] != last:
            self.put("Request",r['status'])

        if 'reason' in r:
            for m in r['reason']:
                self.put(m)

        if 'result' in r:
            size = long(r['size'])
            self.put("Downloading",self._bytename(size))
            done = self._transfer(r['result'],args['target'])
            self.put("Done")
            if done != size:
                raise RuntimeError("Size mismatch: " + str(done) + " and " + str(size))

        self._call('delete',{'request':rid})
        
        if r['status'] == 'aborted':
            raise RuntimeError("Request aborted")
        
    def _transfer(self,url,path):
        result =  urllib.request.urlretrieve(url,path)
        return long(result[1]['content-length'])

    def _bytename(self,size):   
        next = {'':'K','K':'M','M':'G','G':'T','T':'P'}
        l    = ''
        size = size*1.0
        while 1024 < size:
            l = next[l]
            size = size / 1024
        return "%g %sbyte%s" % (size,l,'s')

