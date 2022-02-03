'''
Does monica control

Created on Jan 4, 2012

@author: ban115
'''

import socket
import types
import threading
import datetime
from time import mktime

type_to_monica_map = {int:'int'}

def type_to_monica(value):
    montype = type_to_monica_map.get(type(value))
    if montype is None:
        raise NotImplementedError("Can't convert python type to monica type: %s" % type(value))
    
    return (montype, str(value))

def get_dutc(dutc=None):
    """Returns DUTC in seconds as an integer
    if dutc is None, 34 it returns 34 seconds, which was valid at 23 Jan 2012
    if dutc is an instance of a MonicaClient, it reutnrs the value of mpclock.misc.clock.dUTC
    Otherwise it returns the input value
    """
    if dutc is None:
        return 34
    elif isinstance(dutc, MonicaClient):
        return int(dutc.poll("mpclock.misc.clock.dUTC").value)
    else:
        return dutc

def get_bat(dt=None, dutc=None):
    """Returns the BAT for the given datetime
    If dt is None, returns the Bat for the current UTC time
    If dutc used to get the value of dUTC. See get_dutc() for details.
    
    
    The time is a BAT, expressed as microseconds since MJD = 0. The Java code which takes the current time and returns a valid BAT is:

  private long timeNow()
  {
    // The value of 3506716800000000L is the BAT as at midnight on
    // 1-Jan-1970, which is the base of the time that the system
    // gives us. It has to be adjusted for leap seconds though.
    return (System.currentTimeMillis() * 1000L) + DUTC.get() * 1000000L + 3506716800000000L;
  }

Here, DUTC.get() returns the value of DUTC in seconds. The current value of DUTC is 34 seconds. 
You could hardcode this value or get the current value from MoniCA by looking at the value of the point mpclock.misc.clock.dUTC.
    """
    if dt is None:
        dt = datetime.datetime.utcnow()
        
    dutc = get_dutc(dutc)
    
    time_mills = int(mktime(dt.timetuple()))*1000
    return (time_mills * 1000) + dutc*1000000 + 3506716800000000
    


class MonicaError(Exception):
    def __init__(self, *args):
        Exception.__init__(self, *args)

class MonicaSetError(MonicaError):
    def __init__(self, ctrlpoint, result):
        Exception.__init__(self, (ctrlpoint, result))
        self.ctrlpoint = ctrlpoint
        self.result = result
        
    def __str__(self):
        return "Monica Error: ctrlpoint: %s result: %s" % (self.ctrlpoint, self.result)
    
class MonicaValue(object):
    def __init__(self, ctrlpoint, bat, value):
        self.ctrlpoint = ctrlpoint
        if bat == '?':
            self.bat = None
        else:
            self.bat = int(bat, 0)
        
        if value == '?':
            self.value = None
        else:
            self.value = value
            
        self.details = None
        
    def set_details(self, d):
        self.details = d
        
    def get_description(self):
        if self.details and 'description' in self.details:
            return self.details.description
        else:
            return self.ctrlpoint
        
    def get_unit(self):
        if self.details:
            return self.details.unit
        else:
            return ''
        
    def __str__(self):
        if self.bat:
            bat = hex(self.bat)
        else:
            bat = '?'
            
        if self.details is None:
            s = '\t'.join((self.ctrlpoint, bat, self.value))
        else:
            s = '%s: %s %s (%s BAT:%s)' % (self.get_description(), 
             self.value, 
             self.details.unit,
             self.ctrlpoint, 
             bat)             

        return s
    
    def __repr__(self):
        return str(self)

def format_values(values):
    return '\n'.join(map(str, values))
    
class MonicaValues(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        
    def __str__(self):
        return format_values(list(self.values()))
        
class MonicaDetails(object):
    def __init__(self, ctrlpoint, samplerate, unit, description):
        self.ctrlpoint = ctrlpoint
        self.sample_rate = float(samplerate)
        self.unit = unit.strip('"')
        self.description = description.strip('"')
        
    def __str__(self):
        return '\t'.join(map(str, (self.ctrlpoint, self.sample_rate, self.unit, self.description)))
    
    def __repr__(self):
        return str(self)
    
class MonicaClient(object):
    '''
    Talks to a monica server.
    
    Mopra (default?) port 8051.
    '''


    def __init__(self, hostport, username='kbmonicapy', password='bar', cache_details=True):
        '''
        Constructor
        '''
        self.hostport = hostport
        self.sock = None
        self.username = username
        self.password = password
        self.bufsize = 4096
        self.details_cache = None
        self._lock = threading.Lock()
        if cache_details:
            self.details_cache = {}
        
    def mset(self, ctrlpoint, value):
        """Sets a single control point to a given value"""
        self.connect()
        
        (montype, monvalue) = type_to_monica(value)
        
        cmd = "set\n%s\n%s\n1\n%s\t%s\t%s\n" % (self.username, self.password, ctrlpoint, montype, monvalue)
        
        with self._lock:
            self.sock.sendall(cmd)
            result = self.sock.recv(self.bufsize)
            
        (ctrlpointr, result) = result.split("\t")
        
        if (ctrlpointr.strip() != ctrlpoint):
            raise MonicaError('Protocol error. Set ctrlpoint(%s) != returned ctrlpoint (%s)' % (ctrlpoint, ctrlpointr))
        
        result = result.strip()
        if result != "OK":
            raise MonicaSetError(ctrlpoint, result)
        
        return result
    
    def _send_cmd(self, cmd, parsefunc, ctrlpoints, values=None):
        """Generic method for getting stuff"""
        self.connect()
        
        if type(ctrlpoints) == bytes:
            points = (ctrlpoints,)
        else:
            points = ctrlpoints 
            
        pointstr = "\n".join(points)
        cmd = "%s\n%d\n%s\n" % (cmd, len(points), pointstr)
        if values is None:
            values = {}
        
        with self._lock:
            self.fsock.write(cmd)
            self.fsock.flush()
    
            # points are returned in the order they were requested
            for point in points:
                line = self.fsock.readline()
                if line.startswith("?"):
                    raise MonicaError("Ctrl point %s does not exist:" % point)
                
                bits = line.split('\t')
                bits[-1] = bits[-1].strip()
                rpoint = bits[0]
                if rpoint != point:
                    raise MonicaError("Protocol error. Recieved ctrlpoint (%s) != expected ctrlpoint(%s)" % (rpoint, point))
                
                try:
                    values[rpoint] = (parsefunc(bits))
                except:
                    raise ValueError("Could not parse result for point: %s. Line='%s'" % (point, line))
            
        if type(ctrlpoints) == bytes:
            return values[rpoint]
        else:
            return values
            
    
    def poll(self, ctrlpoints):
        """
        Polls for the given controlpoints and returns MonicaValues
        
        If ctrlpoints is a string, it returns the value for a single controlpoint.
        
        If ctrlpoints is a list of strings, it will return all the controlpoints as a dictionary, keyed by the controlpoint
        
        Raises a MonicaError if the given controlpoint doesn't exist
        
        """
        values = self._send_cmd('poll', lambda bits: MonicaValue(*bits), ctrlpoints, values=MonicaValues())
        
        if self.details_cache is not None:
            # attach the details
            details = self.details(ctrlpoints)
            if isinstance(values, MonicaValue): # it's a one fff:
                values.details = details
            else: # its returned multpile thngs
                for k,d in details.items():
                    values[k].set_details(d)
        
        return values

    def _array_cmd(self, cmd):
        self.connect()

        with self._lock:
            self.fsock.write(cmd)
            self.fsock.flush()
            line = self.fsock.readline()
            ndata = int(line)

            times = []
            data = []
            for i in range(ndata):
                line = self.fsock.readline()
                bits = line.split('\t')
                times.append(int(bits[0], 16))
                data.append(bits[1])

        return times, data


    def since(self, bat, ctrlpoint):
        '''
        returns MonicaValues for a group of control points sinece a given BAT)
        '''

        cmd = 'since\n{bat} {ctrlpoint}\n'.format(hex(bat), ctrlpoint)
        return self._array_cmd(cmd)

    def between(self, bat1, bat2, ctrlpoint):
        cmd = 'between \n{} {} {}\n'.format(hex(bat1), hex(bat2), ctrlpoint)
        return self._array_cmd(cmd)

    
    def _get_details(self, ctrlpoints):
        """queries for all specified details and saves the details to the cache"""
        if len(ctrlpoints) == 0:
            return {}

        details = self._send_cmd('details', lambda bits: MonicaDetails(*bits), ctrlpoints)
        
        if self.details_cache is not None:
            for k,v in details.items():
                self.details_cache[k] = v
                
        return details

        
    def details(self, ctrlpoints):
        """Returns detaisl for the specified control points
        Will lookup the details cache if required
        """
        # returns something that looks ike:
        #site.environment.weather.Temperature 10.0 "C" "Temperature"
        if isinstance(ctrlpoints, bytes):
            ctrlpoints = (ctrlpoints,)
            
        ctrlset = set(ctrlpoints)
        if self.details_cache is not None:
            known_points = set(self.details_cache.keys())
            unknown_points = ctrlset - known_points
            remaining_points = ctrlset - unknown_points
        else:
            unknown_points = ctrlpoints
            remaining_points = []
            
        # _get_details loads the self.details(unknown_points) into the cache    
        details = self._get_details(unknown_points)

        
        # now attach the known details to the current details
        for p in remaining_points:
            details[p] = self.details_cache[p]
        
        return details
        

    def connect(self):
        if self.sock is None:
            self.sock = socket.create_connection(self.hostport)
            self.fsock = self.sock.makefile()
        
    def close(self):
        if self.sock is not None:
            try:
                self.sock.close()
            except:
                pass
            
    def __del__(self):
        self.close()
        
        
