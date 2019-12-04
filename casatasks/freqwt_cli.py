#
# This file was generated using xslt from its XML file
#
# Copyright 2014, Associated Universities Inc., Washington DC
#
import sys
import os
#from casac import *
import casac
import string
import time
import inspect
import gc
import numpy
from casa_stack_manip import stack_frame_find
from odict import odict
from types import * 
from task_freqwt import freqwt
class freqwt_cli_:
    __name__ = "freqwt"
    rkey = None
    i_am_a_casapy_task = None
    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.

    def __init__(self) :
       self.__bases__ = (freqwt_cli_,)
       self.__doc__ = self.__call__.__doc__

       self.parameters={'vis':None, 'specfile':None, 'weightdata':None, 'cutoff':None, }


    def result(self, key=None):
	    #### and add any that have completed...
	    return None


    def __call__(self, vis=None, specfile=None, weightdata=None, cutoff=None, ):

        """Adjust WEIGHTS_SPECTRUM and optionally data by a spectrum

	Detailed Description: 
Adjust WEIGHTS_SPECTRUM and optionally data by a spectrum.

  freqwt loads the the weights as a text file. The first column is frequency in GHz, second column is an amplitude (w).

  freqwt does a first order interpolation if the input specturm and WEIGHTS_SPECTRUM column are on different frequency grids.

  Values of the interpoltion that are less than the cutoff are set to zero.

  freqwt multiplies WEIGHTS_SPECTRUM column by w[i]**2/sum(w**2).

  If weightdata=True, it also divides the DATA column by w[i], where the weights are non-zero.

  
  
	Arguments :
		vis:	Input measurement set
		   Default Value: 

		specfile:	Spectrum file to weight by
		   Default Value: 

		weightdata:	Divide data by the weights as well
		   Default Value: True

		cutoff:	Set weights to zero if spectrum amplitude is below this number
		   Default Value: 0

	Returns: void

	Example :

    Adjust WEIGHTS_SPECTRUM and optionally data by a spectrum

   Keyword arguments:
   
   vis -- name of input image file
   spec -- text file of spectrum - see output of some function. The first column should be frequency in GHz, the second column should be the amplitude

  
        """
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
	#casac = self.__globals__['casac']
	casalog = self.__globals__['casalog']
	casa = self.__globals__['casa']
	#casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = 'freqwt'
        self.__globals__['taskname'] = 'freqwt'
        ###
        self.__globals__['update_params'](func=self.__globals__['taskname'],printtext=False,ipython_globals=self.__globals__)
        ###
        ###
        #Handle globals or user over-ride of arguments
        #
        if type(self.__call__.func_defaults) is NoneType:
            function_signature_defaults={}
	else:
	    function_signature_defaults=dict(zip(self.__call__.func_code.co_varnames[1:],self.__call__.func_defaults))
	useLocalDefaults = False

        for item in function_signature_defaults.iteritems():
                key,val = item
                keyVal = eval(key)
                if (keyVal == None):
                        #user hasn't set it - use global/default
                        pass
                else:
                        #user has set it - use over-ride
			if (key != 'self') :
			   useLocalDefaults = True

	myparams = {}
	if useLocalDefaults :
	   for item in function_signature_defaults.iteritems():
	       key,val = item
	       keyVal = eval(key)
	       exec('myparams[key] = keyVal')
	       self.parameters[key] = keyVal
	       if (keyVal == None):
	           exec('myparams[key] = '+ key + ' = self.itsdefault(key)')
		   keyVal = eval(key)
		   if(type(keyVal) == dict) :
                      if len(keyVal) > 0 :
		         exec('myparams[key] = ' + key + ' = keyVal[len(keyVal)-1][\'value\']')
		      else :
		         exec('myparams[key] = ' + key + ' = {}')
	 
        else :
            print ''

            myparams['vis'] = vis = self.parameters['vis']
            myparams['specfile'] = specfile = self.parameters['specfile']
            myparams['weightdata'] = weightdata = self.parameters['weightdata']
            myparams['cutoff'] = cutoff = self.parameters['cutoff']


	result = None

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['vis'] = vis
        mytmp['specfile'] = specfile
        mytmp['weightdata'] = weightdata
        mytmp['cutoff'] = cutoff
	pathname="file:///Users/ban115/bolton/craft/code/casatasks/"
	trec = casac.casac.utils().torecord(pathname+'freqwt.xml')

        casalog.origin('freqwt')
	try :
          #if not trec.has_key('freqwt') or not casac.casac.utils().verify(mytmp, trec['freqwt']) :
	    #return False

          casac.casac.utils().verify(mytmp, trec['freqwt'], True)
          scriptstr=['']
          saveinputs = self.__globals__['saveinputs']
          if type(self.__call__.func_defaults) is NoneType:
              saveinputs=''
          else:
              saveinputs('freqwt', 'freqwt.last', myparams, self.__globals__,scriptstr=scriptstr)
          tname = 'freqwt'
          spaces = ' '*(18-len(tname))
          casalog.post('\n##########################################'+
                       '\n##### Begin Task: ' + tname + spaces + ' #####')
          if type(self.__call__.func_defaults) is NoneType:
              casalog.post(scriptstr[0]+'\n', 'INFO')
          else :
              casalog.post(scriptstr[1][1:]+'\n', 'INFO')
          result = freqwt(vis, specfile, weightdata, cutoff)
          casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####'+
                       '\n##########################################')

	except Exception, instance:
          if(self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__['__rethrow_casa_exceptions']) :
             raise
          else :
             #print '**** Error **** ',instance
	     tname = 'freqwt'
             casalog.post('An error occurred running task '+tname+'.', 'ERROR')
             pass

        gc.collect()
        return result
#
#
#
    def paramgui(self, useGlobals=True, ipython_globals=None):
        """
        Opens a parameter GUI for this task.  If useGlobals is true, then any relevant global parameter settings are used.
        """
        import paramgui
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )

        if useGlobals:
	    if ipython_globals == None:
                myf=self.__globals__
            else:
                myf=ipython_globals

            paramgui.setGlobals(myf)
        else:
            paramgui.setGlobals({})

        paramgui.runTask('freqwt', myf['_ip'])
        paramgui.setGlobals({})

#
#
#
    def defaults(self, param=None, ipython_globals=None, paramvalue=None, subparam=None):
	if not hasattr(self, "__globals__") or self.__globals__ == None :
           self.__globals__=stack_frame_find( )
        if ipython_globals == None:
            myf=self.__globals__
        else:
            myf=ipython_globals

        a = odict()
        a['vis']  = ''
        a['specfile']  = ''
        a['weightdata']  = True
        a['cutoff']  = 0


### This function sets the default values but also will return the list of
### parameters or the default value of a given parameter
        if(param == None):
                myf['__set_default_parameters'](a)
        elif(param == 'paramkeys'):
                return a.keys()
        else:
            if(paramvalue==None and subparam==None):
               if(a.has_key(param)):
                  return a[param]
               else:
                  return self.itsdefault(param)
            else:
               retval=a[param]
               if(type(a[param])==dict):
                  for k in range(len(a[param])):
                     valornotval='value'
                     if(a[param][k].has_key('notvalue')):
                        valornotval='notvalue'
                     if((a[param][k][valornotval])==paramvalue):
                        retval=a[param][k].copy()
                        retval.pop(valornotval)
                        if(subparam != None):
                           if(retval.has_key(subparam)):
                              retval=retval[subparam]
                           else:
                              retval=self.itsdefault(subparam)
		     else:
                        retval=self.itsdefault(subparam)
               return retval


#
#
    def check_params(self, param=None, value=None, ipython_globals=None):
      if ipython_globals == None:
          myf=self.__globals__
      else:
          myf=ipython_globals
#      print 'param:', param, 'value:', value
      try :
         if str(type(value)) != "<type 'instance'>" :
            value0 = value
            value = myf['cu'].expandparam(param, value)
            matchtype = False
            if(type(value) == numpy.ndarray):
               if(type(value) == type(value0)):
                  myf[param] = value.tolist()
               else:
                  #print 'value:', value, 'value0:', value0
                  #print 'type(value):', type(value), 'type(value0):', type(value0)
                  myf[param] = value0
                  if type(value0) != list :
                     matchtype = True
            else :
               myf[param] = value
            value = myf['cu'].verifyparam({param:value})
            if matchtype:
               value = False
      except Exception, instance:
         #ignore the exception and just return it unchecked
         myf[param] = value
      return value
#
#
    def description(self, key='freqwt', subkey=None):
        desc={'freqwt': 'Adjust WEIGHTS_SPECTRUM and optionally data by a spectrum',
               'vis': 'Input measurement set',
               'specfile': 'Spectrum file to weight by',
               'weightdata': 'Divide data by the weights as well',
               'cutoff': 'Set weights to zero if spectrum amplitude is below this number',

              }

#
# Set subfields defaults if needed
#

        if(desc.has_key(key)) :
           return desc[key]

    def itsdefault(self, paramname) :
        a = {}
        a['vis']  = ''
        a['specfile']  = ''
        a['weightdata']  = True
        a['cutoff']  = 0

        #a = sys._getframe(len(inspect.stack())-1).f_globals

        if a.has_key(paramname) :
	      return a[paramname]
freqwt_cli = freqwt_cli_()
