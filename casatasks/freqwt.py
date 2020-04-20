#
# This file was generated using xslt from its XML file
#
# Copyright 2009, Associated Universities Inc., Washington DC
#
import sys
import os
from  casac import *
import string
from taskinit import casalog
from taskinit import xmlpath
#from taskmanager import tm
import task_freqwt
def freqwt(vis='', specfile='', weightdata=True, cutoff=0):

        """Adjust WEIGHTS_SPECTRUM and optionally data by a spectrum
    Adjust WEIGHTS_SPECTRUM and optionally data by a spectrum

   Keyword arguments:
   
   vis -- name of input image file
   spec -- text file of spectrum - see output of some function. The first column should be frequency in GHz, the second column should be the amplitude

  
        """

#
#    The following is work around to avoid a bug with current python translation
#
        mytmp = {}

        mytmp['vis'] = vis
        mytmp['specfile'] = specfile
        mytmp['weightdata'] = weightdata
        mytmp['cutoff'] = cutoff
	pathname="file:///Users/ban115/bolton/craft/code/casatasks/"
	trec = casac.utils().torecord(pathname+'freqwt.xml')

        casalog.origin('freqwt')
        if trec.has_key('freqwt') and casac.utils().verify(mytmp, trec['freqwt']) :
	    result = task_freqwt.freqwt(vis, specfile, weightdata, cutoff)

	else :
	  result = False
        return result
