#
# User defined tasks setup.
# Generated from buildmytask.
#

import sys
from casa_stack_manip import stack_frame_find

if sys.path[1] != '/Users/ban115/bolton/craft/code/casatasks':
  sys.path.insert(1, '/Users/ban115/bolton/craft/code/casatasks')
from odict import odict
if not globals().has_key('mytasks') :
  mytasks = odict()

mytasks['freqwt'] = 'Adjust WEIGHTS_SPECTRUM and optionally data by a spectrum'

if not globals().has_key('task_location') :
  task_location = odict()

task_location['freqwt'] = '/Users/ban115/bolton/craft/code/casatasks'
myglobals = stack_frame_find( )
tasksum = myglobals['tasksum'] 
for key in mytasks.keys() :
  tasksum[key] = mytasks[key]

from freqwt_cli import  freqwt_cli as freqwt
