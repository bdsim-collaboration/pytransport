"""
A module for converting a TRANSPORT file into gmad for use in BDSIM.
    
To use:
>>> self = pytransport.convert.pytransport()
>>> self.load_file(TRANSPORTfile)
>>> self.convert()
    
Will output:
    
filename.gmad
filename_beam.gmad
filename_components.gmad
filename_options.gmad
filename_samplers.gmad
filename_sequence.gmad
"""

import _General
import elements
import convert
import reader
