"""
pytransport - Royal Holloway utility to manipulate TRANSPORT data and models.

Authors:

 * William Shields
 * Jochem Snuverink

Copyright Royal Holloway, University of London 2023.

"""
try:
    from ._version import version as __version__
    from ._version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")


from . import _General
from . import Compare
from . import Convert
from . import Data
from . import Reader

__all__ = ['Compare',
           'Convert',
           'Data',
           'Reader']
