"""
MDAkit MembCurv
MDAkit for Membrane Curvature
"""

# Add imports here
from .core import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
