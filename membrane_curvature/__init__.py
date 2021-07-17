"""

MDAkit for Membrane Curvature
"""
#
# Released under the GNU Public Licence v3 (GPLv3)
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#

# Add imports here
from .surface import normalized_grid, derive_surface, get_z_surface
from .curvature import mean_curvature, gaussian_curvature
from .base import MembraneCurvature

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
