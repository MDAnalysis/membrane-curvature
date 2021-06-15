

"""
Location of datafiles for Membrane Curvature unit tests
=======================================================

MD simulations files stored in ``data`` sub-directory.

    from mdakit_membcurv.datafiles import *
"""

__all__ = [
    "GRO_",  # Gromacs file of POPC POPE CHOL membrane
    "XTC_"  # Gromacs trajectory of 10 frames.
]

from pkg_resources import resource_filename

GRO_MEMBRANE_PROTEIN = resource_filename(__name__, '../data/test_curvature_abca1.gro')
GRO_PO4 = resource_filename(__name__, '../data/test_curvature_po4_only.gro')
XTC_MEMBRANE_PROTEIN = resource_filename(__name__, '../data/test_curvature_abca1.xtc')
XTC_PO4 = resource_filename(__name__, '../data/test_curvature_po4_only.xtc')
del resource_filename
