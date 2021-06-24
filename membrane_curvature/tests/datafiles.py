

"""
Location of datafiles for Membrane Curvature unit tests
=======================================================

MD simulations files stored in `data` sub-directory.

    from membrane_curvature.datafiles import *
"""

__all__ = [
    "GRO_MEMBRANE_PROTEIN",  # Gromacs file of POPC POPE CHOL membrane
    "XTC_MEMBRANE_PROTEIN"  # Gromacs trajectory of 10 frames.
    "GRO_PO4_SMALL",  # Gromacs file of PO4 beads in POPC POPE membrane with 10 lipids
    "XTC_PO4_SMALL",  # Gromacs traj file of PO4 beacs in POPC POPE membrane with 10 lipids with indexes inverted.
]

from pkg_resources import resource_filename

GRO_MEMBRANE_PROTEIN = resource_filename(__name__, '../data/test_curvature_abca1.gro')
GRO_PO4 = resource_filename(__name__, '../data/test_curvature_po4_only.gro')
XTC_MEMBRANE_PROTEIN = resource_filename(__name__, '../data/test_curvature_abca1.xtc')
XTC_PO4 = resource_filename(__name__, '../data/test_curvature_po4_only.xtc')
GRO_PO4_SMALL = resource_filename(__name__, '../data/test_po4_small.gro')
XTC_PO4_SMALL = resource_filename(__name__, '../data/test_po4_small.xtc')

del resource_filename
