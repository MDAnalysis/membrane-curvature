

"""
Location of datafiles for Membrane Curvature unit tests
=======================================================

MD simulations files stored in `data` sub-directory.

    from membrane_curvature.datafiles import *
"""

__all__ = [
    # From lower to higher complexity
    "GRO_PO4_SMALL",  # Gromacs file of PO4 beads in POPC POPE membrane with 10 lipids
    "XTC_PO4_SMALL",  # Gromacs traj file of PO4 beacs in POPC POPE membrane with 10 lipids with indexes inverted.
    "MEMB_GRO",  # Gromacs file patch of big membrane in squared cell of ~ 240 A dimensions
    "MEMB_XTC",  # Gromacs trajectory file of big membrane.
    "GRO_PO4",  # Gromacs file of PO4 beads in POPC POPE membrane with 914 lipids
    "XTC_PO4",  # Gromacs trajectory of GRO_PO4
    "GRO_MEMBRANE_PROTEIN",  # Gromacs file of POPC POPE CHOL membrane
    "XTC_MEMBRANE_PROTEIN"  # Gromacs trajectory of 10 frames.
    "XTC_MEMBPROT_FIT",  # Gromacs trajectory with rotational and translation fit
    "GRO_MEMBPROT_FIT"  # Gromacs coordinates to load trajectory with fit
]

from pkg_resources import resource_filename

# Membrane protein systems
GRO_MEMBRANE_PROTEIN = resource_filename(__name__, '../data/test_curvature_abca1.gro')
XTC_MEMBRANE_PROTEIN = resource_filename(__name__, '../data/test_curvature_abca1.xtc')
# PO4 beads only
GRO_PO4 = resource_filename(__name__, '../data/test_curvature_po4_only.gro')
XTC_PO4 = resource_filename(__name__, '../data/test_curvature_po4_only.xtc')
# big systems
GRO_PO4_SMALL = resource_filename(__name__, '../data/test_po4_small.gro')
XTC_PO4_SMALL = resource_filename(__name__, '../data/test_po4_small.xtc')
# membrane-only
MEMB_GRO = resource_filename(__name__, '../data/MEMB_traj_short.gro')
MEMB_XTC = resource_filename(__name__, '../data/MEMB_traj_short.xtc')
# membrane-protein
GRO_MEMBPROT_FIT = resource_filename(__name__, '../data/Membrane_protein_fit.gro')
XTC_MEMBPROT_FIT = resource_filename(__name__, '../data/Membrane_protein_fit.xtc')

del resource_filename
