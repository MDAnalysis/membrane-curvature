"""
core.py
MDAkit for Membrane Curvature

Handles the primary functions
"""

from .lib.mods import dict2pickle, curvature, core_fast_leaflet
import sys
import time
import MDAnalysis as mda
import math
from pathlib import Path
__author__ = "Estefania Barreto-Ojeda"
version = 0.1

sys.path.append('lib/')


def main():

    start_time = time.time()

    # 1. Define leaflets and lipid type # Gone in next steps
    leaflets = ['lower', 'upper']


    # 2. Populate universe with coordinates and trajectory
    u = mda.Universe(grofile, trjfile)

    # 3 Set grid: Extract box dimension from MD sim,
    # set grid max width, set number of unit cells
    box_size = u.dimensions[0]
    max_width = box_size * 0.1
    n_cells = math.ceil(max_width / unit_width * 10)

    # 4. Assign lipids in upper and lower leaflet
    selection = u.select_atoms(atoms_upper, atoms_lower)

    # 5. Save pickles zpo4
    z_Ref = np.zeros([n_cells, n_cells])
    dict_z_coords = core_fast_leaflet(u, z_Ref, n_cells, selection, max_width)

    # 6. Calculate curvature
    K, H = curvature(dict_z_coords, leaflets, n_cells)

    dict2pickle(name_ + '_H', H)
    dict2pickle(name_ + '_K', K)

    timer(time.time(), start_time)

    return


if __name__ == '__main__':
    main()
