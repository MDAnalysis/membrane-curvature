"""
core.py
MDAkit for Membrane Curvature

Handles the primary functions
"""

from membrane_curvature.surface import normalized_grid, derive_surface, get_z_surface
from membrane_curvature.curvature import mean_curvature, gaussian_curvature
import time
import MDAnalysis as mda
import math
from pathlib import Path
__author__ = "Estefania Barreto-Ojeda"
version = 0.1


def main():



    # 1. Populate universe with coordinates and trajectory
    u = mda.Universe(topology, trajectory)

    # 2 Set grid: Extract box dimension from MD sim,
    # set grid max width, set number of unit cells
    box_size = u.dimensions[0]
    max_width = box_size * 0.1
    n_cells = math.ceil(max_width / unit_width * 10)

    # 3. Assign lipids in upper and lower leaflet
    selection = u.select_atoms(atoms_upper, atoms_lower)

    # 4. Save pickles zpo4
    z_Ref = np.zeros([n_cells, n_cells])
    z_coords = derive_surface(n_cells, selection, max_width)

    # 5. Calculate curvature
    K = gaussian_curvature(z_coords)
    H = mean_curvature(z_coords)

    return


if __name__ == '__main__':
    main()
