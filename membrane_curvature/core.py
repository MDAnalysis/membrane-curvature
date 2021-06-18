"""
core.py
MDAkit for Membrane Curvature

Handles the primary functions
"""

from .lib.mods import dict2pickle, core_fast, curvature, def_all_beads
import sys
import os
import time
import mdtraj as md
import MDAnalysis as mda
import math
from pathlib import Path
__author__ = "Estefania Barreto-Ojeda"
version = 0.1

sys.path.append('lib/')


def main():

    start_time = time.time()

    # 3. Define leaflets and lipid type
    leaflets = ['lower', 'upper']

    #lipid_types = [ ]
    lipid_types = ['POPC', 'POPE']

    # 4. head indexes
    head_index = [iil, jjl, jju]

    # 5. Assign topology
    topology = md.load(grofile).topology

    # 6. Populate universe with coordinates and trajectory
    u = mda.Universe(grofile, trjfile)

    # 6.1 Set grid: Extract box dimension from MD sim,
    # set grid max width, set number of unit cells
    box_size = u.dimensions[0]
    max_width = box_size * 0.1
    n_cells = math.ceil(max_width / unit_width * 10)

    # 7. Assign lipids in upper and lower leaflet
    lipid_ref_beads = def_all_beads(lipid_types, leaflets,
                                    head_index, topology)

    # 8. Load trajectory and grofile using mdtraj
    traj = md.load(trjfile, top=grofile)

    # 9. Save pickles zpo4
    dict_z_coords = core_fast(traj, jump, n_cells, leaflets, lipid_types,
                              lipid_ref_beads, box_size, max_width, name_)

    # 10. Calculate curvature
    K, H = curvature(dict_z_coords, leaflets, n_cells)

    dict2pickle(name_ + '_H', H)
    dict2pickle(name_ + '_K', K)

    timer(time.time(), start_time)

    return


if __name__ == '__main__':
    main()
