"""
core.py
MDAkit for Membrane Curvature

Handles the primary functions
"""

__author__ = "Estefania Barreto-Ojeda"
version = 0.1

from pathlib import Path
import math
import MDAnalysis as mda
import mdtraj as md
import argparse
import time
import os
import sys
sys.path.append('lib/')
from mods import *

class Grid():
    # instance variables
    def __init__(self, box_width, max_width, unit_cell_width, skip):
        self.box_width = box_width
        self.max_width = max_width
        self.unit_cell_width = unit_cell_width 
        self.skip = skip
        self.n_cells = math.ceil( max_width / unit_cell_width *10 )


class head_indexes:
    # class variables
    leaflets = ['lower', 'upper']

    # instance variables
    def __init__(self, lipid_types, head_index, top):
        self.lipid_types = lipid_types
        self.head_index = head_index
        self.top = top
    
    #@classmethod
    def list_head_beads(self):
        lfs = self.leaflets
        return def_all_beads(self.lipid_types, lfs, self.head_index, self.top)

def main():

    start_time = time.time()

    # 1. Parse input files
    parser = input_options()

    args = parser.parse_args()
    grofile, trjfile = args.f, args.x  
    unit_width, jump = args.uw, args.sk
    output = args.out 

    try:
        name_user = str(args.name)
        print("Prefix assigned by user. Set as {}".format(name_user) )
        
    except NameError:
        name_user = 'system'
        print("Prefix not assigned by user. Set as {}".format(name_user) )
        

    # Define index leaflets
    ii, jj = parse_range(args.io), parse_range(args.ii)
    iil, iiu = def_range_leaflets(args.io, 1)
    jjl, jju = def_range_leaflets(args.ii, 2)


    # 2 --- Set OUTPUT files
    # 2.1 Create output folder
    os.makedirs(os.path.dirname( 'output/' ), exist_ok = True)

    # 2.2 -- Pickle files (head group height)
    prefix = name_user 
    name_ = output + prefix

    ## 3. Define leaflets and lipid type
    leaflets = [ 'lower', 'upper' ]

    #lipid_types = [ ]
    lipid_types = [ 'POPC', 'POPE' ]

    # 4. head indexes
    head_index = [ iil, jjl, jju ]
   
    # 5. Assign topology
    topology = md.load(grofile).topology

    # 6. Populate universe with coordinates and trajectory
    u = mda.Universe(grofile, trjfile)


    # 6.1 Set grid: Extract box dimension from MD sim, 
    # set grid max width, set number of unit cells
    box_size = u.dimensions[0]
    max_width = box_size*0.1
    n_cells = math.ceil( max_width / unit_width *10 )

    # 7. Assign lipids in upper and lower leaflet
    lipid_po4_beads = def_all_beads(lipid_types, leaflets, 
                                    head_index, topology)

    # 8. Load trajectory and grofile using mdtraj
    traj = md.load(trjfile, top=grofile)

    # 9. Save pickles zpo4
    dict_z_coords = core_fast(traj, jump, n_cells, leaflets, lipid_types, 
                         lipid_po4_beads, box_size, max_width, name_)

    # 10. Calculate curvature
    K, H = curvature(dict_z_coords, leaflets, n_cells)

    dict2pickle( name_ + '_H', H)
    dict2pickle( name_ + '_K', K)

    timer(time.time(), start_time)
   
    return 


if __name__ == '__main__':
    main()

