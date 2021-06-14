#from collections import OrderedDict
#from itertools import permutations
import itertools as it
import mdtraj as md
import numpy as np
import argparse
import pickle
import textwrap
import time
import math
import os
import re


class Grid:

    """
    A class to set a grid with n x n elements.

    ...

    Attributes
    ----------
    [Names are self explanatory]

    box_width : <float>
        Box dimensions in x.
    max_width : <int>
        Max width of the grid.
    unit_cell_width(): <int>
        Width of unit cell in grid.


    Methods
    -------
    list_by_index(additional=""):
        Returns list with [inner_leaflet_first_index, upper_leaflet_first_index, upper_leaflet_last_index].
    """

    # instance variables
    def __init__(self, box_width, max_width, unit_cell_width):
        self.box_width = box_width
        self.max_width = max_width
        self.unit_cell_width = unit_cell_width
        self.n_cells = math.ceil(max_width / unit_cell_width * 10)


class define_leaflets:

    """
    A class to define leaflets by index when provided by the user.

    ...

    Attributes
    ----------
    int_inner : str:str
        Index interval of the inner leaflet.
    surname : str:str
        Index interval of the outer leaflet.

    Methods
    -------
    list_by_index():
        Returns list [i for i in range(<index interval>)]].
    """

    def __init__(self, int_inner, int_upper):
        self.int_inner = int_inner
        self.int_upper = int_upper

    def list_by_index(self):
        ii, jj = parse_range(self.int_upper), parse_range(self.int_inner)
        iil, iiu = def_range_leaflets(self.int_upper, 1)
        jjl, jju = def_range_leaflets(self.int_inner, 2)

        return [iil, jjl, jju]


class membrane_index:

    """
    A class to select PO4 beads from defined leaflets.

    ...

    Attributes
    ----------
    lipid_types : <list>
        Index interval of the inner leaflet.
    head_index : <>
        Index interval of the outer leaflet.

    Methods
    -------
    list_head_beads():
        Returns list with [inner_leaflet_first_index, upper_leaflet_first_index, upper_leaflet_last_index].
    """

    # class variables
    leaflets = ['upper', 'lower']

    # instance variables
    def __init__(self, lipid_types, head_index, top):
        self.lipid_types = lipid_types
        self.head_index = head_index
        self.top = top

    # #@classmethod
    def list_head_beads(self):
        return def_po4_beads(self.lipid_types, self.leaflets, self.head_index, self.top)


def parse_range(astr):
    """
    Split range provided by user
    """

    rangeG = set()
    x = [int(i) for i in astr.split(':')]
    rangeG.update(range(x[0], x[1] + 1))

    return list(rangeG)


def def_range_leaflets(nn, ind):
    """
    Define range of index residues
    """
    mm = parse_range(nn)

    return min(mm), max(mm)


def def_all_beads(lipid_types, leaflets, head_list, topology):
    """
    Select reference elements to derive membrane surface.


    Parameters
    ----------
    lipid_types : list.
        List of lipid types in simulation box.
    leaflets : str. Default ["lower", "upper"]
        Leaflets of bilayer.
    head_list : list
        List of indexes by leaflet. Provided by user in -ii and -io.
    topology: obj
        Topology from grofile.
    Returns
    -------
    [ element_1, ..., element_n] : list
        List with elements used as a reference to derive membrane surface.

    """

    dic_all_beads = {lf: {lt: [] for lt in lipid_types} for lf in leaflets}
    print('==== Lipid types in membrane ==== ')
    for lt in lipid_types:
        print('====>', lt)
        dic_all_beads['upper'][lt] = np.concatenate((topology.select('resname ' +
                                                                     lt +
                                                                     ' and index ' +
                                                                     str(head_list[0]) +
                                                                     ' to ' +
                                                                     str(head_list[1]) +
                                                                     ' and name PO4'), topology.select('resname ' +
                                                                                                       lt +
                                                                                                       ' and index ' +
                                                                                                       str(head_list[0]) +
                                                                                                       ' to ' +
                                                                                                       str(head_list[1]) +
                                                                                                       ' and name GM1'))).astype(int).tolist()
        dic_all_beads['lower'][lt] = np.concatenate((topology.select(
            'resname ' + lt + ' and index ' + str(head_list[1] + 1) + ' to ' + str(head_list[2]) + ' and name PO4'),
            topology.select(
            'resname ' + lt + ' and index ' + str(head_list[1] + 1) + ' to ' + str(head_list[2]) + ' and name GM1'))).astype(int).tolist()

        print("upper", len(dic_all_beads['upper'][lt]))
        print("lower", len(dic_all_beads['lower'][lt]))

    return dic_all_beads


def dict2pickle(name, dict_):
    """
    Exports values stored in a dictionary to a pickle file.

    Parameters
    ----------
    name : str,
        name of pickle file.
    dict_ : dict like {index -> value}
        dictionary to export.

    Returns
    -------
        Returns a pickled dictionary.

    """

    with open(name + ".pickle", 'wb') as pk:
        pickle.dump(dict_, pk, protocol=pickle.HIGHEST_PROTOCOL)


def core_fast(traj, jump, n_cells, leaflets, lipid_types, lipid_ref,
              box_size, max_width, prefix):
    """
    Runs core_fast_leaflet for each leaflet

    Parameters
    ----------
    traj: trajectory.
        MD trajectory to analyze.
    jump: int. Default 1.
        Skip <jump> number of frames in calculation.
    n_cells : int.
        number of cells in the grid of size `max_width`.
    leaflets : str. Default ["lower", "upper"]
        Leaflets of bilayer.
    lipid_types : list.
        List of lipid types in simulation box.
    lipid_Ref : list
        List of elements used as reference to define membrane surface.
    box_size : float
        Size of box (x dimension).
    max_width : int.
        Maximum width of simulation box.
    prefix : str. Default `None`
        Name of pickle
    Returns
    -------

    Returns a dictionary with name of leaflet as keys. The resulting values
    are nested dictionaries with tuples of [i,j] as keys, where
    0<i<`max_width`, and 0<j<`max_width`. Each key [i,j] has a list of
    `len=(n_frames)` containing the extracted z_coordinate a given [i,j]
    cell in each frame.

    """
    z_ref = {key1: np.zeros([n_cells, n_cells]) for key1 in leaflets}

    for leaflet in leaflets:
        core_fast_leaflet(z_ref[leaflet], leaflet, traj, jump, n_cells,
                          lipid_types, lipid_ref, max_width)

    dict2pickle(prefix, z_ref)

    return z_ref


def core_fast_leaflet(z_Ref, leaflet, traj, jump, n_cells, lipid_types,
                      lipid_ref, max_width):
    """
    Runs core_fast_leaflet for each leaflet

    Parameters
    ----------
    z_Ref : dict.
        Dictionary to store values in every iteration over frames.
    leaflet : str. {"lower", "upper"}
        Leaflets of bilayer.
    traj: trajectory.
        MD trajectory to analyze.
    jump: int. Default 1.
        Skip <jump> number of frames in calculation.
    n_cells : int.
        number of cells in the grid of size `max_width`.
    lipid_types : list.
        List of lipid types in simulation box.
    lipid_ref : list
        List of reference selection for atoms that define the surface
        of the membrane.
    box_size : float
        Size of box (x dimension).
    max_width : int.
        Maximum width of simulation box.

    Returns
    -------

    Returns a dictionary with tuples of [i,j] as keys, where 0<i<`max_width`,
    and 0<j<`max_width`. The resulting values for each key are a list of
    `len=(n_frames)` containing the extracted z_coordinate a given [i,j] cell
    in each frame.

    """

    grid_count = np.zeros([n_cells, n_cells])

    for frame in range(0, traj.n_frames, jump):
        grid_1 = np.zeros([n_cells, n_cells])
        grid_2 = np.zeros([n_cells, n_cells])

        factor = np.float32(n_cells / max_width)

        for lipid_type in lipid_types:

            for bead in lipid_ref[leaflet][lipid_type]:

                x, y, z = traj.xyz[frame, bead, :]

                l = int(abs(x) * factor)
                m = int(abs(y) * factor)

                try:
                    grid_1[l, m] += z
                    grid_2[l, m] += 1

                except BaseException:
                    pass

        for i, j in it.product(range(n_cells), range(n_cells)):
            if grid_2[i, j] > 0:
                z_Ref[i, j] += grid_1[i, j] / grid_2[i, j]
                grid_count[i, j] += 1

    for i, j in it.product(range(n_cells), range(n_cells)):
        if grid_count[i, j] > 0:
            z_Ref[i, j] /= grid_count[i, j]

        else:
            z_Ref[i, j] = np.nan


def gaussian_curvature(Z):
    """
    Calculate gaussian curvature from Z cloud points.


    Parameters
    ----------
    Z : Numpy array, cloud of points.


    Returns
    -------
    K : 2d-array
        Returns 2-dimensional array object with values of mean curvature.

    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    K = (Zxx * Zyy - (Zxy ** 2)) / (1 + (Zx ** 2) + (Zy ** 2)) ** 2

    return K


def mean_curvature(Z):
    """
    Calculates mean curvature from Z cloud points.


    Parameters
    ----------
    Z : Numpy array, cloud of points.


    Returns
    -------
    H : 2d-array
        Returns 2-dimensional array object with values of gaussian curvature.

    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    H = (Zx**2 + 1) * Zyy - 2 * Zx * Zy * Zxy + (Zy**2 + 1) * Zxx
    H = -H / (2 * (Zx**2 + Zy**2 + 1)**(1.5))

    return H
