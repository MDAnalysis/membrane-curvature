"""
Unit and regression test for the membrane_curvature package.
"""

# Import package, test suite, and other packages as needed
import pytest
import itertools as it
from ..lib.mods import core_fast_leaflet, curvature, mean_curvature, gaussian_curvature, grid_map
import numpy as np
from numpy.testing import assert_almost_equal
import MDAnalysis as mda
from membrane_curvature.tests.datafiles import (GRO_PO4_SMALL, XTC_PO4_SMALL)

# Reference data from datafile
MEMBRANE_CURVATURE_DATA = {
    'z_avg':{
        'small': np.array([[10., 11., 12.], [15., 20, 25], [15., 40, 50]]) },

    'gaussian_curvature': {
        'small': np.array([[-2.19478738e-02, -2.32254318e-03, -5.47176909e-04],
                           [-1.38453218e-01, -1.21945074e-03, -1.35208221e-04],
                           [-9.72884280e-04, -3.94840040e-04, -1.32808172e-04]])},
    
    'mean_curvature': {
        'small': np.array([[0.16037507, 0.04033506, 0.02057139],
                           [0.99647932, 0.14502529, 0.04590719],
                           [0.05019947, 0.26763173, 0.16841648]]) },

    'grid': {'small':
             {'upper': np.array([[15., 15., np.nan], [15., np.nan, np.nan], [15., np.nan, np.nan]]),
              'lower': np.array([[np.nan, np.nan, 12.], [np.nan, 12., 12.], [np.nan, 12., 12.]])}},

    'beads': {'small':
              {'upper': {'POPC': [0, 1, 2, 3]},
               'lower': {'POPC': [4, 5, 6, 7, 8]}}}

}

def test_gaussian_curvature():
    K_test = gaussian_curvature(MEMBRANE_CURVATURE_DATA['z_avg']['small'])
    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_curvature']['small'], K_test):
        assert_almost_equal(k, k_test)


def test_mean_curvature():
    H_test = mean_curvature(MEMBRANE_CURVATURE_DATA['z_avg']['small'])
    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_curvature']['small'], H_test):
        assert_almost_equal(h, h_test)


def test_core_fast_leaflets():
    n_cells, max_width = 3, 3
    z_calc = np.zeros([n_cells, n_cells])
    u = mda.Universe(GRO_PO4_SMALL, XTC_PO4_SMALL)
    selection = u.select_atoms('index 0:3')
    core_fast_leaflet(u, z_calc, n_cells, selection, max_width)
    for z, z_test in zip(MEMBRANE_CURVATURE_DATA['grid']['small']['upper'], z_calc):
        print(z, z_test)
        assert_almost_equal(z, z_test)

@pytest.mark.parametrize('dummy_coordinates, test_mapper, n_cells, max_width', [(
    # dummy coordinates (x,y)
    ((0, 0), (1, 0), (2, 0),
     (0, 1), (1, 1), (2, 1),
     (0, 2), (1, 2), (2, 2)),
    # should map to
    lambda xy: (xy[0], xy[1]),
    3, 3),
])
def test_grid_map_small_9grid(dummy_coordinates, test_mapper, n_cells, max_width):
    factor = np.float32(n_cells / max_width)
    for dummy_coord in dummy_coordinates:
        assert test_mapper(dummy_coord) == grid_map(dummy_coord, factor)

@pytest.mark.parametrize('dummy_coordinates, test_mapper, n_cells, max_width', [(
    # dummy coordinates (x,y)
    ((0, 0), (1, 0), (2, 0), (3, 0), (4, 0),
     (0, 1), (1, 1), (2, 1), (3, 1), (4, 1),
     (0, 2), (1, 2), (2, 2), (3, 2), (4, 2),
     (0, 3), (1, 3), (2, 3), (3, 3), (4, 3),
     (0, 4), (1, 4), (2, 4), (3, 4), (4, 4)),
    # should map to
    lambda xy: (xy[0], xy[1]),
    5, 5),
])
def test_grid_map_25grid(dummy_coordinates, test_mapper, n_cells, max_width):
    factor = np.float32(n_cells / max_width)
    for dummy_coord in dummy_coordinates:
        assert test_mapper(dummy_coord) == grid_map(dummy_coord, factor)

