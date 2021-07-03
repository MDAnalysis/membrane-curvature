"""
Unit and regression test for the membrane_curvature package.
"""

# Import package, test suite, and other packages as needed
import pytest
from ..lib.mods import mean_curvature, gaussian_curvature, avg_unit_cell, derive_surface
import numpy as np
from numpy.testing import assert_almost_equal
import MDAnalysis as mda
from membrane_curvature.tests.datafiles import (GRO_PO4_SMALL, XTC_PO4_SMALL)

# Reference data from datafile
MEMBRANE_CURVATURE_DATA = {

    'z_avg_coords': np.array([[15.52998744, 15.68577027, 15.54460264, 15.42630775, 15.41610875,
                               15.5227999, 15.40030259, 15.47866627, 15.5402739, 15.58982381],
                              [15.7210999, 15.66186774, 15.66852082, 15.5650629, 15.59180416,
                               15.52540419, 15.43453586, 15.31625871, 15.53287837, 15.72669662],
                              [15.65202988, 15.87587967, 15.80443607, 15.2156001, 15.9012679,
                               15.61628637, 15.54388155, 15.51714352, 15.56670715, 15.67794791],
                              [15.79848438, 15.86021965, 15.91864504, 15.69624091, 16.07622673,
                               15.55569959, 15.33016631, 15.51812384, 15.57359306, 15.74745874],
                              [15.99946331, 15.80777305, 15.93547516, 16.06895507, 16.49331617,
                               16.47925099, 15.74151387, 15.63836954, 15.72500532, 15.91917461],
                              [16.00306618, 15.44807954, 14.91708903, 15.62196315, 15.63327226,
                               16.56862831, 16.11996038, 15.89331185, 15.83679604, 15.98079948],
                              [15.76326987, 15.64245421, 15.52017543, 15.47039083, 14.88002689,
                               15.53698847, 15.80894958, 15.80674416, 15.7546208, 15.88995425],
                              [15.70029967, 15.61600719, 15.53771198, 15.47668145, 15.19293903,
                               15.41956353, 15.66099601, 15.71207747, 15.80128901, 15.71542822],
                              [15.68171964, 15.63995413, 15.53009812, 15.21636634, 15.27028742,
                               15.38000841, 15.51385563, 15.64464232, 15.79420718, 15.77794963],
                              [11.9375362, 15.67473274, 15.56693593, 15.28816211, 15.26380259,
                               15.34311786, 15.50101777, 15.5245863, 13.5766693, 15.69553947]]),

    'z_ref': {
        'small': np.array([[10., 11., 12.], [15., 20, 25], [15., 40, 50]]),
        'all': np.array([[15.96339989, 16.36299992, 16.57025003, 16.30183315, 16.14233398,
                        15.94885731, 15.99250078, 16.17940025, 16.08424997, 15.93374944],
                         [15.9995718, 16.48649979, 16.59700012, 16.3276666, 16.26959991,
                        15.77983316, 15.67449999, 15.59950042, 16.05650028, 16.11733341],
                         [15.9301672, 16.04720001, 16.24383338, 16.38975, 16.05666653,
                        15.71950006, 15.7414999, 15.65285724, 15.71783352, 15.91666635],
                         [15.87350019, 16.0994997, 16.45200014, 16.38366667, 15.91100025,
                          15.44099998, 15.55220013, 15.74933386, 15.7957499, 15.93225002],
                         [15.89100003, 16.01559982, 16.45950031, 16.68450022, 16.34674978,
                          16.27950001, 15.97475028, 16.0142498, 16.07933331, 15.96724939],
                         [15.81012511, 15.84583362, 16.09700036, 15.98525, 15.49299908,
                          16.36499977, 16.20639992, 15.8682003, 15.82559967, 15.87400055],
                         [15.73475003, 15.67866707, 15.85220013, 15.60228566, 15.12299967,
                          15.70033328, 15.87920036, 15.80550003, 15.60928576, 15.8010006],
                         [15.79039974, 15.91499996, 15.97549987, 15.80860004, 15.73637486,
                          15.51133362, 15.80240021, 15.78233337, 15.65516663, 15.72000027],
                         [15.8411665, 16.13249969, 16.48759995, 16.25674987, 15.78233369,
                          15.71450011, 15.33062541, 15.99500027, 15.83737516, 15.74500052],
                         [15.82900023, 16.06166649, 16.2973334, 16.43733342, 16.12957178,
                          16.09366608, 15.95349979, 16.22599983, 16.17750025, 16.00225067]])},

    'gaussian_curvature': {
        'small': np.array([[-2.19478738e-02, -2.32254318e-03, -5.47176909e-04],
                           [-1.38453218e-01, -1.21945074e-03, -1.35208221e-04],
                           [-9.72884280e-04, -3.94840040e-04, -1.32808172e-04]]),

        'all': np.array([[-0.00189129, 0.04934495, 0.04660196, -0.00351411, -0.01943406,
                          -0.03893078, -0.0246221, -0.03280081, -0.12642188, -0.03525892],
                        [-0.00693942, 0.03182844, 0.00120618, -0.00352255, 0.00342091,
                            -0.00142247, 0.02465936, 0.04605395, -0.0106348, -0.00313678],
                        [-0.01649198, 0.00168855, -0.02731173, -0.01564413, -0.00973393,
                            0.02368517, -0.00604347, 0.0169452, 0.0113627, -0.00027235],
                        [0.00088378, -0.00464509, 0.0075489, 0.02877688, -0.00328288,
                            0.07988983, 0.02487373, -0.00071568, -0.0050775, -0.02188734],
                        [-0.013283, -0.01107469, 0.02580236, 0.0847512, -0.07205011,
                            -0.05251712, -0.03977956, -0.03174133, 0.0151017, 0.00278634],
                        [-0.00688401, -0.01030927, -0.03964658, -0.01066038, 0.01942349,
                            0.01216024, 0.05635031, -0.0138591, -0.00547026, -0.02372161],
                        [0.00552771, 0.00111084, -0.0688243, 0.0081551, 0.14126933,
                            -0.01476609, -0.00715425, -0.0059002, 0.02781192, 0.00485525],
                        [-0.02954713, -0.01709626, -0.01171343, -0.00766876, -0.01780511,
                            -0.04009694, -0.00779307, -0.04735893, 0.00799721, -0.00961057],
                        [-0.0033638, 0.009781, 0.06724205, 0.00795326, 0.00034145,
                            0.02682387, 0.10108879, -0.01423631, -0.01802192, -0.00922448],
                        [-0.00313899, -0.00418259, 0.03487913, -0.04456535, -0.00768992,
                            -0.00642677, 0.0254065, -0.01830984, -0.00904487, -0.01182518]])},

    'mean_curvature': {
        'small': np.array([[0.16037507, 0.04033506, 0.02057139],
                           [0.99647932, 0.14502529, 0.04590719],
                           [0.05019947, 0.26763173, 0.16841648]]),

        'all': np.array([[0.06396822, 0.2254189, 0.22399141, 0.02498465, 0.05090628,
                          -0.09999741, -0.10476551, -0.05358351, 0.12631834, 0.09827992],
                         [0.09151143, 0.19973338, 0.13808081, 0.04574428, 0.10289627,
                          -0.04137634, -0.18304112, -0.2249239, 0.04760685, 0.16729716],
                         [-0.01779114, -0.04597901, -0.00778294, 0.09933345, -0.03402726,
                          -0.16382893, -0.09622266, -0.15120851, -0.13753756, -0.04029748],
                         [-0.0369991, 0.00233226, 0.16181129, 0.18932471, -0.01349983,
                          -0.28269057, -0.16712506, 0.01942325, -0.03927295, -0.03811465],
                         [-0.06309625, -0.00820666, 0.16715942, 0.29301601, 0.21669994,
                          0.07393958, -0.0129063, 0.04319332, 0.14495082, 0.07021294],
                         [-0.05875299, -0.04393517, 0.0837562, -0.05789939, -0.19497179,
                          0.25517884, 0.25387131, -0.03625653, -0.03502722, -0.01136375],
                         [-0.08988405, -0.10299823, -0.04910499, -0.21748747, -0.41463502,
                          -0.03913769, 0.1765791, -0.03554145, -0.1677006, -0.10516181],
                         [0.0095539, 0.0364163, 0.00168944, -0.06394463, -0.04456537,
                          -0.25037463, -0.03814847, -0.02531541, -0.11902046, -0.10177806],
                         [0.00184591, 0.12102329, 0.28902913, 0.0966804, -0.03156109,
                          -0.16964647, -0.336664, 0.03280685, 0.03212416, -0.08340905],
                         [0.01478105, 0.08136444, 0.23413597, 0.1472945, -0.06672947,
                          -0.09468121, -0.21140388, 0.03506431, 0.03308529, -0.01943328]])},

    'grid': {'small':
             {'upper': np.array([[15., 15., np.nan], [15., np.nan, np.nan], [15., np.nan, np.nan]]),
              'lower': np.array([[np.nan, np.nan, 12.], [np.nan, 12., 12.], [np.nan, 12., 12.]])}},

    'beads': {'small':
              {'upper': {'POPC': [0, 1, 2, 3]},
               'lower': {'POPC': [4, 5, 6, 7, 8]}}}

}


@pytest.fixture
def small_grofile():
    u = mda.Universe(GRO_PO4_SMALL)
    sel = u.select_atoms('name PO4')
    return sel


def test_gaussian_curvature_small():
    K_test = gaussian_curvature(MEMBRANE_CURVATURE_DATA['z_ref']['small'])
    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_curvature']['small'], K_test):
        assert_almost_equal(k, k_test)


def test_mean_curvature_small():
    H_test = mean_curvature(MEMBRANE_CURVATURE_DATA['z_ref']['small'])
    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_curvature']['small'], H_test):
        assert_almost_equal(h, h_test)


def test_gaussian_curvature_all():
    K_test = gaussian_curvature(MEMBRANE_CURVATURE_DATA['z_ref']['all'])
    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_curvature']['all'], K_test):
        assert_almost_equal(k, k_test)


def test_mean_curvature_all():
    H_test = mean_curvature(MEMBRANE_CURVATURE_DATA['z_ref']['all'])
    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_curvature']['all'], H_test):
        assert_almost_equal(h, h_test)


@pytest.mark.parametrize('n_cells, grid_z_coords, unit', [(
    # number of cells
    3,
    # z position in grid
    np.array([[10., 10., 10.],
              [10., 10., 10.],
              [10., 10., 10.]]),
    # number of beads per grid
    [[1., 1., 1.],
     [1., 1., 1.],
     [1., 1., 1.]])])
def test_avg_unit_cell_identity(n_cells, grid_z_coords, unit):
    z_ref = np.zeros([n_cells, n_cells])
    unit = np.ones([n_cells, n_cells])
    unit_cell = avg_unit_cell(z_ref, grid_z_coords, unit)
    for z_position, cell in zip(unit_cell, range(n_cells)):
        assert z_position[cell] == 10.


@pytest.mark.parametrize('n_cells, grid_z_coords', [(
    # number of cells
    3,
    # z position in grid
    np.array([[10., 20., 30.],
              [10., 20., 30.],
              [10., 20., 30.]]
             ))])
def test_avg_unit_cell_identity_other_values(n_cells, grid_z_coords):
    z_ref = np.zeros([n_cells, n_cells])
    unit = np.ones([n_cells, n_cells])
    z_avg = avg_unit_cell(z_ref, grid_z_coords, unit)
    for cell in range(n_cells):
        assert_almost_equal(z_avg[cell], grid_z_coords[cell])


@pytest.mark.parametrize('n_cells, grid_z_coords, grid_norm, grid_avg', [(
    # number of cells
    3,
    # sum of z coordinate in grid
    np.array([[10., 10., 10.],
              [10., 10., 10.],
              [10., 10., 10.]]),
    # grid number of beads per unit
    np.array([[2., 1., 1.],
              [1., 2., 1.],
              [1., 1., 2.]]),
    # avg z coordinate in grid
    np.array([[5., 10., 10.],
              [10., 5., 10.],
              [10., 10., 5.]]),
)])
def test_avg_unit_cell_more_beads(n_cells, grid_z_coords, grid_norm, grid_avg):
    z_ref = np.zeros([n_cells, n_cells])
    unit_cell = avg_unit_cell(z_ref, grid_z_coords, grid_norm)
    for cell in range(n_cells):
        assert_almost_equal(unit_cell[cell], grid_avg[cell])


@pytest.mark.parametrize('n_cells, max_width, z_avg', [(
    # number of cells
    3, 3,
    # z ref of surface
    # z coordinate in grid
    np.array(([150., 150., 120.],
              [150., 120., 120.],
              [150., 120., 120.])))])
def test_derive_surface(small_grofile, n_cells, max_width, z_avg):
    surf = derive_surface(n_cells, small_grofile, max_width)
    for cell, calc_cell in zip(z_avg, surf):
        assert_almost_equal(cell, calc_cell)


@pytest.mark.parametrize('n_cells, max_width', [(3, 3)])
@pytest.mark.parametrize('dummy_array, z_avg', [(
    # dummy coordinates [x,y,z]
    np.array([[0., 0., 150.],
              [10., 0., 150.],
              [20., 0., 150.],
              [0., 10., 150.],
              [10., 10., 120.],
              [20., 10., 120.],
              [0., 20., 120.],
              [10., 20., 120.],
              [20., 20., 120.]]),
    # z_avg in grid
    np.array(([150., 150., 120.],
              [150., 120., 120.],
              [150., 120., 120.]))
)])
def test_derive_surface_from_numpy(dummy_array, n_cells, max_width, z_avg):
    u = mda.Universe(dummy_array, n_atoms=len(dummy_array))
    selection = u.select_atoms('index 0:9')
    surf = derive_surface(n_cells, selection, max_width)
    for cell, calc_cell in zip(z_avg, surf):
        assert_almost_equal(cell, calc_cell)


@pytest.mark.parametrize('n_cells, max_width', [(3, 3)])
@pytest.mark.parametrize('dummy_array, z_avg', [(
    # dummy coordinates [x,y,z]
    np.array([[0., 0., 150.],
              [0., 0., 150.],
              [20., 0., 150.],
              [0., 0., 150.],
              [10., 10., 150.],
              [20., 10., 150.],
              [0., 20., 150.],
              [10., 20., 150.],
              [20., 20., 150.]]),
    # z_avg in grid
    np.array(([150., np.nan, 150.],
              [np.nan, 150., 150.],
              [150., 150., 150.]))
)])
def test_derive_surface_from_numpy_with_nans(dummy_array, n_cells, max_width, z_avg):
    u = mda.Universe(dummy_array, n_atoms=len(dummy_array))
    selection = u.select_atoms('index 0:9')
    surf = derive_surface(n_cells, selection, max_width)
    for cell, calc_cell in zip(z_avg, surf):
        assert_almost_equal(cell, calc_cell)
