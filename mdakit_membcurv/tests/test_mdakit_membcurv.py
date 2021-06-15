"""
Unit and regression test for the mdakit_membcurv package.
"""

# Import package, test suite, and other packages as needed
import pickle

import pytest
import sys
import math
import os
import mdtraj as md  # This will be gone after refactoring
import itertools as it
from ..lib.mods import core_fast_leaflet, curvature, mean_curvature, gaussian_curvature, dict2pickle
import numpy as np
from numpy.testing import assert_almost_equal
import MDAnalysis as mda
from mdakit_membcurv.tests.datafiles import GRO_MEMBRANE_PROTEIN, XTC_MEMBRANE_PROTEIN, GRO_PO4, XTC_PO4

# Reference data from datafile
MEMBRANE_CURVATURE_DATA = {
    'mean_curvature': np.array([[0.06396822, 0.2254189, 0.22399141, 0.02498465, 0.05090628,
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
                                 -0.09468121, -0.21140388, 0.03506431, 0.03308529, -0.01943328]]),

    'gaussian_curvature': np.array([[-0.00189129, 0.04934495, 0.04660196, -0.00351411, -0.01943406,
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
                                     -0.00642677, 0.0254065, -0.01830984, -0.00904487, -0.01182518]]),

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

    'mean_from_z_avg': {'upper': np.array([[0.12555075, 0.02135105, 0.04160992, -0.01648831, 0.00093022,
                                            0.05210213, -0.01413144, 0.42855804, 0.11933335, -0.15700576],
                                           [0.0119322, -0.03312442, 0.16867713, 0.11186921, -0.06924967,
                                            -0.17337403, -0.03075403, 0.00998621, 0.03106606, 0.12265291],
                                           [0.01316543, 0.00331553, 0.021467, 0.14910619, -0.02331827,
                                            -0.40607687, 0.00254864, 0.09637312, 0.00504627, 0.45395702],
                                           [0.07435368, 0.01111333, -0.26369175, 0.02026796, 0.15654274,
                                            -0.04698233, -0.07173794, -0.00261281, -0.11980999, -0.1344619],
                                           [-0.04629355, -0.00581553, 0.03905765, 0.20953977, 0.3918796,
                                            0.0310723, -0.42269884, -0.16012767, -0.01203511, -0.04650698],
                                           [-0.00449059, 0.01816231, -0.03467836, -0.09131014, 0.23924521,
                                            0.3350186, -0.07902011, -0.16677442, -0.02237877, -0.01612266],
                                           [-0.03812502, -0.00993137, -0.05389885, -0.23281703, -0.03940271,
                                            0.22679189, 0.15163259, 0.00348903, -0.03725647, 0.23348912],
                                           [-0.10594603, -0.14290049, -0.04246863, -0.05403972, -0.14319528,
                                            -0.02398386, 0.06424594, 0.03606966, 0.03191477, -0.03725069],
                                           [-0.00526536, -0.04075257, -0.04120521, -0.03985792, -0.03383821,
                                            -0.03392527, -0.04174206, 0.31247241, 0.19253582, -0.27770588],
                                           [0.04900292, 0.03659575, -0.03416543, -0.0529405, 0.00695103,
                                            0.01313807, -0.03327879, 0.01243288, 0.04821536, -0.09293237]])},

    'gaussian_from_z_avg': {'upper': np.array([[-2.39194800e-02, -3.61634194e-03, -1.08310985e-02,
                                                -4.04262352e-02, -1.07953930e-01, -1.58940765e-03,
                                                -5.14815931e-02, -4.08411708e-03, -1.73969187e-01,
                                                -1.27348350e-02],
                                               [-1.80881429e-02, -8.37058009e-04, 1.55266256e-02,
                                                7.96503053e-03, -8.70989095e-02, 2.62876632e-02,
                                                -5.03604863e-02, -7.13194031e-04, -8.88446427e-01,
                                                -1.94381867e-01],
                                               [-6.40974462e-03, -9.66907198e-03, 1.30718905e-04,
                                                -5.47821599e-02, -1.76480851e-02, 1.54886560e-01,
                                                -1.36958801e-02, 2.75305642e-03, -3.89289835e-03,
                                                -2.05699114e-02],
                                               [-2.13816311e-02, -3.60467810e-03, 6.62439558e-02,
                                                -9.62066579e-03, -2.83642243e-02, -6.27955986e-02,
                                                -5.31456175e-02, -8.46882808e-03, 1.37664910e-02,
                                                1.68963909e-02],
                                               [-2.57613213e-03, -5.75152336e-03, -3.04009492e-03,
                                                4.38720368e-02, 8.95669062e-02, -2.11138887e-03,
                                                1.21192352e-01, 2.35164083e-02, -1.32594642e-02,
                                                -8.60185561e-03],
                                               [-6.55374523e-03, -7.43313172e-03, -4.29013834e-02,
                                                -1.93519324e-03, 2.36387158e-02, 2.59636175e-02,
                                                4.16371248e-03, -5.25782227e-03, -2.71611977e-03,
                                                3.10811328e-05],
                                               [-5.32693655e-03, -2.85760622e-03, -5.49568585e-04,
                                                2.04948874e-02, -1.76805356e-02, -1.36782305e-02,
                                                -3.16246302e-02, -5.62996471e-03, -1.44329381e-04,
                                                -3.73865110e-02],
                                               [6.22342127e-03, 1.94289676e-02, 2.82196188e-05,
                                                2.69411991e-03, -1.30866347e-02, -2.52293095e-02,
                                                -7.22431515e-03, -6.73201764e-03, -2.51286123e-01,
                                                -3.22793339e-01],
                                               [-2.24613530e-02, 7.82002286e-04, -3.98911117e-04,
                                                5.29635968e-04, -4.19699155e-03, -1.30479296e-02,
                                                -2.38928651e-04, 4.24992903e-02, 1.70628313e-02,
                                                -4.85138439e-02],
                                               [-1.94220106e-02, -2.13212988e-04, 1.05435416e-03,
                                                1.30710748e-03, -3.49076390e-03, -1.28226749e-02,
                                                -1.41168166e-02, -1.08870941e-02, -1.21323346e+00,
                                                -1.55778010e-01]])}
}


@pytest.fixture()
def universe():
    u = mda.Universe(GRO_MEMBRANE_PROTEIN, XTC_MEMBRANE_PROTEIN)
    return u


@pytest.fixture()
def ref_beads():
    u = mda.Universe(GRO_PO4)
    ref_beads = u.atoms[0:455]
    return ref_beads


@pytest.fixture()
def Z_cloud(ref_beads):
    n_cells, max_width = 10, 19

    # set grid_1. z coordinates for [i,j] cells.
    # set grid_2. Count number of beads populating each [i,j] cell.
    # set z_ref. Final values of average z per [i,j] cell.
    grid_1, grid_2, z_ref = [np.zeros([n_cells, n_cells]) for i in range(3)]

    # factor used to map (x,y) to [i,j]
    factor = np.float32(n_cells / max_width)

    # iterate over the bead of reference
    for atom in ref_beads:

        # extract positions in nm.
        x, y, z = atom.position / 10

        # map (x,y) to [i,j]
        cell_l = int(abs(x) * factor)
        cell_m = int(abs(y) * factor)

        # sum value of z coordinate
        grid_1[cell_l, cell_m] += z
        grid_2[cell_l, cell_m] += 1

    # To calculate average, iterate over each cell in the grid:
    for i, j in it.product(range(n_cells), range(n_cells)):
        # if the element [i,j] is not empty
        if grid_2[i, j] > 0:
            # then calculate the average of z:
            # grid_1 has the sum of all the z coordinates.
            # grid_2 counted how many beads were in that grid.
            z_ref[i, j] += grid_1[i, j] / grid_2[i, j]
        else:
            # if there are not beads in that cell, store a nan.
            z_ref[i, j] = np.nan

    return z_ref


@pytest.fixture()
def mdtraj_po4():
    # trajectory PO4 beads only imported using mdtraj, gone after refactoring
    mdtraj = md.load(XTC_PO4, top=GRO_PO4)
    return mdtraj


@pytest.fixture(scope="session")
def output(tmpdir_factory):
    tmp = tmpdir_factory.mktemp("output")
    return tmp


@pytest.fixture()
def md_ref_beads():
    # reference beads using mdtraj, gone after refactoring. Select upper leaflet
    topology = md.load(GRO_PO4).topology
    md_ref_beads = {'upper': {'POPC': topology.select('resname POPC and index 1 to 457').astype(int).tolist()}}
    return md_ref_beads


def test_mdakit_membcurv_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mdakit_membcurv" in sys.modules


def test_md_topology():  # line 92, gone after refactoring
    top = md.load(GRO_MEMBRANE_PROTEIN).topology
    assert top.n_atoms == 18891


def test_box_size(universe):  # line 99
    box_size = universe.dimensions[0]
    assert_almost_equal(box_size, 184.31013, decimal=5)


def test_n_cells(universe):  # line 101
    unit_width = 20
    max_width = universe.dimensions[0] * 0.1
    n_cells = math.ceil(max_width / unit_width * 10)
    print(n_cells)
    assert n_cells == 10.


def test_dict_to_pickle(output):  # line 117-118
    name = 'test_pickle_output'
    dict_ = {'A': 1, 'B': 2, 'C': 3}
    with output.as_cwd():
        dict2pickle(name, dict_)
        unpickled = pickle.load(open(name + '.pickle', 'rb'))
        assert dict_ == unpickled


def test_gaussian_curvature(Z_cloud):  # line 391 of mods
    K_test = gaussian_curvature(Z_cloud)
    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_curvature'], K_test):
        assert_almost_equal(k, k_test)


def test_mean_curvature(Z_cloud):  # line 417 of mods
    H_test = mean_curvature(Z_cloud)
    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_curvature'], H_test):
        assert_almost_equal(h, h_test)


def test_core_fast_leaflet(md_ref_beads, mdtraj_po4):
    jump = 1
    n_cells = 10
    max_width = 19
    z_calc = np.zeros([n_cells, n_cells])
    core_fast_leaflet(z_calc, "upper", mdtraj_po4, jump, n_cells, ["POPC"], md_ref_beads, max_width)
    for z, z_test in zip(MEMBRANE_CURVATURE_DATA['z_avg_coords'], z_calc):
        assert_almost_equal(z, z_test)


def test_curvature(md_ref_beads, mdtraj_po4):
    jump = 1
    n_cells = 10
    max_width = 19
    z_test = np.zeros([n_cells, n_cells])
    lf = 'upper'
    core_fast_leaflet(z_test, lf, mdtraj_po4, jump, n_cells, ["POPC"], md_ref_beads, max_width)
    z_ref = {lf: z_test}
    H_test, K_test = curvature(z_ref, [lf], n_cells)

    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_from_z_avg'][lf], H_test[lf]):
        assert_almost_equal(h, h_test)

    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_from_z_avg'][lf], K_test[lf]):
        assert_almost_equal(k, k_test)
