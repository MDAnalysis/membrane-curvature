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
from ..lib.mods import *

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal
import MDAnalysis as mda
from mdakit_membcurv.tests.datafiles import GRO_, XTC_, GRO_PO4

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
                                     -0.00642677, 0.0254065, -0.01830984, -0.00904487, -0.01182518]])

}


@pytest.fixture()
def universe():
    u = mda.Universe(GRO_, XTC_)
    return u


@pytest.fixture()
def ref_beads():
    u = mda.Universe(GRO_PO4)
    ref_beads = u.atoms[0:455]
    return ref_beads


@pytest.fixture()
def Z_cloud(ref_beads):
    n_cells, max_width = 10, 19
    grid_count, grid_1, grid_2, z_ref = [np.zeros([n_cells, n_cells]) for i in range(4)]
    factor = np.float32(n_cells / max_width)

    for atom in ref_beads:
        x, y, z = atom.position / 10

        l = int(abs(x) * factor)
        m = int(abs(y) * factor)

        try:
            grid_1[l, m] += z
            grid_2[l, m] += 1

        except BaseException:
            pass

    for i, j in it.product(range(n_cells), range(n_cells)):
        if grid_2[i, j] > 0:
            z_ref[i, j] += grid_1[i, j] / grid_2[i, j]
            grid_count[i, j] += 1

    for i, j in it.product(range(n_cells), range(n_cells)):
        if grid_count[i, j] > 0:
            z_ref[i, j] /= grid_count[i, j]

    return z_ref


@pytest.fixture()
def mdtraj():
    u = mda.Universe(GRO_, XTC_)
    return u


@pytest.fixture()
def leaflets():
    return ['lower', 'upper']


@pytest.fixture()
def lipid_types():
    return ['POPC', 'POPE']


def test_mdakit_membcurv_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mdakit_membcurv" in sys.modules


def test_output_folder():  # line 76, gone after refactoring
    tmp_dir = 'tmp_output/'
    os.makedirs(os.path.dirname(tmp_dir), exist_ok=True)
    assert os.path.exists(tmp_dir)
    os.rmdir(tmp_dir)


def test_md_topology():  # line 92, gone after refactoring
    top = md.load(GRO_).topology
    assert top.n_atoms == 18891


def test_box_size(universe):  # line 99
    box_size = universe.dimensions[0]
    assert_almost_equal(box_size, 184.31013, decimal=5)


def test_max_width(universe):  # line 100
    max_width = universe.dimensions[0] * 0.1
    assert max_width == 18.431013488769533


def test_n_cells(universe):  # line 101
    unit_width = 20
    max_width = universe.dimensions[0] * 0.1
    n_cells = math.ceil(max_width / unit_width * 10)
    print(n_cells)
    assert n_cells == 10.


def test_dict_to_pickle():  # line 117-118
    name = 'test_pickle_output'
    dict_ = {'A': 1, 'B': 2, 'C': 3}
    dict2pickle(name, dict_)
    unpickled = pickle.load(open(name + '.pickle', 'rb'))
    assert dict_ == unpickled
    os.remove(name + '.pickle')


def test_md_traj():  # line 108, gone after refactoring
    traj = md.load(XTC_, top=GRO_)
    assert traj.n_frames == 11


def test_gaussian_curvature(Z_cloud):  # line 391 of mods
    K_test = gaussian_curvature(Z_cloud)
    for k, k_test in zip(MEMBRANE_CURVATURE_DATA['gaussian_curvature'], K_test):
        assert_almost_equal(k, k_test)


def test_mean_curvature(Z_cloud):  # line 417 of mods
    H_test = mean_curvature(Z_cloud)
    for h, h_test in zip(MEMBRANE_CURVATURE_DATA['mean_curvature'], H_test):
        assert_almost_equal(h, h_test)
