"""Unit tests for binning nearest-neighbor surface functions."""

import MDAnalysis as mda
import numpy as np
import pytest
from numpy.testing import assert_allclose

from membrane_curvature.base import MembraneCurvature
from membrane_curvature.binning_nearest_surface import get_z_surface_nearest, lipid_center_positions


@pytest.fixture
def box_100():
    return np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])


@pytest.fixture
def grid_2x2():
    return dict(n_x_bins=2, n_y_bins=2)


@pytest.fixture
def bilayer_universe(box_100):
    """Two single-bead lipids: upper at (90, 10, 80), lower at (10, 10, 20)."""
    universe = mda.Universe.empty(
        2,
        trajectory=True,
        n_residues=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 0],
    )
    universe.atoms.positions = np.array([[90.0, 10.0, 80.0], [10.0, 10.0, 20.0]])
    universe.dimensions = box_100
    universe.trajectory[0]
    return universe


@pytest.fixture
def xy_nearest_universe():
    """Lipid at (5, 5, 105) is xy-nearest to corner (0, 0), not (30, 30, 96)."""
    universe = mda.Universe.empty(
        2,
        trajectory=True,
        n_residues=2,
        atom_resindex=[0, 1],
        residue_segindex=[0, 0],
    )
    universe.atoms.positions = np.array([[5.0, 5.0, 105.0], [30.0, 30.0, 96.0]])
    universe.dimensions = np.array([100.0, 100.0, 233.0, 90.0, 90.0, 90.0])
    universe.trajectory[0]
    return universe


@pytest.fixture
def two_bead_residue_universe(box_100):
    """One residue with two atoms for mass-weighted center tests."""
    universe = mda.Universe.empty(
        2,
        trajectory=True,
        n_residues=1,
        atom_resindex=[0, 0],
        residue_segindex=[0],
    )
    universe.atoms.positions = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 10.0]])
    universe.add_TopologyAttr('masses')
    universe.atoms[0].mass = 1.0
    universe.atoms[1].mass = 3.0
    universe.dimensions = box_100
    return universe


@pytest.fixture
def nearest_surface(box_100, grid_2x2):
    def _call(atoms, *, box=None, grid_origin='box', x_range=None, y_range=None):
        return get_z_surface_nearest(
            atoms,
            box=box_100 if box is None else box,
            grid_origin=grid_origin,
            x_range=x_range,
            y_range=y_range,
            **grid_2x2,
        )

    return _call


def test_get_z_surface_nearest_upper_leaflet_mic(bilayer_universe, nearest_surface):
    z_surface, grid_axes = nearest_surface(bilayer_universe.atoms[[0]])
    assert_allclose(z_surface[0, 0], 80.0)
    assert_allclose(z_surface[1, 1], 80.0)
    assert_allclose(grid_axes['dx'], 50.0)
    assert_allclose(grid_axes['dy'], 50.0)


def test_get_z_surface_nearest_lower_leaflet(bilayer_universe, nearest_surface):
    z_surface, _grid_axes = nearest_surface(bilayer_universe.atoms[[1]])
    assert_allclose(z_surface[0, 0], 20.0)


def test_get_z_surface_nearest_uses_xy_distance_only(xy_nearest_universe, nearest_surface):
    z_surface, _grid_axes = nearest_surface(xy_nearest_universe.atoms)
    assert_allclose(z_surface[0, 0], 105.0)


def test_get_z_surface_nearest_requires_box(bilayer_universe, grid_2x2):
    with pytest.raises(ValueError, match='box is required'):
        get_z_surface_nearest(bilayer_universe.atoms[[0]], box=None, **grid_2x2)


def test_get_z_surface_nearest_range_override(bilayer_universe, nearest_surface):
    _z_surface, grid_axes = nearest_surface(
        bilayer_universe.atoms[[0]],
        x_range=(10.0, 90.0),
        y_range=(20.0, 80.0),
    )
    assert_allclose(grid_axes['dx'], 40.0)
    assert_allclose(grid_axes['dy'], 30.0)


def test_get_z_surface_nearest_lipid_bbox_extent(xy_nearest_universe, nearest_surface):
    _z_surface, grid_axes = nearest_surface(xy_nearest_universe.atoms, grid_origin='lipid_bbox')
    assert_allclose(grid_axes['left_x'], 5.0)
    assert_allclose(grid_axes['left_y'], 5.0)
    assert_allclose(grid_axes['dx'], 12.5)
    assert_allclose(grid_axes['dy'], 12.5)


@pytest.mark.parametrize(
    'masses, expected',
    [
        ([1.0, 3.0], [[7.5, 0.0, 7.5]]),  # mass-weighted center of residue
        ([0.0, 0.0], [[5.0, 0.0, 5.0]]),  # all-zero masses -> geometric center
    ],
)
def test_lipid_center_positions(two_bead_residue_universe, masses, expected):
    two_bead_residue_universe.atoms.masses = masses
    assert_allclose(lipid_center_positions(two_bead_residue_universe.atoms), expected)


def test_lipid_center_positions_without_masses(bilayer_universe):
    centers = lipid_center_positions(bilayer_universe.atoms[[0]])
    assert_allclose(centers, bilayer_universe.atoms[[0]].positions)


def test_membrane_curvature_binning_nearest_runs(bilayer_universe, grid_2x2):
    mc = MembraneCurvature(
        bilayer_universe,
        select='index 0',
        surface_method='binning_nearest',
        grid_origin='box',
        fft_filter=None,
        **grid_2x2,
    ).run()
    assert mc.results.z_surface.shape == (1, 2, 2)
    assert_allclose(mc.results.z_surface[0, 0, 0], 80.0)


@pytest.mark.parametrize(
    'kwargs, match',
    [
        (dict(wrap=True), r"wrap=True is not valid when surface_method='binning_nearest'"),
        (
            dict(select='index 0', grid_origin='lipid_bbox', fft_filter='auto'),
            r"grid_origin='lipid_bbox' cannot be combined with fft_filter",
        ),
        (
            dict(select='index 0', grid_origin='lipid_bbox', padding=True),
            r"grid_origin='lipid_bbox' cannot be combined with padding",
        ),
        (dict(grid_origin='invalid'), r"grid_origin must be 'box' or 'lipid_bbox'"),
    ],
)
def test_membrane_curvature_binning_nearest_raises(bilayer_universe, kwargs, match):
    with pytest.raises(ValueError, match=match):
        MembraneCurvature(bilayer_universe, surface_method='binning_nearest', **kwargs)


def test_membrane_curvature_binning_nearest_resolves_fft_filter(bilayer_universe, grid_2x2):
    mc = MembraneCurvature(
        bilayer_universe,
        select='index 0',
        surface_method='binning_nearest',
        grid_origin='box',
        fft_filter='auto',
        **grid_2x2,
    )
    assert mc._fft_q_bounds is not None
    assert mc._fft_q_bounds[0] == 0.0


def test_membrane_curvature_binning_nearest_with_padding(bilayer_universe, grid_2x2):
    mc = MembraneCurvature(
        bilayer_universe,
        select='index 0',
        surface_method='binning_nearest',
        grid_origin='box',
        fft_filter=None,
        padding=True,
        edge_pad_bins=1,
        **grid_2x2,
    ).run()
    assert mc.padding is True
    assert mc.results.z_surface.shape == (1, 2, 2)
    assert mc.results.mean.shape == (1, 2, 2)
    assert mc.results.gaussian.shape == (1, 2, 2)
    assert np.any(np.isfinite(mc.results.mean))


def test_membrane_curvature_binning_nearest_lipid_bbox_runs(xy_nearest_universe, grid_2x2):
    mc = MembraneCurvature(
        xy_nearest_universe,
        select='all',
        surface_method='binning_nearest',
        grid_origin='lipid_bbox',
        fft_filter=None,
        **grid_2x2,
    ).run()
    assert mc.results.z_surface.shape == (1, 2, 2)
    assert np.all(np.isfinite(mc.results.z_surface))
    assert np.all(np.isfinite(mc.results.mean))
