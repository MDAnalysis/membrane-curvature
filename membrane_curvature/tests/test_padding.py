"""
Unit tests for edge padding: validators, helpers, and curvature on padded grids.
"""

import MDAnalysis as mda
import numpy as np
import pytest
from numpy.testing import assert_allclose

from membrane_curvature.base import MembraneCurvature
from membrane_curvature.binning_surface import get_z_surface
from membrane_curvature.curvature import (
    curvature_from_primary_with_edge_pad,
    curvature_with_edge_pad,
    gaussian_curvature_monge,
    mean_curvature_monge,
)
from membrane_curvature.padding import (
    clip_padded_grid,
    get_z_surface_padded,
    get_z_surface_with_edge_pad,
    padded_grid_spec,
    tile_xy_buffer,
)
from membrane_curvature.padding_validators import (
    validate_edge_pad_bins,
    validate_orthorhombic,
)

BOX_XY = 100.0
X_RANGE = (0.0, BOX_XY)
Y_RANGE = (0.0, BOX_XY)
SURFACE_BINS = 20
SURFACE_AMPLITUDE = 10.0
# Second derivatives reach this many cells in from each edge.
ARTIFACT_DEPTH = 2


@pytest.fixture()
def analytic_surface():
    """
    Height at bin centres, with reference curvature from exact derivatives.
    """
    spacing = BOX_XY / SURFACE_BINS
    k = 2 * np.pi / BOX_XY
    centres = (np.arange(SURFACE_BINS) + 0.5) * spacing
    x, y = np.meshgrid(centres, centres, indexing='ij')

    sin_x, cos_x = np.sin(k * x), np.cos(k * x)
    sin_y, cos_y = np.sin(k * y), np.cos(k * y)

    z = SURFACE_AMPLITUDE * sin_x * cos_y
    fx = SURFACE_AMPLITUDE * k * cos_x * cos_y
    fy = -SURFACE_AMPLITUDE * k * sin_x * sin_y
    fxx = -SURFACE_AMPLITUDE * k**2 * sin_x * cos_y
    fyy = -SURFACE_AMPLITUDE * k**2 * sin_x * cos_y
    fxy = -SURFACE_AMPLITUDE * k**2 * cos_x * sin_y

    reference = {
        'mean': mean_curvature_monge(fx, fy, fxx, fyy, fxy),
        'gaussian': gaussian_curvature_monge(fx, fy, fxx, fyy, fxy),
    }
    return z, spacing, reference


@pytest.fixture()
def analytic_positions(analytic_surface):
    """
    Atoms sit at bin centres to reproduce the binned height field exactly.
    """
    z, spacing, _ = analytic_surface
    centres = (np.arange(SURFACE_BINS) + 0.5) * spacing
    x, y = np.meshgrid(centres, centres, indexing='ij')
    return np.column_stack([x.ravel(), y.ravel(), z.ravel()])


@pytest.fixture()
def edge_atom_positions():
    """Atoms near left/right edges and one in the box centre."""
    return np.array(
        [
            [5.0, 50.0, 1.0],
            [95.0, 50.0, 2.0],
            [50.0, 50.0, 3.0],
        ]
    )


@pytest.fixture()
def pad_kwargs():
    return dict(
        n_x_bins=2,
        n_y_bins=2,
        x_range=X_RANGE,
        y_range=Y_RANGE,
        edge_pad_bins=1,
    )


@pytest.fixture()
def make_universe():
    """Factory for a single-frame orthorhombic universe from positions."""

    def factory(positions):
        universe = mda.Universe(positions, n_atoms=positions.shape[0])
        universe.dimensions = [BOX_XY, BOX_XY, BOX_XY, 90.0, 90.0, 90.0]
        return universe

    return factory


@pytest.fixture()
def run_binning(make_universe):
    """Factory that runs MembraneCurvature binning on given positions."""

    def factory(positions, **kwargs):
        params = dict(
            select='all',
            surface_method='binning',
            fft_filter=None,
            n_x_bins=SURFACE_BINS,
            n_y_bins=SURFACE_BINS,
        )
        params.update(kwargs)
        return MembraneCurvature(make_universe(positions), **params).run()

    return factory


@pytest.fixture()
def edge_interior_error_ratio():
    """
    Factory for worst border error divided by worst interior error.

    A correctly padded field has a ratio near 1. An unpadded field is much
    larger because only the border carries one-sided difference error.
    """

    def factory(field, reference):
        error = np.abs(field - reference)
        border = np.ones(error.shape, dtype=bool)
        border[ARTIFACT_DEPTH:-ARTIFACT_DEPTH, ARTIFACT_DEPTH:-ARTIFACT_DEPTH] = False
        interior = error[ARTIFACT_DEPTH:-ARTIFACT_DEPTH, ARTIFACT_DEPTH:-ARTIFACT_DEPTH]
        return error[border].max() / interior.max()

    return factory


class TestPaddingValidators:
    def test_accepts_orthorhombic(self):
        validate_orthorhombic([10.0, 20.0, 30.0, 90.0, 90.0, 90.0])

    @pytest.mark.parametrize(
        'dimensions, match',
        [
            ([10.0, 20.0, 30.0, 90.0, 90.0, 60.0], r'orthorhombic'),
            ([10.0, 20.0, 30.0], r'dimensions must include'),
        ],
    )
    def test_rejects_invalid_dimensions(self, dimensions, match):
        with pytest.raises(ValueError, match=match):
            validate_orthorhombic(dimensions)

    @pytest.mark.parametrize('edge_pad_bins', [1, 2, np.int64(3)])
    def test_accepts_valid_edge_pad_bins(self, edge_pad_bins):
        assert validate_edge_pad_bins(edge_pad_bins, n_x_bins=10, n_y_bins=10) == int(edge_pad_bins)

    @pytest.mark.parametrize('edge_pad_bins', [0, -1, 1.5, '2'])
    def test_rejects_invalid_edge_pad_bins(self, edge_pad_bins):
        with pytest.raises(ValueError, match=r'edge_pad_bins must be an int >= 1'):
            validate_edge_pad_bins(edge_pad_bins, n_x_bins=10, n_y_bins=10)

    @pytest.mark.parametrize(
        'edge_pad_bins, n_x_bins, n_y_bins',
        [
            (5, 3, 10),
            (5, 10, 3),
        ],
    )
    def test_rejects_edge_pad_bins_larger_than_bin_counts(self, edge_pad_bins, n_x_bins, n_y_bins):
        with pytest.raises(ValueError, match=r'edge_pad_bins must be <= n_x_bins'):
            validate_edge_pad_bins(edge_pad_bins, n_x_bins=n_x_bins, n_y_bins=n_y_bins)

    def test_warns_when_edge_pad_bins_larger_than_four(self):
        with pytest.warns(UserWarning, match=r'edge_pad_bins > 4'):
            assert validate_edge_pad_bins(5, n_x_bins=10, n_y_bins=10) == 5


class TestPaddedGridSpec:
    @pytest.mark.parametrize(
        'n_x_bins, n_y_bins, x_range, y_range, edge_pad_bins, expected',
        [
            (
                10,
                20,
                (0.0, 100.0),
                (0.0, 200.0),
                2,
                {
                    'n_x_bins_pad': 14,
                    'n_y_bins_pad': 24,
                    'dx': 10.0,
                    'dy': 10.0,
                    'pad_x': 20.0,
                    'pad_y': 20.0,
                    'x_range_pad': (-20.0, 120.0),
                    'y_range_pad': (-20.0, 220.0),
                    'crop': (2, 2),
                },
            ),
            (
                2,
                2,
                X_RANGE,
                Y_RANGE,
                1,
                {
                    'n_x_bins_pad': 4,
                    'n_y_bins_pad': 4,
                    'dx': 50.0,
                    'dy': 50.0,
                    'pad_x': 50.0,
                    'pad_y': 50.0,
                    'x_range_pad': (-50.0, 150.0),
                    'y_range_pad': (-50.0, 150.0),
                    'crop': (1, 1),
                },
            ),
        ],
    )
    def test_padded_grid_spec(self, n_x_bins, n_y_bins, x_range, y_range, edge_pad_bins, expected):
        spec = padded_grid_spec(n_x_bins, n_y_bins, x_range, y_range, edge_pad_bins)
        assert spec['n_x_bins_pad'] == expected['n_x_bins_pad']
        assert spec['n_y_bins_pad'] == expected['n_y_bins_pad']
        assert spec['crop'] == expected['crop']
        assert_allclose(spec['dx'], expected['dx'])
        assert_allclose(spec['dy'], expected['dy'])
        assert_allclose(spec['pad_x'], expected['pad_x'])
        assert_allclose(spec['pad_y'], expected['pad_y'])
        assert_allclose(spec['x_range_pad'], expected['x_range_pad'])
        assert_allclose(spec['y_range_pad'], expected['y_range_pad'])


class TestTileXyBuffer:
    def test_adds_edge_images(self, edge_atom_positions):
        tiled = tile_xy_buffer(
            edge_atom_positions,
            box_length_x=BOX_XY,
            box_length_y=BOX_XY,
            pad_x=10.0,
            pad_y=10.0,
        )
        # primary + left image of right-edge atom + right image of left-edge atom
        assert tiled.shape[0] == 5
        assert np.any(np.isclose(tiled[:, 0], -5.0))
        assert np.any(np.isclose(tiled[:, 0], 105.0))

    def test_unwrapped_matches_wrapped(self):
        wrapped = np.array(
            [
                [5.0, 50.0, 1.0],
                [95.0, 50.0, 2.0],
            ]
        )
        unwrapped = wrapped.copy()
        unwrapped[1, 0] += BOX_XY

        a = tile_xy_buffer(wrapped, BOX_XY, BOX_XY, pad_x=10.0, pad_y=10.0)
        b = tile_xy_buffer(unwrapped, BOX_XY, BOX_XY, pad_x=10.0, pad_y=10.0)
        order_a = np.lexsort(a.T[::-1])
        order_b = np.lexsort(b.T[::-1])
        assert_allclose(a[order_a], b[order_b])

    def test_empty_positions(self):
        tiled = tile_xy_buffer(np.empty((0, 3)), BOX_XY, BOX_XY, pad_x=10.0, pad_y=10.0)
        assert tiled.shape == (0, 3)

    @pytest.mark.parametrize(
        'pad_x, pad_y',
        [
            (-50.0, 10.0),
            (10.0, -50.0),
        ],
    )
    def test_empty_when_domain_collapsed(self, pad_x, pad_y):
        positions = np.array([[50.0, 50.0, 1.0]])
        tiled = tile_xy_buffer(positions, BOX_XY, BOX_XY, pad_x=pad_x, pad_y=pad_y)
        assert tiled.shape == (0, 3)


class TestClipPaddedGrid:
    def test_clips_to_primary(self):
        padded = np.arange(16, dtype=float).reshape(4, 4)
        clipped = clip_padded_grid(padded, edge_pad_bins=1)
        assert_allclose(clipped, np.array([[5.0, 6.0], [9.0, 10.0]]))


class TestGetZSurfacePadded:
    def test_empty_coordinates(self, pad_kwargs):
        z_padded = get_z_surface_padded(np.empty((0, 3)), **pad_kwargs)
        assert z_padded.shape == (4, 4)
        assert np.all(np.isnan(z_padded))

    def test_bins_mean_z(self):
        coordinates = np.array(
            [
                [5.0, 5.0, 10.0],
                [5.0, 5.0, 20.0],  # same primary bin (0, 0) and mean z = 15
                [15.0, 15.0, 5.0],  # primary bin (1, 1)
            ]
        )
        z_surface = get_z_surface_with_edge_pad(
            coordinates,
            n_x_bins=2,
            n_y_bins=2,
            x_range=(0.0, 20.0),
            y_range=(0.0, 20.0),
            edge_pad_bins=1,
        )
        assert_allclose(z_surface[0, 0], 15.0)
        assert_allclose(z_surface[1, 1], 5.0)
        assert np.isnan(z_surface[0, 1])
        assert np.isnan(z_surface[1, 0])

    @pytest.mark.parametrize(
        'box_length_x, box_length_y',
        [
            (None, None),
            (BOX_XY, BOX_XY),
        ],
    )
    def test_places_atom_on_padded_grid(self, box_length_x, box_length_y, pad_kwargs):
        coordinates = np.array([[25.0, 25.0, 7.0]])
        z_padded = get_z_surface_padded(
            coordinates,
            box_length_x=box_length_x,
            box_length_y=box_length_y,
            **pad_kwargs,
        )
        assert z_padded.shape == (4, 4)
        # Padded x domain is [-50, 150). Atom at x=25 maps to bin floor((25+50)/50) = 1.
        assert_allclose(z_padded[1, 1], 7.0)

    def test_atom_tiling_fills_buffer_where_wrap_pad_copies_nan(self):
        """
        Atom tiling and ``np.pad(..., mode='wrap')`` are not interchangeable.

        An atom just outside the primary cell still contributes to the padded
        domain via its lattice image. Binning the primary grid alone and then
        wrap-padding copies empty (NaN) cells into the buffer instead.
        """
        n_bins = 4
        edge_pad_bins = 1
        # Outside [0, L) along x axis. Image at x=10 still lands in the padded domain.
        coordinates = np.array([[110.0, 50.0, 8.0]])
        grid_kwargs = dict(
            n_x_bins=n_bins,
            n_y_bins=n_bins,
            x_range=(0.0, BOX_XY),
            y_range=(0.0, BOX_XY),
        )

        z_tiled_padded = get_z_surface_padded(
            coordinates,
            edge_pad_bins=edge_pad_bins,
            box_length_x=BOX_XY,
            box_length_y=BOX_XY,
            **grid_kwargs,
        )
        with pytest.warns(UserWarning, match=r'exceed'):
            z_primary = get_z_surface(coordinates, **grid_kwargs)
        z_wrap_padded = np.pad(z_primary, edge_pad_bins, mode='wrap')

        assert np.isnan(z_primary).all()
        # Image at x=10 is in the padded bin floor((10+25)/25) = 1
        # and y=50 is in bin 3.
        assert_allclose(z_tiled_padded[1, 3], 8.0)
        assert np.isnan(z_wrap_padded).all()


class TestCurvatureEdgePad:
    def test_curvature_with_edge_pad_flat_surface(self):
        z_padded = np.full((4, 4), 5.0)
        z_surface, mean, gaussian = curvature_with_edge_pad(
            z_padded,
            dx=1.0,
            dy=1.0,
            edge_pad_bins=1,
        )
        assert z_surface.shape == (2, 2)
        assert mean.shape == (2, 2)
        assert gaussian.shape == (2, 2)
        assert_allclose(z_surface, 5.0)
        assert_allclose(mean, 0.0, atol=1e-12)
        assert_allclose(gaussian, 0.0, atol=1e-12)

    def test_curvature_from_primary_with_edge_pad_flat_surface(self):
        z_surface = np.full((3, 3), 2.0)
        z_out, mean, gaussian = curvature_from_primary_with_edge_pad(
            z_surface,
            dx=1.0,
            dy=1.0,
            edge_pad_bins=1,
        )
        assert_allclose(z_out, z_surface)
        assert_allclose(mean, 0.0, atol=1e-12)
        assert_allclose(gaussian, 0.0, atol=1e-12)

    @pytest.mark.parametrize('edge_pad_bins', [0, -1])
    def test_curvature_from_primary_rejects_non_positive_pad(self, edge_pad_bins):
        with pytest.raises(ValueError, match=r'edge_pad_bins must be >= 1'):
            curvature_from_primary_with_edge_pad(
                np.zeros((3, 3)),
                dx=1.0,
                dy=1.0,
                edge_pad_bins=edge_pad_bins,
            )

    @pytest.mark.parametrize('key', ['mean', 'gaussian'])
    def test_wrap_pad_removes_edge_artifact_on_analytic_surface(self, analytic_surface, key, edge_interior_error_ratio):
        """Wrap-pad of a filled primary height field equalises border and interior error."""
        z, spacing, reference = analytic_surface
        _, mean, gaussian = curvature_from_primary_with_edge_pad(z, spacing, spacing, edge_pad_bins=ARTIFACT_DEPTH)
        field = mean if key == 'mean' else gaussian
        ratio = edge_interior_error_ratio(field, reference[key])
        assert ratio == pytest.approx(1.0, abs=0.1)


class TestPaddingAnalysis:
    """End to end padding behaviour through :class:`MembraneCurvature`."""

    def test_binned_surface_reproduces_analytic_height(self, analytic_surface, analytic_positions, run_binning):
        """Atoms at bin centres bin back to the analytic surface exactly."""
        z, _, _ = analytic_surface
        analysis = run_binning(analytic_positions, padding=True, edge_pad_bins=ARTIFACT_DEPTH)
        assert_allclose(analysis.results.z_surface[0], z)

    @pytest.mark.parametrize('attribute, key', [('mean', 'mean'), ('gaussian', 'gaussian')])
    def test_padding_removes_edge_artifact(
        self, analytic_surface, analytic_positions, attribute, key, run_binning, edge_interior_error_ratio
    ):
        _, _, reference = analytic_surface

        unpadded = getattr(run_binning(analytic_positions, padding=False).results, attribute)[0]
        padded = getattr(
            run_binning(analytic_positions, padding=True, edge_pad_bins=ARTIFACT_DEPTH).results, attribute
        )[0]

        assert edge_interior_error_ratio(unpadded, reference[key]) > 4.0
        assert edge_interior_error_ratio(padded, reference[key]) == pytest.approx(1.0, abs=0.1)

    @pytest.mark.parametrize('attribute', ['z_surface', 'mean', 'gaussian'])
    def test_padding_recovers_atoms_outside_the_box(self, analytic_positions, attribute, run_binning):
        """
        Tiling supplies lateral PBC, so displacing atoms by a box length is a no-op.

        This is what padding does that a wrap-pad of the primary height field
        cannot: the buffer is filled from image coordinates, so atoms outside
        the primary cell still land in a bin instead of being dropped.
        """
        drifted = analytic_positions.copy()
        outside = drifted[:, 0] > 0.8 * BOX_XY
        drifted[outside, 0] += BOX_XY
        assert outside.sum() > 0

        reference = run_binning(analytic_positions, wrap=False, padding=True, edge_pad_bins=ARTIFACT_DEPTH)
        displaced = run_binning(drifted, wrap=False, padding=True, edge_pad_bins=ARTIFACT_DEPTH)

        assert_allclose(
            getattr(displaced.results, attribute)[0],
            getattr(reference.results, attribute)[0],
        )

    def test_atoms_outside_the_box_are_dropped_without_padding(self, analytic_positions, run_binning):
        """Contrast with :meth:`test_padding_recovers_atoms_outside_the_box`."""
        drifted = analytic_positions.copy()
        drifted[drifted[:, 0] > 0.8 * BOX_XY, 0] += BOX_XY

        with pytest.warns(UserWarning):
            displaced = run_binning(drifted, wrap=False, padding=False)

        assert np.isnan(displaced.results.z_surface[0]).any()
