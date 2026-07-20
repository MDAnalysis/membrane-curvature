"""
Unit tests for edge padding helpers and curvature on padded grids.
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
    _expand_inclusive_ranges,
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


@pytest.fixture
def padding_settings():
    return {
        'box': 100.0,
        'n_bins': 20,
        'amplitude': 10.0,
        'edge_pad_bins': 2,
    }


@pytest.fixture
def analytic(padding_settings):
    """Sinusoid height at bin centres, matching atom positions, and exact curvatures."""
    box = padding_settings['box']
    n_bins = padding_settings['n_bins']
    amp = padding_settings['amplitude']
    spacing = box / n_bins
    k = 2 * np.pi / box

    centres = (np.arange(n_bins) + 0.5) * spacing
    x, y = np.meshgrid(centres, centres, indexing='ij')
    sx, cx = np.sin(k * x), np.cos(k * x)
    sy, cy = np.sin(k * y), np.cos(k * y)

    z = amp * sx * cy
    fx = amp * k * cx * cy
    fy = -amp * k * sx * sy
    fxx = fyy = -amp * k**2 * sx * cy
    fxy = -amp * k**2 * cx * sy

    return {
        'z': z,
        'spacing': spacing,
        'positions': np.column_stack([x.ravel(), y.ravel(), z.ravel()]),
        'mean': mean_curvature_monge(fx, fy, fxx, fyy, fxy),
        'gaussian': gaussian_curvature_monge(fx, fy, fxx, fyy, fxy),
    }


def border_interior_ratio(field, reference, pad):
    error = np.abs(field - reference)
    border = np.ones(error.shape, dtype=bool)
    border[pad:-pad, pad:-pad] = False
    return error[border].max() / error[pad:-pad, pad:-pad].max()


def run_binning(positions, box, n_bins, **kwargs):
    universe = mda.Universe(positions, n_atoms=positions.shape[0])
    universe.dimensions = [box, box, box, 90.0, 90.0, 90.0]
    return MembraneCurvature(
        universe,
        select='all',
        surface_method='binning',
        fft_filter=None,
        n_x_bins=n_bins,
        n_y_bins=n_bins,
        **kwargs,
    ).run()


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

    @pytest.mark.parametrize('edge_pad_bins', [1, np.int64(3)])
    def test_accepts_valid_edge_pad_bins(self, edge_pad_bins):
        assert validate_edge_pad_bins(edge_pad_bins, n_x_bins=10, n_y_bins=10) == int(edge_pad_bins)

    @pytest.mark.parametrize('edge_pad_bins', [0, -1, 1.5, '2'])
    def test_rejects_invalid_edge_pad_bins(self, edge_pad_bins):
        with pytest.raises(ValueError, match=r'edge_pad_bins must be an int >= 1'):
            validate_edge_pad_bins(edge_pad_bins, n_x_bins=10, n_y_bins=10)

    @pytest.mark.parametrize('n_x_bins, n_y_bins', [(3, 10), (10, 3)])
    def test_rejects_edge_pad_bins_larger_than_bin_counts(self, n_x_bins, n_y_bins):
        with pytest.raises(ValueError, match=r'edge_pad_bins must be <= n_x_bins'):
            validate_edge_pad_bins(4, n_x_bins=n_x_bins, n_y_bins=n_y_bins)

    def test_warns_when_edge_pad_bins_larger_than_four(self):
        with pytest.warns(UserWarning, match=r'edge_pad_bins > 4'):
            assert validate_edge_pad_bins(5, n_x_bins=10, n_y_bins=10) == 5


class TestPaddingHelpers:
    def test_expand_inclusive_ranges_empty(self):
        item, values = _expand_inclusive_ranges(np.array([1]), np.array([0]))
        assert item.shape == (0,)
        assert values.shape == (0,)

    def test_padded_grid_spec(self):
        spec = padded_grid_spec(10, 20, (0.0, 100.0), (0.0, 200.0), 2)
        assert spec['n_x_bins_pad'] == 14
        assert spec['n_y_bins_pad'] == 24
        assert spec['crop'] == (2, 2)
        assert_allclose(spec['dx'], 10.0)
        assert_allclose(spec['dy'], 10.0)
        assert_allclose(spec['pad_x'], 20.0)
        assert_allclose(spec['pad_y'], 20.0)
        assert_allclose(spec['x_range_pad'], (-20.0, 120.0))
        assert_allclose(spec['y_range_pad'], (-20.0, 220.0))

    def test_tile_adds_edge_images(self):
        positions = np.array(
            [
                [5.0, 50.0, 1.0],
                [95.0, 50.0, 2.0],
                [50.0, 50.0, 3.0],
            ]
        )
        tiled = tile_xy_buffer(positions, 100.0, 100.0, pad_x=10.0, pad_y=10.0)
        assert tiled.shape[0] == 5
        assert np.any(np.isclose(tiled[:, 0], -5.0))
        assert np.any(np.isclose(tiled[:, 0], 105.0))

    def test_tile_unwrapped_matches_wrapped(self):
        wrapped = np.array([[5.0, 50.0, 1.0], [95.0, 50.0, 2.0]])
        unwrapped = wrapped.copy()
        unwrapped[1, 0] += 100.0

        a = tile_xy_buffer(wrapped, 100.0, 100.0, pad_x=10.0, pad_y=10.0)
        b = tile_xy_buffer(unwrapped, 100.0, 100.0, pad_x=10.0, pad_y=10.0)
        assert_allclose(a[np.lexsort(a.T[::-1])], b[np.lexsort(b.T[::-1])])

    def test_tile_empty_positions(self):
        tiled = tile_xy_buffer(np.empty((0, 3)), 100.0, 100.0, pad_x=10.0, pad_y=10.0)
        assert tiled.shape == (0, 3)

    def test_clip_padded_grid(self):
        padded = np.arange(16, dtype=float).reshape(4, 4)
        assert_allclose(clip_padded_grid(padded, 1), [[5.0, 6.0], [9.0, 10.0]])

    def test_z_surface_empty_coordinates(self):
        z_padded = get_z_surface_padded(
            np.empty((0, 3)),
            n_x_bins=2,
            n_y_bins=2,
            x_range=(0.0, 100.0),
            y_range=(0.0, 100.0),
            edge_pad_bins=1,
        )
        assert z_padded.shape == (4, 4)
        assert np.all(np.isnan(z_padded))

    def test_z_surface_bins_mean_z(self):
        coordinates = np.array(
            [
                [5.0, 5.0, 10.0],
                [5.0, 5.0, 20.0],
                [15.0, 15.0, 5.0],
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

    def test_tiling_fills_buffer_where_wrap_pad_is_nan(self):
        # Atom outside [0, L); image at x=10 lands in the padded domain.
        coordinates = np.array([[110.0, 50.0, 8.0]])
        grid_kwargs = dict(
            n_x_bins=4,
            n_y_bins=4,
            x_range=(0.0, 100.0),
            y_range=(0.0, 100.0),
        )
        z_tiled = get_z_surface_padded(
            coordinates,
            edge_pad_bins=1,
            box_length_x=100.0,
            box_length_y=100.0,
            **grid_kwargs,
        )
        with pytest.warns(UserWarning, match=r'outside'):
            z_primary = get_z_surface(coordinates, **grid_kwargs)

        assert_allclose(z_tiled[1, 3], 8.0)
        assert np.isnan(np.pad(z_primary, 1, mode='wrap')).all()

    def test_flat_surface_zero_curvature(self):
        z, mean, gaussian = curvature_with_edge_pad(
            np.full((4, 4), 5.0),
            dx=1.0,
            dy=1.0,
            edge_pad_bins=1,
        )
        assert_allclose(z, 5.0)
        assert_allclose(mean, 0.0, atol=1e-12)
        assert_allclose(gaussian, 0.0, atol=1e-12)

        z, mean, gaussian = curvature_from_primary_with_edge_pad(
            np.full((3, 3), 2.0),
            dx=1.0,
            dy=1.0,
            edge_pad_bins=1,
        )
        assert_allclose(z, 2.0)
        assert_allclose(mean, 0.0, atol=1e-12)
        assert_allclose(gaussian, 0.0, atol=1e-12)

    @pytest.mark.parametrize('edge_pad_bins', [0, -1])
    def test_rejects_non_positive_pad(self, edge_pad_bins):
        with pytest.raises(ValueError, match=r'edge_pad_bins must be >= 1'):
            curvature_from_primary_with_edge_pad(
                np.zeros((3, 3)),
                dx=1.0,
                dy=1.0,
                edge_pad_bins=edge_pad_bins,
            )

    @pytest.mark.parametrize('key', ['mean', 'gaussian'])
    def test_wrap_pad_equalises_border_error(self, key, analytic, padding_settings):
        pad = padding_settings['edge_pad_bins']
        _, mean, gaussian = curvature_from_primary_with_edge_pad(
            analytic['z'], analytic['spacing'], analytic['spacing'], pad
        )
        field = mean if key == 'mean' else gaussian
        ratio = border_interior_ratio(field, analytic[key], pad)
        assert ratio == pytest.approx(1.0, abs=0.1)


class TestPaddingAnalysis:
    @pytest.mark.parametrize('attr', ['mean', 'gaussian'])
    def test_padding_removes_edge_artifact(self, attr, analytic, padding_settings):
        box = padding_settings['box']
        n_bins = padding_settings['n_bins']
        pad = padding_settings['edge_pad_bins']
        positions = analytic['positions']

        unpadded = getattr(run_binning(positions, box, n_bins, padding=False).results, attr)[0]
        padded = getattr(
            run_binning(positions, box, n_bins, padding=True, edge_pad_bins=pad).results,
            attr,
        )[0]

        assert border_interior_ratio(unpadded, analytic[attr], pad) > 4.0
        assert border_interior_ratio(padded, analytic[attr], pad) == pytest.approx(1.0, abs=0.1)

    @pytest.mark.parametrize('attr', ['z_surface', 'mean', 'gaussian'])
    def test_padding_recovers_atoms_outside_box(self, attr, analytic, padding_settings):
        box = padding_settings['box']
        n_bins = padding_settings['n_bins']
        pad = padding_settings['edge_pad_bins']
        positions = analytic['positions']

        drifted = positions.copy()
        drifted[drifted[:, 0] > 0.8 * box, 0] += box

        reference = run_binning(positions, box, n_bins, wrap=False, padding=True, edge_pad_bins=pad)
        displaced = run_binning(drifted, box, n_bins, wrap=False, padding=True, edge_pad_bins=pad)
        assert_allclose(
            getattr(displaced.results, attr)[0],
            getattr(reference.results, attr)[0],
        )

    def test_atoms_outside_box_dropped_without_padding(self, analytic, padding_settings):
        box = padding_settings['box']
        n_bins = padding_settings['n_bins']
        drifted = analytic['positions'].copy()
        drifted[drifted[:, 0] > 0.8 * box, 0] += box

        with pytest.warns(UserWarning):
            displaced = run_binning(drifted, box, n_bins, wrap=False, padding=False)
        assert np.isnan(displaced.results.z_surface[0]).any()
