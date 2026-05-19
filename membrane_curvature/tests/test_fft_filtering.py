import MDAnalysis as mda
import numpy as np
import pytest
from numpy.testing import assert_allclose

from membrane_curvature.binning_surface import get_z_surface
from membrane_curvature.curvature import mean_curvature
from membrane_curvature.fft_filtering import (
    _height_for_fft,
    _squared_radius_grid,
    apply_fft_filter,
    nyquist_q,
    resolve_fft_filter,
)
from membrane_curvature.base import MembraneCurvature


@pytest.fixture()
def universe_dummy_wrap():
    positions = np.array(
        [
            [300.0, 0.0, 110.0],
            [100.0, 0.0, 150.0],
            [200.0, 0.0, 150.0],
            [0.0, 100.0, 150.0],
            [100.0, 100.0, 150.0],
            [-100.0, 100.0, 150.0],
            [300.0, 200.0, 110.0],
            [100.0, 200.0, 150.0],
            [200.0, 200.0, 150.0],
        ]
    )
    universe = mda.Universe(positions, n_atoms=9)
    universe.dimensions = [300, 300, 300, 90.0, 90.0, 90.0]
    return universe


def test_resolve_fft_filter_q():
    assert resolve_fft_filter({'q': (0.0, 0.5)}, 10.0, 10.0) == (0.0, 0.5)


def test_resolve_fft_filter_q_list_pair():
    assert resolve_fft_filter({'q': [0.0, 0.3]}, 10.0, 10.0) == (0.0, 0.3)


@pytest.mark.parametrize('dx, dy', [(0.0, 1.0), (1.0, -2.0)])
def test_nyquist_q_rejects_non_positive_bin_widths(dx, dy):
    with pytest.raises(ValueError, match='Bin widths must be positive'):
        nyquist_q(dx, dy)


def test_resolve_fft_filter_auto():
    dx, dy = 10.0, 20.0
    q_nyq = nyquist_q(dx, dy)
    assert resolve_fft_filter('auto', dx, dy) == (0.0, 0.5 * q_nyq)


@pytest.mark.parametrize(
    'fft_filter, match',
    [
        ({'q': (0.0,)}, 'must be a pair'),
        ({'r': (0.0, 0.3)}, 'dict of the form'),
        ({'q_high_frac': 0.5}, 'dict of the form'),
        ({'q': (0.5, 0.1)}, 'q_low must be <= q_high'),
        ({'q': (0.3, 0.3)}, 'must differ'),
        ({'q': (-0.1, 0.5)}, 'non-negative'),
        (True, 'dict of the form'),
        ('invalid', 'dict of the form'),
        (None, 'dict of the form'),
    ],
)
def test_resolve_fft_filter_errors(fft_filter, match):
    with pytest.raises(ValueError, match=match):
        resolve_fft_filter(fft_filter, 10.0, 10.0)


@pytest.mark.parametrize('n_bins', [3, 4, 5, 9, 32, 100])
def test_squared_radius_grid_matches_fftfreq(n_bins):
    dx = 2.5
    radius_sq = _squared_radius_grid(n_bins, n_bins, dx, dx)
    q_true = np.abs(np.fft.fftfreq(n_bins, d=dx) * (2.0 * np.pi))
    q_impl = np.sqrt(radius_sq[:, 0])
    assert_allclose(q_impl, q_true, rtol=0, atol=1e-12)


def test_resolve_fft_filter_warns_above_nyquist():
    dx, dy = 10.0, 10.0
    q_nyq = nyquist_q(dx, dy)
    with pytest.warns(UserWarning, match='exceeds Nyquist'):
        resolve_fft_filter({'q': (0.0, q_nyq * 2.0)}, dx, dy)


def test_resolve_fft_filter_warns_positive_q_low():
    with pytest.warns(UserWarning, match='mean height'):
        resolve_fft_filter({'q': (0.1, 0.5)}, 10.0, 10.0)


def test_resolve_fft_filter_rejects_q_low_above_nyquist():
    dx, dy = 10.0, 10.0
    q_nyq = nyquist_q(dx, dy)
    with pytest.raises(ValueError, match='exceeds the Nyquist limit'):
        resolve_fft_filter({'q': (q_nyq * 1.1, q_nyq * 2.0)}, dx, dy)


def test_constant_height_unchanged_by_low_pass():
    height = np.full((8, 8), 42.0)
    filtered = apply_fft_filter(height, 1.0, 1.0, (0.0, 0.5))
    assert_allclose(filtered, 42.0, rtol=1e-10)


def test_nan_bins_preserved():
    height = np.array([[1.0, np.nan], [3.0, 4.0]])
    filtered = apply_fft_filter(height, 1.0, 1.0, (0.0, 0.5))
    assert np.isnan(filtered[0, 1])
    assert not np.isnan(filtered[1, 0])


def test_height_for_fft_fills_with_finite_mean():
    height = np.array([[2.0, np.nan], [4.0, np.nan]])
    filled = _height_for_fft(height, np.isnan(height))
    assert filled[0, 1] == 3.0
    assert filled[1, 1] == 3.0


def test_height_for_fft_all_nan_uses_zero_fill():
    height = np.full((2, 2), np.nan)
    filled = _height_for_fft(height, np.isnan(height))
    assert_allclose(filled, 0.0)


def test_resolve_fft_filter_manual_and_auto_bounds():
    manual = (0.1, 0.2)
    assert resolve_fft_filter('auto', 10.0, 10.0) == (0.0, 0.5 * nyquist_q(10.0, 10.0))
    assert resolve_fft_filter({'q': manual}, 10.0, 10.0) == manual


def test_high_frequency_removed():
    x = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    y = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    xx, yy = np.meshgrid(x, y, indexing='ij')
    height = np.sin(3 * xx) + np.sin(12 * yy)
    dx = dy = 2 * np.pi / 32

    filtered = apply_fft_filter(height, dx, dy, (0.0, 0.15))
    residual = height - filtered
    assert np.std(residual) > np.std(filtered)


@pytest.fixture
def binning_3x3_fft_mc(universe_dummy_wrap):
    """Binning run with manual ``fft_filter`` on a 3x3 grid (single-frame universe)."""
    x_range = (0, 300)
    y_range = (0, 300)
    dx = (x_range[1] - x_range[0]) / 3
    q_high = 0.5 * nyquist_q(dx, dx)
    mc = MembraneCurvature(
        universe_dummy_wrap,
        n_x_bins=3,
        n_y_bins=3,
        surface_method='binning',
        fft_filter={'q': (0.0, q_high)},
    ).run()
    wrapped = universe_dummy_wrap.atoms.copy()
    wrapped.wrap()
    height = get_z_surface(wrapped.positions, n_x_bins=3, n_y_bins=3, x_range=x_range, y_range=y_range)
    return mc, q_high, height


def test_membrane_curvature_average_z_surface_is_filtered_time_mean(binning_3x3_fft_mc):
    mc, q_high, height = binning_3x3_fft_mc
    expected = apply_fft_filter(height, mc.dx, mc.dy, (0.0, q_high))
    assert_allclose(mc.results.average_z_surface, expected)


def test_membrane_curvature_per_frame_z_surface_not_fft_filtered(binning_3x3_fft_mc):
    mc, _q_high, height = binning_3x3_fft_mc
    assert_allclose(mc.results.z_surface[0], height)


def test_membrane_curvature_binning_filter_changes_average_surface(universe_dummy_wrap):
    mc_raw = MembraneCurvature(
        universe_dummy_wrap, n_x_bins=12, n_y_bins=12, surface_method='binning', fft_filter=None
    ).run()
    dx = 300.0 / 12
    q_high = 0.4 * nyquist_q(dx, dx)
    mc_filtered = MembraneCurvature(
        universe_dummy_wrap, n_x_bins=12, n_y_bins=12, surface_method='binning', fft_filter={'q': (0.0, q_high)}
    ).run()
    assert not np.allclose(mc_raw.results.average_z_surface, mc_filtered.results.average_z_surface, equal_nan=True)


def test_average_curvature_order_filter_then_curvature():
    """Curvature of filtered time-average differs from time-average of per-frame curvatures."""
    rng = np.random.default_rng(0)
    n_bins = 8
    dx = 1.0
    q_high = 0.4 * nyquist_q(dx, dx)
    heights = [rng.normal(size=(n_bins, n_bins)), rng.normal(size=(n_bins, n_bins))]
    z_avg = np.mean(np.stack(heights), axis=0)
    correct = mean_curvature(apply_fft_filter(z_avg, dx, dx, (0.0, q_high)), dx, dx)
    wrong = np.mean(
        [mean_curvature(apply_fft_filter(h, dx, dx, (0.0, q_high)), dx, dx) for h in heights],
        axis=0,
    )
    assert not np.allclose(correct, wrong)


@pytest.fixture
def two_frame_binning_auto_mc():
    """Two-frame binning run with default ``fft_filter='auto'`` on a 2x2 grid."""
    xy = np.array([[25.0, 25.0], [75.0, 25.0], [25.0, 75.0], [75.0, 75.0]])
    z_frame0 = np.array([10.0, 12.0, 11.0, 13.0])
    z_frame1 = np.array([14.0, 10.0, 13.0, 11.0])
    universe = mda.Universe(np.array([np.column_stack([xy, z_frame0]), np.column_stack([xy, z_frame1])]))
    universe.dimensions = [100.0, 100.0, 30.0, 90.0, 90.0, 90.0]
    return MembraneCurvature(
        universe,
        n_x_bins=2,
        n_y_bins=2,
        surface_method='binning',
        fft_filter='auto',
        wrap=False,
    ).run()


def test_membrane_curvature_average_z_surface_is_fft_of_time_mean(two_frame_binning_auto_mc):
    mc = two_frame_binning_auto_mc
    q_bounds = resolve_fft_filter(mc.fft_filter, mc.dx, mc.dy)
    z_avg = np.nanmean(mc.results.z_surface, axis=0)
    expected = apply_fft_filter(z_avg, mc.dx, mc.dy, q_bounds)
    assert_allclose(mc.results.average_z_surface, expected, rtol=1e-10, equal_nan=True)


def test_membrane_curvature_average_mean_from_filtered_time_average(two_frame_binning_auto_mc):
    mc = two_frame_binning_auto_mc
    q_bounds = resolve_fft_filter(mc.fft_filter, mc.dx, mc.dy)
    z_avg = np.nanmean(mc.results.z_surface, axis=0)
    z_filtered = apply_fft_filter(z_avg, mc.dx, mc.dy, q_bounds)
    assert_allclose(mc.results.average_mean, mean_curvature(z_filtered, mc.dx, mc.dy), rtol=1e-10, equal_nan=True)


def test_membrane_curvature_fourier_rejects_fft_filter_dict(universe_dummy_wrap):
    with pytest.raises(ValueError, match='fft_filter dict is only allowed'):
        MembraneCurvature(
            universe_dummy_wrap,
            surface_method='fourier',
            fourier_m=1,
            fourier_n=1,
            fft_filter={'q': (0.0, 0.2)},
        )


def test_membrane_curvature_fourier_without_fft_filter(universe_dummy_wrap):
    mc = MembraneCurvature(
        universe_dummy_wrap, n_x_bins=3, n_y_bins=3, surface_method='fourier', fourier_m=1, fourier_n=1
    ).run()
    assert mc.fft_filter is None


def test_membrane_curvature_binning_default_fft_filter_auto(universe_dummy_wrap):
    mc = MembraneCurvature(universe_dummy_wrap, n_x_bins=3, n_y_bins=3, surface_method='binning').run()
    assert mc.fft_filter == 'auto'
    q_nyq = nyquist_q(mc.dx, mc.dy)
    q_bounds = resolve_fft_filter(mc.fft_filter, mc.dx, mc.dy)
    assert q_bounds == (0.0, 0.5 * q_nyq)


def test_membrane_curvature_binning_fft_filter_none(universe_dummy_wrap):
    mc = MembraneCurvature(universe_dummy_wrap, n_x_bins=3, n_y_bins=3, surface_method='binning', fft_filter=None).run()
    assert mc.fft_filter is None
    z_avg = np.nanmean(mc.results.z_surface, axis=0)
    assert_allclose(mc.results.average_z_surface, z_avg, equal_nan=True)
