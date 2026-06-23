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


@pytest.mark.parametrize(
    'fft_filter, expected',
    [
        ({'q': (0.0, 0.5)}, (0.0, 0.5)),
        ({'q': [0.0, 0.3]}, (0.0, 0.3)),
        ({'q': (0, 1)}, (0, 1)),
    ],
)
def test_resolve_fft_filter_q_variants(fft_filter, expected):
    assert resolve_fft_filter(fft_filter, 10.0, 10.0) == expected


@pytest.mark.parametrize('dx, dy', [(0.0, 1.0), (1.0, -2.0)])
def test_nyquist_q_rejects_non_positive_bin_widths(dx, dy):
    with pytest.raises(ValueError, match='Bin widths must be positive'):
        nyquist_q(dx, dy)


def test_resolve_fft_filter_auto():
    dx, dy = 10.0, 20.0
    q_nyq = nyquist_q(dx, dy)
    assert_allclose(resolve_fft_filter('auto', dx, dy), (0.0, 0.5 * q_nyq))


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
    assert_allclose(filled, 0.0, rtol=0)


def test_height_for_fft_does_not_mutate_input():
    height = np.array([[2.0, np.nan], [4.0, np.nan]])
    original = height.copy()
    _height_for_fft(height, np.isnan(height))
    assert_allclose(height, original, equal_nan=True)


def test_height_for_fft_no_nan_returns_input():
    height = np.array([[1.0, 2.0], [3.0, 4.0]])
    filled = _height_for_fft(height, np.isnan(height))
    assert filled is height


def test_resolve_fft_filter_manual_and_auto_bounds():
    manual = (0.1, 0.2)
    assert_allclose(resolve_fft_filter('auto', 10.0, 10.0), (0.0, 0.5 * nyquist_q(10.0, 10.0)))
    assert resolve_fft_filter({'q': manual}, 10.0, 10.0) == manual


def test_high_frequency_removed():
    """Low-pass filter keeps low q modes and removes high q modes."""
    x = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    y = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    xx, yy = np.meshgrid(x, y, indexing='ij')
    dx = dy = 2 * np.pi / 32
    low_freq = np.sin(xx)
    high_freq = np.sin(10 * yy)
    height = low_freq + high_freq
    filtered = apply_fft_filter(height, dx, dy, (0.0, 2.5))
    residual = height - filtered

    assert_allclose(filtered, low_freq, rtol=1e-10, atol=1e-10)
    assert_allclose(residual, high_freq, rtol=1e-10, atol=1e-10)


@pytest.fixture
def binning_3x3_fft_mc(universe_dummy_wrap):
    """Binning run with manual ``fft_filter`` on a 3x3 grid (single-frame universe)."""
    n_x_bins = n_y_bins = 3
    x_range = (0, 300)
    y_range = (0, 300)
    dx = (x_range[1] - x_range[0]) / n_x_bins
    dy = (y_range[1] - y_range[0]) / n_y_bins
    q_bounds = (0.0, 0.5 * nyquist_q(dx, dy))
    mc = MembraneCurvature(
        universe_dummy_wrap,
        n_x_bins=n_x_bins,
        n_y_bins=n_y_bins,
        surface_method='binning',
        fft_filter={'q': q_bounds},
    ).run()
    wrapped = universe_dummy_wrap.atoms.copy()
    wrapped.wrap()
    height = get_z_surface(
        wrapped.positions,
        n_x_bins=n_x_bins,
        n_y_bins=n_y_bins,
        x_range=x_range,
        y_range=y_range,
    )
    expected_average_z_surface = apply_fft_filter(height, dx, dy, q_bounds)
    return mc, dx, dy, q_bounds, height, expected_average_z_surface


def test_membrane_curvature_average_z_surface_is_filtered_time_mean(binning_3x3_fft_mc):
    mc, dx, dy, _q_bounds, _height, expected = binning_3x3_fft_mc
    assert mc.dx == dx
    assert mc.dy == dy
    assert_allclose(mc.results.average_z_surface, expected, rtol=1e-10, equal_nan=True)


def test_membrane_curvature_per_frame_z_surface_not_fft_filtered(binning_3x3_fft_mc):
    mc, _dx, _dy, _q_bounds, height, _expected = binning_3x3_fft_mc
    assert_allclose(mc.results.z_surface[0], height, rtol=1e-10, equal_nan=True)


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
    n_bins = 3
    dx = dy = universe_dummy_wrap.dimensions[0] / n_bins
    expected_q_bounds = resolve_fft_filter('auto', dx, dy)
    mc = MembraneCurvature(universe_dummy_wrap, n_x_bins=n_bins, n_y_bins=n_bins, surface_method='binning').run()
    assert mc.fft_filter == 'auto'
    assert mc.dx == dx
    assert mc.dy == dy
    assert_allclose(mc._fft_q_bounds, expected_q_bounds)


def test_membrane_curvature_binning_fft_filter_none(universe_dummy_wrap):
    # 3x3 binning of wrapped dummy positions - hardcoded expected surface
    # at z=110, all other occupied bins at z=150 as in universe_dummy_wrap
    expected = np.array(
        [
            [110.0, 150.0, 110.0],
            [150.0, 150.0, 150.0],
            [150.0, 150.0, 150.0],
        ]
    )

    mc = MembraneCurvature(universe_dummy_wrap, n_x_bins=3, n_y_bins=3, surface_method='binning', fft_filter=None).run()

    assert mc.fft_filter is None
    assert_allclose(mc.results.average_z_surface, expected, rtol=0)
