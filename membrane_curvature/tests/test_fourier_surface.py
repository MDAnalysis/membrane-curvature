"""
Unit tests for Fourier surface functions.
"""

from dataclasses import dataclass

import numpy as np
import pytest
import MDAnalysis as mda
from numpy.testing import assert_almost_equal

from membrane_curvature.fourier_surface import (
    _bin_centre_mesh,
    _build_fourier_matrix,
    _compute_wavevector,
    _harmonic_height_and_phi_derivatives,
    _eval_fourier_surface,
    _fourier_fit_from_atoms,
    _solve_design_least_squares_svd,
    n_fourier_parameters,
    fourier_mode_list,
    fourier_height_from_atoms,
    fourier_height_derivatives_from_atoms,
)


@dataclass(frozen=True)
class EvalFourierSurfaceCommonInputs:
    X_rel: np.ndarray
    Y_rel: np.ndarray
    A00: float
    A: float
    coeffs: dict
    Lx: float
    Ly: float


# This fixture is a dummy universe to check Fourier surface for
# a 2x2 grid with heights z(0,0)=10; z(0,1)=9; z(1,0)=9; z(1,1)=8
#   +-------+
#   |10 | 9 |
#   +-------+ # The z mean or A_00 is 9
#   | 9 | 8 |
#   +-------+
@pytest.fixture()
def dummy_fourier_universe():
    a = np.array(
        [
            [0.0, 0.0, 10.0],  # z(0,0)=10;
            [0.0, 1.0, 9.0],  # z(0,1)=9;
            [1.0, 0.0, 9.0],  # z(1,0)=9;
            [1.0, 1.0, 8.0],  # z(1,1)=8
        ]
    )
    u = mda.Universe(a, n_atoms=4)
    u.dimensions = [2.0, 2.0, 10.0, 90.0, 90.0, 90.0]
    return u


@pytest.fixture()
def dummy_fourier_fit(dummy_fourier_universe):
    """Fixture to obtain _box_length_x, _box_length_y, mean_term,
    and _coeffs from dummy Fourier universe"""
    positions = dummy_fourier_universe.atoms.positions.copy()
    return _fourier_fit_from_atoms(
        positions,
        x_range=(0.0, 2.0),
        y_range=(0.0, 2.0),
        M=1,
        N=0,
    )


@pytest.fixture()
def eval_fourier_surface_common_inputs(dummy_fourier_universe):
    x_positions = dummy_fourier_universe.atoms.positions[:, 0]
    y_positions = dummy_fourier_universe.atoms.positions[:, 1]
    return EvalFourierSurfaceCommonInputs(
        X_rel=x_positions.reshape(2, 2).astype(np.float64),
        Y_rel=y_positions.reshape(2, 2).astype(np.float64),
        A00=9.0,
        A=0.5,
        coeffs={(1, 0): (0.5, 0.0)},
        Lx=2.0,
        Ly=2.0,
    )


@pytest.fixture()
def dummy_two_mode_fourier_system(dummy_fourier_universe):
    positions = dummy_fourier_universe.atoms.positions
    x_positions = positions[:, 0]
    y_positions = positions[:, 1]
    z_values = positions[:, 2]
    modes = [(1, 0), (0, 1)]
    design_matrix = _build_fourier_matrix(
        x_positions,
        y_positions,
        2.0,
        2.0,
        modes,
    )
    return design_matrix, z_values


@pytest.fixture()
def full_rank_design_system():
    """Full-rank design with a known exact least-squares solution."""
    design_matrix = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float64,
    )
    theta_true = np.array([2.0, 3.0], dtype=np.float64)
    targets = design_matrix @ theta_true
    return design_matrix, targets, theta_true


@pytest.fixture()
def rank_deficient_design_system():
    """Duplicate columns -> rank 1; minimum-norm solution splits mean(targets) equally."""
    design_matrix = np.ones((3, 2), dtype=np.float64)
    targets = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    expected_theta = np.full(2, targets.mean() / 2.0)
    return design_matrix, targets, expected_theta


def test_calculate_a00_2x2_grid(dummy_fourier_universe, dummy_fourier_fit):
    """Test A_00 matches np.mean of z coords"""
    expected_A00 = 9  # equivalent to np.mean(positions[:, 2])
    _box_length_x, _box_length_y, calculated_A_00, _coeffs = dummy_fourier_fit
    assert_almost_equal(calculated_A_00, expected_A00)


def test_fourier_matrix_2x2_grid(dummy_fourier_universe):
    """Design matrix should match the expected [1, cos, sin] columns."""
    x_positions = dummy_fourier_universe.atoms.positions[:, 0]
    y_positions = dummy_fourier_universe.atoms.positions[:, 1]

    modes = [(1, 0)]
    box_length_x = 2.0
    box_length_y = 2.0

    design_matrix = _build_fourier_matrix(x_positions, y_positions, box_length_x, box_length_y, modes)

    expected_matrix = np.array(
        [
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
        ],
        dtype=np.float64,
    )

    assert_almost_equal(design_matrix, expected_matrix)


@pytest.mark.parametrize(
    'M, N, expected_n_modes',
    [
        (0, 0, 0),
        (2, 2, 12),
    ],
)
def test_fourier_mode_list_size(M, N, expected_n_modes):
    assert len(fourier_mode_list(M, N)) == expected_n_modes


@pytest.mark.parametrize(
    'M, N, expected_n_parameters',
    [
        (0, 0, 1),
        (2, 2, 25),
    ],
)
def test_n_fourier_parameters(M, N, expected_n_parameters):
    assert n_fourier_parameters(M, N) == expected_n_parameters


@pytest.mark.parametrize('M, N', [(-1, 0), (0, -1), (-1, -1)])
def test_fourier_fit_negative_mode_indices(dummy_fourier_universe, M, N):
    positions = dummy_fourier_universe.atoms.positions.copy()
    with pytest.raises(ValueError, match='M and N must be non-negative'):
        _fourier_fit_from_atoms(
            positions,
            x_range=(0.0, 2.0),
            y_range=(0.0, 2.0),
            M=M,
            N=N,
        )


def test_fourier_matrix_mismatched_xy_shapes(dummy_fourier_universe):
    x_positions = dummy_fourier_universe.atoms.positions[:, 0]
    y_positions = dummy_fourier_universe.atoms.positions[:3, 1]
    with pytest.raises(ValueError, match='x_positions and y_positions must have the same shape'):
        _build_fourier_matrix(x_positions, y_positions, 2.0, 2.0, [(1, 0)])


# Here we progressively test 3 different cases:
# Single mode (1,0): varies with x, constant in y
# Single mode (0,1): varies with y, constant in x
# Combined modes: both effects together.
@pytest.mark.parametrize(
    'modes, expected_matrix',
    [
        (
            [(1, 0)],  # single mode, varies in x
            np.array(
                [
                    [1.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, -1.0, 0.0],
                ],
                dtype=np.float64,
            ),
        ),
        (
            [(0, 1)],  # single mode, varies in y
            np.array(
                [
                    [1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, -1.0, 0.0],
                ],
                dtype=np.float64,
            ),
        ),
    ],
)
def test_fourier_matrix_single_modes_parametrized(dummy_fourier_universe, modes, expected_matrix):
    x_positions = dummy_fourier_universe.atoms.positions[:, 0]
    y_positions = dummy_fourier_universe.atoms.positions[:, 1]
    design_matrix = _build_fourier_matrix(x_positions, y_positions, 2.0, 2.0, modes)
    assert_almost_equal(design_matrix, expected_matrix)


def test_fourier_matrix_two_modes_reconstructs_2x2_grid(dummy_fourier_universe):
    """Test combining modes (0,1) and (1,0)"""
    x_positions = dummy_fourier_universe.atoms.positions[:, 0]
    y_positions = dummy_fourier_universe.atoms.positions[:, 1]
    z_values = dummy_fourier_universe.atoms.positions[:, 2]

    modes = [(1, 0), (0, 1)]
    design_matrix = _build_fourier_matrix(x_positions, y_positions, 2.0, 2.0, modes)
    expected_matrix = np.array(
        [
            [1.0, 1.0, 0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0, -1.0, 0.0],
            [1.0, -1.0, 0.0, 1.0, 0.0],
            [1.0, -1.0, 0.0, -1.0, 0.0],
        ],
        dtype=np.float64,
    )
    theta = np.array([9.0, 0.5, 0.0, 0.5, 0.0], dtype=np.float64)
    assert_almost_equal(design_matrix, expected_matrix)
    assert_almost_equal(design_matrix @ theta, z_values)


@pytest.mark.parametrize(
    'mode_x, mode_y, expected_kx, expected_ky',
    [
        (1, 0, np.pi, 0.0),
        (0, 1, 0.0, np.pi),
    ],
)
def test_compute_wavevector(mode_x, mode_y, expected_kx, expected_ky):
    """Computes wavevector and check T=2 with modes (1, 0) and (0,1)
    Solution gives kx=pi and ky=0"""
    kx_base = 2.0 * np.pi / 2.0
    ky_base = 2.0 * np.pi / 2.0
    kx, ky = _compute_wavevector(mode_x, mode_y, kx_base, ky_base)

    assert_almost_equal(kx, expected_kx)
    assert_almost_equal(ky, expected_ky)


@pytest.mark.parametrize('A, B', [(2.0, 3.0), (0.5, 0.5)])
def test_harmonic_height_and_phi_derivatives(dummy_fourier_universe, A, B):
    """Parametrizes for non-symmetric values, and known values for
    the dummy Fourier universe"""
    x_positions = dummy_fourier_universe.atoms.positions[:, 0]
    phase = np.pi * x_positions
    cos_phase = np.cos(phase)
    sin_phase = np.sin(phase)

    h, h_phi, h_phiphi = _harmonic_height_and_phi_derivatives(
        A=A,
        B=B,
        cos_phase=cos_phase,
        sin_phase=sin_phase,
    )

    expected_h = A * cos_phase + B * sin_phase
    expected_h_phi = -A * sin_phase + B * cos_phase
    expected_h_phiphi = -A * cos_phase - B * sin_phase
    assert_almost_equal(h, expected_h)
    assert_almost_equal(h_phi, expected_h_phi)
    assert_almost_equal(h_phiphi, expected_h_phiphi)


def test_solve_design_least_squares_svd_targets_shape_mismatch(full_rank_design_system):
    design_matrix, targets, _ = full_rank_design_system
    with pytest.raises(ValueError, match='targets must have shape \\(design_matrix.shape\\[0\\],\\)'):
        _solve_design_least_squares_svd(design_matrix, targets[:2])


def test_solve_design_least_squares_svd_empty_design_matrix_returns_zeros():
    design_matrix = np.zeros((0, 3), dtype=np.float64)
    targets = np.zeros(0, dtype=np.float64)
    theta, rank, singular_values = _solve_design_least_squares_svd(design_matrix, targets)
    assert rank == 0
    assert singular_values.size == 0
    assert_almost_equal(theta, np.zeros(3))


def test_solve_design_least_squares_svd_recovers_known_theta(full_rank_design_system):
    design_matrix, targets, theta_true = full_rank_design_system
    theta, rank, _ = _solve_design_least_squares_svd(design_matrix, targets)
    assert rank == design_matrix.shape[1]
    assert_almost_equal(theta, theta_true)


def test_solve_design_least_squares_svd_exact_fit_zero_residual(full_rank_design_system):
    design_matrix, targets, theta_true = full_rank_design_system
    theta, _, _ = _solve_design_least_squares_svd(design_matrix, targets)
    residual = design_matrix @ theta - targets
    assert_almost_equal(residual, np.zeros_like(targets))
    assert_almost_equal(theta, theta_true)


def test_solve_design_least_squares_svd_singular_values_match_numpy(full_rank_design_system):
    design_matrix, targets, _ = full_rank_design_system
    _, _, singular_values = _solve_design_least_squares_svd(design_matrix, targets)
    _, expected_singular_values, _ = np.linalg.svd(design_matrix, full_matrices=False)
    assert_almost_equal(singular_values, expected_singular_values)


def test_solve_design_least_squares_svd_rank_matches_truncation_mask(full_rank_design_system):
    design_matrix, targets, _ = full_rank_design_system
    rcond = 1e-12
    _, rank, singular_values = _solve_design_least_squares_svd(
        design_matrix,
        targets,
        rcond=rcond,
    )
    cutoff = rcond * singular_values[0]
    expected_rank = int(np.count_nonzero(singular_values > cutoff))
    assert rank == expected_rank


def test_solve_design_least_squares_svd_coerces_inputs_to_float64(full_rank_design_system):
    design_matrix, targets, theta_true = full_rank_design_system
    design_int = design_matrix.astype(np.int64)
    targets_int = targets.astype(np.int64)
    theta, rank, singular_values = _solve_design_least_squares_svd(design_int, targets_int)
    assert theta.dtype == np.float64
    assert singular_values.dtype == np.float64
    assert rank == design_matrix.shape[1]
    assert_almost_equal(theta, theta_true)


def test_solve_design_least_squares_svd_rank_deficient_minimum_norm(rank_deficient_design_system):
    design_matrix, targets, expected_theta = rank_deficient_design_system
    theta, rank, _ = _solve_design_least_squares_svd(design_matrix, targets, rcond=None)
    assert rank == 1
    assert rank < design_matrix.shape[1]
    assert_almost_equal(theta, expected_theta)


def test_solve_design_least_squares_svd_large_rcond_returns_zero_coefficients(rank_deficient_design_system):
    design_matrix, targets, _ = rank_deficient_design_system
    theta, rank, singular_values = _solve_design_least_squares_svd(
        design_matrix,
        targets,
        rcond=1.0,
    )
    assert rank == 0
    assert_almost_equal(theta, np.zeros(design_matrix.shape[1]))
    assert singular_values.size > 0


def test_solve_design_least_squares_svd_matches_lstsq_reference(dummy_two_mode_fourier_system):
    design_matrix, targets = dummy_two_mode_fourier_system
    theta_svd, _, _ = _solve_design_least_squares_svd(design_matrix, targets, rcond=None)
    theta_lstsq, _, _, _ = np.linalg.lstsq(design_matrix, targets, rcond=None)
    assert_almost_equal(theta_svd, theta_lstsq)


def test_fourier_warning_rank_deficient(dummy_fourier_universe):
    positions = dummy_fourier_universe.atoms.positions.copy()
    warning_regex = r'rank-deficient or underdetermined'

    with pytest.warns(UserWarning, match=warning_regex):
        _fourier_fit_from_atoms(
            positions,
            x_range=(0.0, 2.0),
            y_range=(0.0, 2.0),
            M=1,
            N=0,
        )


def test_fourier_height_from_atoms_range_width_non_positive():
    positions = np.array([[0.0, 0.0, 1.0]], dtype=np.float64)
    with pytest.raises(ValueError, match='x_range and y_range must have positive width'):
        fourier_height_from_atoms(
            positions,
            x_range=(0.0, 0.0),
            y_range=(0.0, 10.0),
            n_x_bins=4,
            n_y_bins=4,
            M=0,
            N=0,
        )


def test_eval_fourier_surface_derivatives_zeros(dummy_fourier_universe):
    x_positions = dummy_fourier_universe.atoms.positions[:, 0]
    y_positions = dummy_fourier_universe.atoms.positions[:, 1]
    X_rel = x_positions.reshape(2, 2).astype(np.float64)
    Y_rel = y_positions.reshape(2, 2).astype(np.float64)
    Z, fx, fy, fxx, fyy, fxy = _eval_fourier_surface(
        X_rel=X_rel,
        Y_rel=Y_rel,
        Lx=2.0,
        Ly=2.0,
        A00=9.0,
        coeffs={},
        derivatives=True,
    )
    assert_almost_equal(Z, np.full_like(X_rel, 9.0))
    assert_almost_equal(fx, np.zeros_like(X_rel))
    assert_almost_equal(fy, np.zeros_like(X_rel))
    assert_almost_equal(fxx, np.zeros_like(X_rel))
    assert_almost_equal(fyy, np.zeros_like(X_rel))
    assert_almost_equal(fxy, np.zeros_like(X_rel))


def test_eval_surface_derivatives_single_mode_x(eval_fourier_surface_common_inputs):
    """checks derivatives in _eval_fourier_surface: wrong sign in fx/fxx,
    swapping kx and ky, or wrong accumulation target for fxy and fxx"""
    inputs = eval_fourier_surface_common_inputs
    Z, fx, fy, fxx, fyy, fxy = _eval_fourier_surface(
        inputs.X_rel,
        inputs.Y_rel,
        inputs.Lx,
        inputs.Ly,
        inputs.A00,
        inputs.coeffs,
        derivatives=True,
    )
    phase = np.pi * inputs.X_rel
    expected_Z = inputs.A00 + inputs.A * np.cos(phase)
    expected_fx = -inputs.A * np.sin(phase) * np.pi
    expected_fy = np.zeros_like(inputs.X_rel)
    expected_fxx = -inputs.A * np.cos(phase) * (np.pi**2)
    expected_fyy = np.zeros_like(inputs.X_rel)
    expected_fxy = np.zeros_like(inputs.X_rel)
    assert_almost_equal(Z, expected_Z)
    assert_almost_equal(fx, expected_fx)
    assert_almost_equal(fy, expected_fy)
    assert_almost_equal(fxx, expected_fxx)
    assert_almost_equal(fyy, expected_fyy)
    assert_almost_equal(fxy, expected_fxy)


def test_bin_centre_mesh_returns_expected_centres():
    X_rel, Y_rel = _bin_centre_mesh(Lx=2.0, Ly=4.0, n_x_bins=2, n_y_bins=4)

    expected_X = np.array(
        [
            [0.5, 0.5, 0.5, 0.5],
            [1.5, 1.5, 1.5, 1.5],
        ],
        dtype=np.float64,
    )
    expected_Y = np.array(
        [
            [0.5, 1.5, 2.5, 3.5],
            [0.5, 1.5, 2.5, 3.5],
        ],
        dtype=np.float64,
    )

    assert X_rel.shape == (2, 4)
    assert Y_rel.shape == (2, 4)
    assert_almost_equal(X_rel, expected_X)
    assert_almost_equal(Y_rel, expected_Y)


@pytest.mark.parametrize('n_x_bins, n_y_bins', [(0, 2), (2, 0), (-1, 2), (2.5, 2)])
def test_bin_centre_mesh_error_bad_bin_counts(n_x_bins, n_y_bins):
    with pytest.raises(ValueError, match=r'_bin_centre_mesh requires positive bin counts') as exc_info:
        _bin_centre_mesh(Lx=2.0, Ly=4.0, n_x_bins=n_x_bins, n_y_bins=n_y_bins)
    assert isinstance(exc_info.value.__cause__, ValueError)


@pytest.mark.parametrize('derivatives', [False, True])
def test_eval_surface_derivatives_flag(eval_fourier_surface_common_inputs, derivatives):
    inputs = eval_fourier_surface_common_inputs

    out = _eval_fourier_surface(
        inputs.X_rel,
        inputs.Y_rel,
        inputs.Lx,
        inputs.Ly,
        inputs.A00,
        inputs.coeffs,
        derivatives=derivatives,
    )

    if derivatives:
        assert isinstance(out, tuple)
        assert len(out) == 6
        Z, fx, fy, fxx, fyy, fxy = out
        assert Z.shape == inputs.X_rel.shape
        assert fx.shape == inputs.X_rel.shape
        assert fy.shape == inputs.X_rel.shape
        assert fxx.shape == inputs.X_rel.shape
        assert fyy.shape == inputs.X_rel.shape
        assert fxy.shape == inputs.X_rel.shape
    else:
        assert isinstance(out, tuple)
        assert len(out) == 1
        assert out[0].shape == inputs.X_rel.shape


def test_fourier_height_from_atoms_calls_eval_without_derivatives(dummy_fourier_universe):
    positions = dummy_fourier_universe.atoms.positions.copy()

    z_grid = fourier_height_from_atoms(
        positions,
        x_range=(0.0, 2.0),
        y_range=(0.0, 2.0),
        n_x_bins=2,
        n_y_bins=2,
        M=1,
        N=0,
    )

    assert z_grid.shape == (2, 2)
    expected_surface = np.array([[9.0, 9.0], [9.0, 9.0]], dtype=np.float64)
    assert_almost_equal(z_grid, expected_surface)


def test_fourier_height_derivatives_from_atoms(dummy_fourier_universe):
    positions = dummy_fourier_universe.atoms.positions.copy()
    z_grid, fx, fy, fxx, fyy, fxy = fourier_height_derivatives_from_atoms(
        positions,
        x_range=(0.0, 2.0),
        y_range=(0.0, 2.0),
        n_x_bins=2,
        n_y_bins=2,
        M=1,
        N=0,
    )
    Lx, Ly, A00, coeffs = _fourier_fit_from_atoms(
        positions,
        x_range=(0.0, 2.0),
        y_range=(0.0, 2.0),
        M=1,
        N=0,
    )
    X_rel, Y_rel = _bin_centre_mesh(Lx, Ly, n_x_bins=2, n_y_bins=2)
    z_ref, fx_ref, fy_ref, fxx_ref, fyy_ref, fxy_ref = _eval_fourier_surface(
        X_rel, Y_rel, Lx, Ly, A00, coeffs, derivatives=True
    )

    assert z_grid.shape == (2, 2)
    assert fx.shape == (2, 2)
    assert fy.shape == (2, 2)
    assert fxx.shape == (2, 2)
    assert fyy.shape == (2, 2)
    assert fxy.shape == (2, 2)
    assert_almost_equal(z_grid, z_ref)
    assert_almost_equal(fx, fx_ref)
    assert_almost_equal(fy, fy_ref)
    assert_almost_equal(fxx, fxx_ref)
    assert_almost_equal(fyy, fyy_ref)
    assert_almost_equal(fxy, fxy_ref)


@pytest.mark.parametrize(
    'bad_positions',
    [
        np.empty((0, 3)),
        np.zeros(5),
        np.zeros((2, 2)),
    ],
)
def test_public_entry_point_error_bad_positions(bad_positions):
    with pytest.raises(ValueError, match='positions must'):
        fourier_height_from_atoms(
            bad_positions,
            x_range=(0.0, 1.0),
            y_range=(0.0, 1.0),
            n_x_bins=2,
            n_y_bins=2,
            M=0,
            N=0,
        )


@pytest.mark.parametrize('func', [fourier_height_from_atoms, fourier_height_derivatives_from_atoms])
def test_public_entry_points_error_nonpositive_bins(dummy_fourier_universe, func):
    positions = dummy_fourier_universe.atoms.positions.copy()
    with pytest.raises(ValueError, match='n_x_bins and n_y_bins must be positive'):
        func(
            positions,
            x_range=(0.0, 2.0),
            y_range=(0.0, 2.0),
            n_x_bins=0,
            n_y_bins=2,
            M=1,
            N=0,
        )
