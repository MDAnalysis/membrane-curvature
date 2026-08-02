r"""

---------------
Fourier Surface
---------------

.. versionadded:: 2.0.0

With ``surface_method='fourier'``, :class:`~membrane_curvature.base.MembraneCurvature` derives a surface by
fitting a truncated 2D Fourier series to the :math:`z` coordinates of the reference atoms. The truncation is
set by the maximum mode indices along :math:`x` and :math:`y` with ``fourier_m`` and ``fourier_n``.

The fitted parameters are the Fourier coefficients: the mean term :math:`A_{00}`, and a cosine and a sine
amplitude, :math:`A_{mn}` and :math:`B_{mn}`, for each retained mode. ``MembraneCurvature`` evaluates the
basis functions at the atom positions to build a design matrix, with one row per atom and one column per
coefficient, and solves the resulting linear least squares system for the coefficients.

The fitted series is evaluated at the centre of every bin on a regular ``n_x_bins x n_y_bins`` grid to form the
height field. Because a Fourier series is a continuous function, every bin gets a value regardless of how many
atoms fall inside it, hence avoiding ``NaN`` values. Mean and Gaussian curvature are then calculated on the same
grid.

.. note::

    The Fourier method computes partial derivatives analytically from the fitted series.
    Unlike the binning methods, it does not rely on finite differences and therefore avoids artifacts
    at the edges of the box.

.. warning::

    Use ``fourier_m = fourier_n = 2``, the default values, unless you have a clear reason to resolve
    shorter wavelengths.

For more details on the Fourier method, refer to the :ref:`fourier_method` section in the
:ref:`algorithm` page.

References
----------

The Fourier method in MembraneCurvature uses the Fourier surface modeling [CAG2009]_
and molecular Fourier shape descriptors [JMG1988]_ as conceptual references:

.. [CAG2009] Shen et al., *Fourier method for large-scale surface modeling and registration*,
   Computers & Graphics (2009), doi: `10.1016/j.cag.2009.03.002`_.

.. [JMG1988] Leicester et al., *Description of molecular surface shape using Fourier descriptors*,
   Journal of Molecular Graphics (1988), doi: `10.1016/0263-7855(88)85008-2`_.

.. _`10.1016/j.cag.2009.03.002`: https://doi.org/10.1016/j.cag.2009.03.002
.. _`10.1016/0263-7855(88)85008-2`: https://doi.org/10.1016/0263-7855(88)85008-2

Functions
---------
"""

import warnings
from typing import Literal, overload

import numpy as np

from .fourier_validators import (
    _coerce_positions,
    validate_positions,
    validate_positive_bin_counts,
    validate_positive_domain_widths,
)


# Step 1: non-redundant Fourier modes
def fourier_mode_list(M: int, N: int) -> list[tuple[int, int]]:
    """
    Non-redundant 2D Fourier modes for a real-valued function on a periodic domain.

    Modes are :math:`(m, n)` where:

    - :math:`m=1..M`, :math:`n=-N..N` (all n retained, no symmetry to exploit)
    - :math:`m=0`, :math:`n=1..N` only (:math:`F(0,-n) = F(0,n)^*` makes n<0 redundant)

    The :math:`(0,0)` mean is excluded and handled separately as :math:`A_{00}`.

    Parameters
    ----------
    M : int
        Maximum positive mode index in :math:`x`.
    N : int
        Maximum mode index magnitude in :math:`y`.

    Returns
    -------
    modes : list of tuple of int
        Non-redundant ``(m, n)`` indices, excluding ``(0, 0)``.

    Raises
    ------
    ValueError
        If ``M`` or ``N`` is negative.
    """
    if M < 0 or N < 0:
        raise ValueError('M and N must be non-negative')

    # Modes where m > 0: no conjugate symmetry to exploit, so keep all n in -N..N.
    modes_with_m_positive = [(m, n) for m in range(1, M + 1) for n in range(-N, N + 1)]

    # Modes where m == 0: F(0, -n) = F(0, n)^* makes n < 0 redundant, and
    # n == 0 is the mean term A00 (handled separately), so keep only n in 1..N.
    modes_with_m_equal_zero = [(0, n) for n in range(1, N + 1)]

    return modes_with_m_positive + modes_with_m_equal_zero


def n_fourier_parameters(M: int, N: int) -> int:
    """
    Total number of scalar Fourier parameters for the truncated series.

    Uses the mode list from :func:`fourier_mode_list`: one scalar parameter for the mean
    term :math:`A_{00}` and two parameters per non-zero mode (cosine and sine
    coefficients). The total number of parameters is ``n = 1 + 2 * len(fourier_mode_list(M, N))``.

    Parameters
    ----------
    M : int
        Maximum positive mode index in :math:`x` (see :func:`fourier_mode_list`).
    N : int
        Maximum mode index magnitude in :math:`y`.

    Returns
    -------
    n : int
        Total number of scalar Fourier parameters ``1 + 2 * n_modes``.
    """
    n_modes = len(fourier_mode_list(M, N))
    return 1 + 2 * n_modes


# Step 2: wave numbers
def _compute_wavevector(mode_x, mode_y, kx_base, ky_base):
    r"""
    Wavevector components for a single mode.

    Builds :math:`(k_x, k_y) = (m\, k_{x,\mathrm{base}}, n\, k_{y,\mathrm{base}})`
    used in the phase :math:`k_x x + k_y y` for that mode.

    Parameters
    ----------
    mode_x : int
        Mode index :math:`m`.
    mode_y : int
        Mode index :math:`n`.
    kx_base : float
        Reciprocal step :math:`2\pi / L_x`.
    ky_base : float
        Reciprocal step :math:`2\pi / L_y`.

    Returns
    -------
    kx : float
        ``mode_x * kx_base``.
    ky : float
        ``mode_y * ky_base``.
    """
    return mode_x * kx_base, mode_y * ky_base


# Step 3: build design matrix
def _build_fourier_matrix(
    x_positions: np.ndarray,
    y_positions: np.ndarray,
    box_length_x: float,
    box_length_y: float,
    fourier_modes: list[tuple[int, int]],
) -> np.ndarray:
    r"""
    Build the least-squares design matrix from the Fourier basis functions.

    Each row is an atom position :math:`(x_i, y_i)` and each column is a basis function
    evaluated at that position: a constant column for :math:`A_{00}`, then a
    :math:`\cos(k \cdot r)` and :math:`\sin(k \cdot r)` pair per mode, with wavevectors
    :math:`k = (2\pi m / L_x, 2\pi n / L_y)`.

    Parameters
    ----------
    x_positions : ndarray
        :math:`x` coordinates of points within one periodic box (same shape as ``y_positions``).
    y_positions : ndarray
        :math:`y` coordinates of points within one periodic box.
    box_length_x : float
        Period length :math:`L_x`.
    box_length_y : float
        Period length :math:`L_y`.
    fourier_modes : list of tuple of int
        Non-redundant ``(m, n)`` indices, e.g. from :func:`fourier_mode_list`.

    Returns
    -------
    design_matrix : ndarray, shape (n_points, 1 + 2 * n_modes)
        Basis evaluated at each point. ``dtype`` is ``numpy.float64``.

    Raises
    ------
    ValueError
        If ``x_positions`` and ``y_positions`` do not have the same shape.
    """
    x_positions = np.asarray(x_positions, dtype=np.float64)
    y_positions = np.asarray(y_positions, dtype=np.float64)
    if x_positions.shape != y_positions.shape:
        raise ValueError('x_positions and y_positions must have the same shape')

    n_points = x_positions.shape[0]
    n_modes = len(fourier_modes)

    design_matrix = np.empty((n_points, 1 + 2 * n_modes), dtype=np.float64)
    # Column 0: constant offset for A00
    design_matrix[:, 0] = 1.0

    kx_base = 2.0 * np.pi / box_length_x
    ky_base = 2.0 * np.pi / box_length_y

    column_index = 1

    # Add a plane wave with wavevector k_mn - all atom positions
    for mode_x, mode_y in fourier_modes:
        kx, ky = _compute_wavevector(mode_x, mode_y, kx_base, ky_base)
        phase = kx * x_positions + ky * y_positions
        design_matrix[:, column_index] = np.cos(phase)
        design_matrix[:, column_index + 1] = np.sin(phase)
        column_index += 2

    return design_matrix


def _solve_design_least_squares_svd(
    design_matrix: np.ndarray,
    targets: np.ndarray,
    rcond: float | None = None,
) -> tuple[np.ndarray, int, np.ndarray]:
    r"""
    Least-squares solution with :func:`numpy.linalg.svd` and singular-value truncation.

    Solves ``min ||design_matrix @ theta - targets||_2`` with a relative singular-value
    cutoff ``rcond`` (values :math:`s \le rcond \cdot s_{\max}` are treated as
    zero). When the system is rank-deficient, returns the minimum Euclidean-norm
    coefficient vector among least-squares minimizers.

    Parameters
    ----------
    design_matrix : ndarray, shape (n_rows, n_columns)
        Design matrix with rows corresponding to atom positions and columns
        corresponding to the basis functions.
    targets : ndarray, shape (n_rows,)
        Target heights.
    rcond : float or None, optional
        Relative cutoff for small singular values. ``None`` uses a heuristic
        cutoff based on machine precision and the problem size.

    Returns
    -------
    theta : ndarray, shape (n_columns,)
        Fitted coefficients, one per column of the design matrix. The number of columns
        is the number of Fourier parameters.
    rank : int
        Number of singular values above the cutoff.
    singular_values : ndarray
        Singular values of ``design_matrix``.

    Raises
    ------
    ValueError
        If ``targets`` does not have shape ``(n_rows,)``.
    """
    design_matrix = np.asarray(design_matrix, dtype=np.float64)
    targets = np.asarray(targets, dtype=np.float64)
    n_rows, n_columns = design_matrix.shape
    if targets.shape != (n_rows,):
        raise ValueError('targets must have shape (design_matrix.shape[0],)')

    if rcond is None:
        rcond = np.finfo(np.float64).eps * max(n_rows, n_columns)

    unitary_array, singular_values, vh = np.linalg.svd(design_matrix, full_matrices=False)

    if singular_values.size == 0:
        return np.zeros(n_columns, dtype=np.float64), 0, singular_values

    cutoff = rcond * singular_values[0]
    keep = singular_values > cutoff
    rank = int(np.count_nonzero(keep))
    if rank == 0:
        return np.zeros(n_columns, dtype=np.float64), 0, singular_values

    projection = np.dot(unitary_array[:, keep].T, targets)
    coefficients_on_span = projection / singular_values[keep]
    theta = vh[keep].T @ coefficients_on_span
    return theta, rank, singular_values


# Step 4: least squares for Fourier coefficients
def _unpack_coefficients(
    theta: np.ndarray, modes: list[tuple[int, int]]
) -> tuple[float, dict[tuple[int, int], tuple[float, float]]]:
    """
    Split the least-squares coefficient vector into the mean and per-mode amplitudes.

    Index ``0`` is :math:`A_{00}`. Remaining entries are cosine/sine pairs in ``modes`` order.

    Parameters
    ----------
    theta : ndarray, shape (n_coefficients,)
        Least-squares solution vector: mean at index 0, then pairs ``(A_{mn}, B_{mn})``
        in the same order as ``modes``.
    modes : list of tuple of int
        Mode list from :func:`fourier_mode_list` (or equivalent ordering).

    Returns
    -------
    A00 : float
        Constant (mean) term.
    coeffs : dict
        Mapping ``(m, n) -> (A, B)`` cosine and sine coefficients per mode.
    """
    A00 = float(theta[0])
    pairs = theta[1:].reshape(-1, 2)
    coeffs = {mode: (float(a), float(b)) for mode, (a, b) in zip(modes, pairs)}
    return A00, coeffs


def _fourier_fit_from_atoms(
    positions: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    M: int,
    N: int,
    rcond: float | None = None,
) -> tuple[float, float, float, dict[tuple[int, int], tuple[float, float]]]:
    r"""
    Fit Fourier coefficients to atom heights, without evaluating on a grid.

    Builds the design matrix, solves for :math:`A_{00}` and the mode amplitudes, and
    returns the domain sizes and coefficients needed to evaluate the surface on a grid.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atom coordinates in the same length unit as the box. The third
        column (``positions[:, 2]``) is the height :math:`z`.
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    M, N : int
        Maximum Fourier mode indices (see :func:`fourier_mode_list`).
    rcond : float or None, optional
        Relative cutoff for small singular values in :func:`_solve_design_least_squares_svd`.

    Returns
    -------
    Lx : float
        Period length :math:`L_x = x_{\max} - x_{\min}`.
    Ly : float
        Period length :math:`L_y = y_{\max} - y_{\min}`.
    A00 : float
        Mean height coefficient :math:`A_{00}`.
    coeffs : dict[tuple[int, int], tuple[float, float]]
        Mapping ``(m, n) -> (A, B)`` cosine and sine coefficients per mode.

    Raises
    ------
    ValueError
        If ``positions`` cannot be converted to a float64 array or has the wrong shape,
        if ``x_range`` or ``y_range`` have non-positive width, or if ``M`` or ``N`` are
        negative.

    Warns
    -----
    UserWarning
        If the least-squares system is rank-deficient or underdetermined.
    """
    positions = _coerce_positions(positions)
    validate_positions(positions)
    validate_positive_domain_widths(x_range, y_range)
    x0, y0 = x_range[0], y_range[0]
    Lx = float(x_range[1] - x_range[0])
    Ly = float(y_range[1] - y_range[0])

    modes = fourier_mode_list(M, N)
    n_atom = len(positions)

    x_rel = np.mod(positions[:, 0] - x0, Lx)
    y_rel = np.mod(positions[:, 1] - y0, Ly)
    z_at = positions[:, 2]

    design_matrix = _build_fourier_matrix(x_rel, y_rel, Lx, Ly, modes)
    theta, rank, _singular_values = _solve_design_least_squares_svd(design_matrix, z_at, rcond=rcond)
    n_columns = int(design_matrix.shape[1])
    if rank < n_columns:
        warnings.warn(
            'Fourier least-squares system is rank-deficient or underdetermined: '
            f'effective SVD rank is {rank} for {n_columns} parameters and {n_atom} atoms '
            f'(M={M}, N={N}). Coefficients may be non-unique or ill-conditioned.',
            UserWarning,
            stacklevel=2,
        )
    A00, coeffs = _unpack_coefficients(theta, modes)
    return Lx, Ly, A00, coeffs


def _harmonic_height_and_phi_derivatives(
    A: float,
    B: float,
    cos_phase: np.ndarray,
    sin_phase: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Height and :math:`\phi`-derivatives for one mode :math:`h(\phi)=A\cos\phi+B\sin\phi`.

    Derivatives are taken with respect to the phase :math:`\phi`, not with respect to
    :math:`x` and :math:`y`. :func:`_eval_fourier_surface` turns them into spatial
    derivatives with the chain rule, multiplying by :math:`k_x` and :math:`k_y`.

    Parameters
    ----------
    A, B : float
        Cosine and sine coefficients for this mode.
    cos_phase : ndarray
        :math:`\cos\phi`, same shape as ``sin_phase``.
    sin_phase : ndarray
        :math:`\sin\phi`.

    Returns
    -------
    h : ndarray
        Contribution :math:`A\cos\phi + B\sin\phi`.
    h_phi : ndarray
        First derivative :math:`\partial h/\partial \phi`.
    h_phiphi : ndarray
        Second derivative :math:`\partial^2 h/\partial \phi^2`.
    """
    h = A * cos_phase + B * sin_phase
    h_phi = -A * sin_phase + B * cos_phase
    h_phiphi = -A * cos_phase - B * sin_phase
    return h, h_phi, h_phiphi


@overload
def _eval_fourier_surface(
    X_rel: np.ndarray,
    Y_rel: np.ndarray,
    Lx: float,
    Ly: float,
    A00: float,
    coeffs: dict[tuple[int, int], tuple[float, float]],
    *,
    derivatives: Literal[False],
) -> tuple[np.ndarray]: ...


@overload
def _eval_fourier_surface(
    X_rel: np.ndarray,
    Y_rel: np.ndarray,
    Lx: float,
    Ly: float,
    A00: float,
    coeffs: dict[tuple[int, int], tuple[float, float]],
    *,
    derivatives: Literal[True],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]: ...


def _eval_fourier_surface(
    X_rel: np.ndarray,
    Y_rel: np.ndarray,
    Lx: float,
    Ly: float,
    A00: float,
    coeffs: dict[tuple[int, int], tuple[float, float]],
    *,
    derivatives: bool,
) -> tuple[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Evaluate the fitted Fourier sum on a relative-coordinate grid.

    Reconstructs :math:`Z(x,y) = A_{00} + \sum_{(m,n)} A_{mn}\cos\phi_{mn}
    + B_{mn}\sin\phi_{mn}` with :math:`\phi_{mn} = k_x x + k_y y`. When ``derivatives``
    is ``True``, also accumulates the analytic first and second partial derivatives of
    the same sum, used for Monge-gauge curvature.

    Parameters
    ----------
    X_rel : ndarray
        Relative :math:`x` coordinates on the periodic domain (same shape as ``Y_rel``).
    Y_rel : ndarray
        Relative :math:`y` coordinates on the periodic domain.
    Lx : float
        Period length in :math:`x` (same units as ``X_rel``).
    Ly : float
        Period length in :math:`y`.
    A00 : float
        Mean term :math:`A_{00}`.
    coeffs : dict[tuple[int, int], tuple[float, float]]
        Mode amplitudes ``(A, B)`` from :func:`_fourier_fit_from_atoms` / :func:`_unpack_coefficients`.
    derivatives : bool
        If ``False``, return ``(Z,)``. If ``True``, return ``(Z, fx, fy, fxx, fyy, fxy)``.

    Returns
    -------
    tuple of ndarray
        If ``derivatives`` is ``False``, the 1-tuple ``(Z,)`` where ``Z`` matches the shape of
        ``X_rel``. If ``derivatives`` is ``True``, the 6-tuple with each grid matching ``X_rel``.

    Notes
    -----
    The ``derivatives=False`` branch returns ``(Z,)``, a 1-tuple and not a bare
    ``ndarray``, to keep a single return type (``tuple[ndarray, ...]``)
    regardless of the flag. Callers therefore unwrap with ``[0]``.
    """
    # TODO: the 1-tuple return when derivatives=False is awkward (callers do
    # [0] to unwrap). Future refactor options, in preference order:
    #   (a) return a bare ndarray for derivatives=False and a 6-tuple for
    #       derivatives=True; the existing @overload pair already expresses
    #       this union and removes the [0] unwrap at every call site;
    #   (b) split into two functions (no flag, no overloads, no union),
    #       e.g. _eval_fourier_height and _eval_fourier_height_with_derivatives;
    #   (c) return a NamedTuple/dataclass for field-name access (Z, fx, fy,
    #       fxx, fyy, fxy) while still supporting tuple unpacking.
    twopi = 2.0 * np.pi
    X_rel = np.asarray(X_rel, dtype=np.float64)
    Y_rel = np.asarray(Y_rel, dtype=np.float64)
    Z = np.full_like(X_rel, A00, dtype=np.float64)
    if derivatives:
        fx = np.zeros_like(X_rel, dtype=np.float64)
        fy = np.zeros_like(X_rel, dtype=np.float64)
        fxx = np.zeros_like(X_rel, dtype=np.float64)
        fyy = np.zeros_like(X_rel, dtype=np.float64)
        fxy = np.zeros_like(X_rel, dtype=np.float64)

    for (m, n), (A, B) in coeffs.items():
        kx = twopi * m / Lx
        ky = twopi * n / Ly
        phase = kx * X_rel + ky * Y_rel
        cos_phase = np.cos(phase)
        sin_phase = np.sin(phase)
        if derivatives:
            z_c, h_phi, h_phiphi = _harmonic_height_and_phi_derivatives(A, B, cos_phase, sin_phase)
            Z += z_c
            fx += h_phi * kx
            fy += h_phi * ky
            fxx += h_phiphi * (kx * kx)
            fyy += h_phiphi * (ky * ky)
            fxy += h_phiphi * (kx * ky)
        else:
            Z += A * cos_phase + B * sin_phase

    if derivatives:
        return Z, fx, fy, fxx, fyy, fxy
    return (Z,)


def _bin_centre_mesh(Lx: float, Ly: float, n_x_bins: int, n_y_bins: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Relative-coordinate mesh at bin centres.

    Bin centres at ``(i + 1/2) L / n_bins`` for each axis, compatible with the
    periodic domain ``[0, L)``.

    Parameters
    ----------
    Lx : float
        Period length in :math:`x`.
    Ly : float
        Period length in :math:`y`.
    n_x_bins : int
        Number of bins along :math:`x`.
    n_y_bins : int
        Number of bins along :math:`y`.

    Returns
    -------
    X_rel : ndarray
        Mesh of relative :math:`x` coordinates, shape ``(n_x_bins, n_y_bins)``.
    Y_rel : ndarray
        Mesh of relative :math:`y` coordinates, same shape as ``X_rel``.

    Raises
    ------
    ValueError
        If ``n_x_bins`` or ``n_y_bins`` is not a positive integer. The public entry points
        already validate this with
        :func:`~membrane_curvature.fourier_validators.validate_positive_bin_counts`.
    """
    # Note: the wrapper stays here (not in fourier_validators) because the
    # message names this specific helper. Validators are caller-agnostic; this
    # caller-specific context belongs with the caller!
    try:
        validate_positive_bin_counts(n_x_bins, n_y_bins)
    except ValueError as exc:
        raise ValueError(
            '_bin_centre_mesh requires positive bin counts. reached with '
            f'n_x_bins={n_x_bins}, n_y_bins={n_y_bins}. '
            'Callers must validate via validate_positive_bin_counts first.'
        ) from exc
    xc = (np.arange(n_x_bins, dtype=np.float64) + 0.5) * Lx / n_x_bins
    yc = (np.arange(n_y_bins, dtype=np.float64) + 0.5) * Ly / n_y_bins
    X_rel, Y_rel = np.meshgrid(xc, yc, indexing='ij')
    return X_rel, Y_rel


# Steps 1 to 5, or 1 to 6 if derivatives are included
def fourier_height_from_atoms(
    positions: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    n_x_bins: int,
    n_y_bins: int,
    M: int,
    N: int,
    rcond: float | None = None,
) -> np.ndarray:
    """
    Fourier height field on bin centres.

    Fits the Fourier coefficients to the atom heights with :func:`_fourier_fit_from_atoms`,
    then evaluates the fitted series at bin centres.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atom coordinates in the same length unit as the box. The third
        column (``positions[:, 2]``) is the height :math:`z`.
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    n_x_bins, n_y_bins : int
        Output grid size at bin centres.
    M, N : int
        Maximum Fourier mode indices, see :func:`fourier_mode_list`.
    rcond : float or None, optional
        Relative cutoff for small singular values in
        :func:`_solve_design_least_squares_svd`.

    Returns
    -------
    Z : ndarray, shape (n_x_bins, n_y_bins)
        Fitted surface height.

    Raises
    ------
    ValueError
        If ``n_x_bins`` or ``n_y_bins`` is not a positive integer, if ``positions`` cannot
        be converted to a float64 array or has the wrong shape, if ``x_range`` or
        ``y_range`` have non-positive width, or if ``M`` or ``N`` are negative.

    Warns
    -----
    UserWarning
        If the least-squares system is rank-deficient or underdetermined.
    """
    validate_positive_bin_counts(n_x_bins, n_y_bins)
    Lx, Ly, A00, coeffs = _fourier_fit_from_atoms(
        positions,
        x_range,
        y_range,
        M,
        N,
        rcond=rcond,
    )
    X_rel, Y_rel = _bin_centre_mesh(Lx, Ly, n_x_bins, n_y_bins)
    return _eval_fourier_surface(X_rel, Y_rel, Lx, Ly, A00, coeffs, derivatives=False)[0]


def fourier_height_derivatives_from_atoms(
    positions: np.ndarray,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    n_x_bins: int,
    n_y_bins: int,
    M: int,
    N: int,
    rcond: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Fourier height and analytic partial derivatives on bin centres.

    Same pipeline as :func:`fourier_height_from_atoms`, and also evaluates the analytic
    first and second partial derivatives of the fitted series.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atom coordinates. The third column (``positions[:, 2]``) is
        the height :math:`z`.
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    n_x_bins, n_y_bins : int
        Output grid size at bin centres.
    M, N : int
        Maximum Fourier mode indices, see :func:`fourier_mode_list`.
    rcond : float or None, optional
        Relative cutoff for small singular values in
        :func:`_solve_design_least_squares_svd`.

    Returns
    -------
    Z, fx, fy, fxx, fyy, fxy : ndarray
        Each array has shape ``(n_x_bins, n_y_bins)``.

    Raises
    ------
    ValueError
        If ``n_x_bins`` or ``n_y_bins`` is not a positive integer, if ``positions`` cannot
        be converted to a float64 array or has the wrong shape, if ``x_range`` or
        ``y_range`` have non-positive width, or if ``M`` or ``N`` are negative.

    Warns
    -----
    UserWarning
        If the least-squares system is rank-deficient or underdetermined.
    """
    validate_positive_bin_counts(n_x_bins, n_y_bins)
    Lx, Ly, A00, coeffs = _fourier_fit_from_atoms(
        positions,
        x_range,
        y_range,
        M,
        N,
        rcond=rcond,
    )
    X_rel, Y_rel = _bin_centre_mesh(Lx, Ly, n_x_bins, n_y_bins)
    return _eval_fourier_surface(X_rel, Y_rel, Lx, Ly, A00, coeffs, derivatives=True)
