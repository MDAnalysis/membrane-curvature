r"""

---------------
Fourier Surface
---------------

.. versionadded:: 2.0.0

.. note::

    This is an alternative to binning and calculating :func:`numpy.gradient`
    to derive surfaces and compute mean and Gaussian curvature on a
    simulation box with periodic boundary conditions in :math:`x` and
    :math:`y`.

This module fits a smooth periodic surface to atomistic or bead height data using
a truncated 2D Fourier series. The workflow includes six steps:

1. Choose a non-redundant set of Fourier modes:
    - :func:`fourier_mode_list`
    - :func:`n_fourier_parameters`
2. Compute the corresponding wave numbers:
    - :func:`_compute_wavevector`
3. Build a design matrix from cosine and sine basis functions
    - :func:`_build_fourier_matrix`
4. Solve the least-squares problem for the Fourier coefficients
    - :func:`_fourier_fit_from_atoms`
    - :func:`_unpack_coefficients`
5. Reconstruct the continuous surface on a grid:
    - :func:`_bin_centre_mesh`
    - :func:`_eval_fourier_surface`
6. Evaluate partial derivatives analytically:
    - :func:`_eval_fourier_surface` with ``derivatives=True``, or
    - :func:`fourier_height_derivatives_from_atoms`


Public entry points :func:`fourier_height_from_atoms` and
:func:`fourier_height_derivatives_from_atoms` chain these steps on bin centres.

Note that the use of a non-redundant set of Fourier modes is justified by the
conjugate symmetry of Fourier coefficients for real-valued functions:

.. math::

    F(-m, -n) = F(m, n)^*

We therefore keep only the non-redundant half of the mode grid. The :math:`m = 0`
row has the additional symmetry

.. math::

    F(0, -n) = F(0, n)^*

so we keep only :math:`n \ge 0` there and store the mean term :math:`A_{00}`
separately.

A diagram of a non-redundant Fourier mode grid looks like::

    n
    ^
    N  |  ✓  |  ✓  ✓  ✓ ... ✓  |
    :  |  ✓  |  ✓  ✓  ✓ ... ✓  |
    2  |  ✓  |  ✓  ✓  ✓ ... ✓  |
    1  |  ✓  |  ✓  ✓  ✓ ... ✓  |
    0  | A00 |  ✓  ✓  ✓ ... ✓  |   <- right block includes n=0; left cell is A00
    -1 |  x  |  ✓  ✓  ✓ ... ✓  |
    -2 |  x  |  ✓  ✓  ✓ ... ✓  |
    :  |  x  |  ✓  ✓  ✓ ... ✓  |
    -N |  x  |  ✓  ✓  ✓ ... ✓  |
       +-----+-------------------+--> m
          0     1   2   3  ...  M
          ^           ^
          |           |
      m_zero_modes  m_positive_modes
      (n=1..N only) (all n, -N..N)

where the ✓ symbol denotes that the wave includes 2 parameters each:
:math:`cos`, and :math:`sin` coefficients, and where the x symbols
denotes that the wave is excluded to avoid redundancy by
:math:`F(0,-n) = F(0,n)^*`, only occurs at :math:`m=0`.

The associated wavevector for each mode is

    .. math::

        \mathbf{k}_{mn} = \left(\frac{2\pi m}{L_x}, \frac{2\pi n}{L_y}\right)


.. important::

    The simulation box defines the fundamental periodic domain and therefore
    sets the largest allowed wavelength in the system.
    The surface is represented as a sum of periodic waves that are compatible
    with these boundary conditions, meaning all basis functions must fit
    exactly within the box without discontinuities at the boundaries.

    As a result, only discrete Fourier modes are allowed, corresponding to
    integer wave numbers:

    .. math::

        k_x = \frac{2\pi m}{L_x}, \qquad k_y = \frac{2\pi n}{L_y}

    with :math:`m, n \in \mathbb{Z}`.

    The box size :math:`L_x` and :math:`L_y` therefore determines the spacing of
    allowed wavevectors, while higher modes represent progressively shorter
    wavelengths.


References
----------

The Fourier surface representation follows standard Fourier descriptor and
Fourier surface-model approaches for geometric shape modeling [CAG2009]_ [JMG1988]_.

.. [CAG2009] Rhouma et al., *Fourier method for large scale surface modeling and registration*,
   Computers & Graphics (2009), doi: `10.1016/j.cag.2009.03.002`_.

.. [JMG1988] *Description of molecular surface shape using Fourier descriptors*,
   Journal of Molecular Graphics (1988), doi: `10.1016/0263-7855(88)85008-2`_.

.. _`10.1016/j.cag.2009.03.002`: https://doi.org/10.1016/j.cag.2009.03.002
.. _`10.1016/0263-7855(88)85008-2`: https://doi.org/10.1016/0263-7855(88)85008-2

Functions
---------
"""

import warnings

import numpy as np


# Step 1: non-redundant Fourier modes
def fourier_mode_list(M: int, N: int) -> list[tuple[int, int]]:
    """
    Non-redundant 2D Fourier modes for a real-valued function on a periodic domain.

    Implements step 1 of the module workflow (non-redundant mode grid).

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
    """
    if M < 0 or N < 0:
        raise ValueError('M and N must be non-negative')

    # m > 0: no symmetry, keep all n
    m_positive_modes = [(m, n) for m in range(1, M + 1) for n in range(-N, N + 1)]

    # m = 0: drop n < 0 for non-redundant, n = 0 (stored as A00)
    m_zero_modes = [(0, n) for n in range(1, N + 1)]

    return m_positive_modes + m_zero_modes


def n_fourier_parameters(M: int, N: int) -> int:
    """
    Total number of scalar Fourier parameters for the truncated series.

    Uses the same mode list as step 1 (:func:`fourier_mode_list`): one scalar
    parameter for the mean term :math:`A_{00}` and two parameters per non-zero mode
    (cosine and sine coefficients), so ``n = 1 + 2 * len(fourier_mode_list(M, N))``.

    .. warning::

        Use :math:`M = N = 2` unless you have a clear reason to resolve shorter wavelengths;
        always ensure your leaflet has many more atoms than Fourier parameters
        (:func:`~membrane_curvature.fourier_surface.n_fourier_parameters`), and increase
        :math:`M` and :math:`N` only until curvature stops changing systematically and
        starts chasing noise.

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
    Wavevector components for a single mode (step 2 of the module workflow).

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
    Build the least-squares design matrix: basis functions evaluated at sample points.

    Implements step 3 (and uses :func:`_compute_wavevector` for step 2 inside the loop).
    Each row is an atom position; columns are ``1``, then ``cos(k\cdot r)``, ``sin(k\cdot r)``
    pairs for each mode, so that fitting recovers :math:`z(x,y)` at the atoms.

    In the matrix, each row corresponds to a spatial point :math:`(x_i, y_i)`, and each
    column corresponds to a basis function:

    .. math::

        [1,
         cos(k·r), sin(k·r),
         cos(k·r), sin(k·r),
         ...]

    where the wave vectors are given by :math:`k=(2\pi m/L_x, 2\pi n/L_y)`.

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
        Basis evaluated at each point; ``dtype`` is ``numpy.float64``.
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


# Step 4: least squares for Fourier coefficients
def _unpack_coefficients(
    theta: np.ndarray, modes: list[tuple[int, int]]
) -> tuple[float, dict[tuple[int, int], tuple[float, float]]]:
    """
    Split the least-squares coefficient vector into the mean and per-mode amplitudes.

    Used after step 4 (:func:`numpy.linalg.lstsq` in :func:`_fourier_fit_from_atoms`).
    Index ``0`` is :math:`A_{00}`; remaining entries are cosine/sine pairs in ``modes`` order.

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
    Fit Fourier coefficients to atom heights (steps 1-4; no bin grid).

    Builds the design matrix, solves for :math:`A_{00}` and the mode amplitudes, and
    returns domain sizes and coefficients for subsequent grid evaluation (steps 5-6).

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atom coordinates; column 2 is height :math:`z`.
    x_range, y_range : tuple of float
        ``(min, max)`` extents of the periodic domain.
    M, N : int
        Maximum Fourier mode indices (see :func:`fourier_mode_list`).
    rcond : float or None, optional
        Passed to :func:`numpy.linalg.lstsq`.

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

    Warns
    -----
    UserWarning
        If the design matrix is rank-deficient or underdetermined (same condition as
        :func:`fourier_height_from_atoms` and :func:`fourier_height_derivatives_from_atoms`).
    """
    x0, y0 = x_range[0], y_range[0]
    Lx = float(x_range[1] - x_range[0])
    Ly = float(y_range[1] - y_range[0])
    if Lx <= 0 or Ly <= 0:
        raise ValueError('x_range and y_range must have positive width')

    modes = fourier_mode_list(M, N)
    n_atom = len(positions)

    x_rel = np.mod(positions[:, 0] - x0, Lx)
    y_rel = np.mod(positions[:, 1] - y0, Ly)
    z_at = positions[:, 2]

    design = _build_fourier_matrix(x_rel, y_rel, Lx, Ly, modes)
    theta, _residuals, rank, _singular_values = np.linalg.lstsq(design, z_at, rcond=rcond)
    n_columns = int(design.shape[1])
    if rank < n_columns:
        warnings.warn(
            'Fourier least-squares system is rank-deficient or underdetermined: '
            f'lstsq rank is {rank} for {n_columns} parameters and {n_atom} atoms '
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
    Height and :math:`\phi`-derivatives for one mode
    :math:`h(\phi)=A\cos\phi+B\sin\phi`.

    Used internally by :func:`_eval_fourier_surface` to accumulate
    ``fx``, ``fy``, ``fxx``, ``fyy``, ``fxy`` using the chain rule.
    For example, :math:`\partial_x h = (\partial_\phi h)\, \partial_x \phi`.

    Parameters
    ----------
    A, B : float
        Cosine and sine coefficients for this mode.
    cos_phase : ndarray
        ``cos(\phi)``, same shape as ``sin_phase``.
    sin_phase : ndarray
        ``sin(\phi)``.

    Notes
    -----
    This helper returns derivatives with respect to phase :math:`\phi`
    (``h_phi``, ``h_phiphi``), not spatial derivatives.
    Spatial derivatives (``fx``, ``fy``, ``fxx``, ``fyy``, ``fxy``) are
    accumulated in :func:`_eval_fourier_surface` via the chain rule
    using :math:`k_x` and :math:`k_y`.

    See Also
    --------
    :func:`_eval_fourier_surface`
        Consumes these intermediate terms to accumulate ``fx``, ``fy``,
        ``fxx``, ``fyy``, and ``fxy``.

    Returns
    -------
    h : ndarray
        Contribution :math:`A\cos\phi + B\sin\phi`.
    h_phi : ndarray
        Intermediate :math:`\partial h/\partial \phi` used for
        ``fx`` and ``fy`` accumulation.
    h_phiphi : ndarray
        Intermediate :math:`\partial^2 h/\partial \phi^2` used for
        ``fxx``, ``fyy``, ``fxy`` accumulation.
    """
    h = A * cos_phase + B * sin_phase
    h_phi = -A * sin_phase + B * cos_phase
    h_phiphi = -A * cos_phase - B * sin_phase
    return h, h_phi, h_phiphi


def _eval_fourier_surface(
    X_rel: np.ndarray,
    Y_rel: np.ndarray,
    Lx: float,
    Ly: float,
    A00: float,
    coeffs: dict[tuple[int, int], tuple[float, float]],
    *,
    derivatives: bool,
) -> np.ndarray | tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Evaluate the fitted Fourier sum on a relative-coordinate grid (steps 5-6).

    Step 5: reconstruct :math:`Z(x,y) = A_{00} + \sum_{(m,n)} A_{mn}\cos\phi_{mn}
    + B_{mn}\sin\phi_{mn}` with :math:`\phi_{mn} = k_x x + k_y y`.
    Step 6 (optional): accumulate first and second partial derivatives for Monge-gauge
    curvature using the analytic derivatives of the same sum.

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
        If ``False``, return only ``Z``. If ``True``, return ``Z`` and five derivative grids
        ``fx``, ``fy``, ``fxx``, ``fyy``, ``fxy``.

    Returns
    -------
    ndarray or tuple of ndarray
        If ``derivatives`` is ``False``, the height field ``Z`` with the same shape as
        ``X_rel``. If ``derivatives`` is ``True``, ``(Z, fx, fy, fxx, fyy, fxy)`` with each
        array matching the shape of ``X_rel``.
    """
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
    return Z


def _bin_centre_mesh(Lx: float, Ly: float, n_x_bins: int, n_y_bins: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Relative-coordinate mesh at bin centres for surface reconstruction (step 5).

    Bin centres lie at ``(i + 1/2) L / n_bins`` for each axis, compatible with the
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
    """
    xc = (np.arange(n_x_bins, dtype=np.float64) + 0.5) * Lx / n_x_bins
    yc = (np.arange(n_y_bins, dtype=np.float64) + 0.5) * Ly / n_y_bins
    X_rel, Y_rel = np.meshgrid(xc, yc, indexing='ij')
    return X_rel, Y_rel


# Steps 1 to 5; or 1 to 6 if derivatives are included
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
    Fourier height field on bin centres (module workflow steps 1-5).

    Runs steps 1-4 via :func:`_fourier_fit_from_atoms`, then step 5 via
    :func:`_bin_centre_mesh` and :func:`_eval_fourier_surface` with ``derivatives=False``.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atom coordinates in the same length unit as the box; column 2 is height.
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    n_x_bins, n_y_bins : int
        Output grid size at bin centres.
    M, N : int
        Maximum Fourier mode indices; see :func:`fourier_mode_list`.
    rcond : float or None, optional
        Passed to :func:`numpy.linalg.lstsq`.

    Returns
    -------
    Z : ndarray, shape (n_x_bins, n_y_bins)
        Fitted surface height.

    Warns
    -----
    UserWarning
        If the least-squares system is rank-deficient or underdetermined.
    """
    Lx, Ly, A00, coeffs = _fourier_fit_from_atoms(
        np.asarray(positions, dtype=np.float64),
        x_range,
        y_range,
        M,
        N,
        rcond=rcond,
    )
    X_rel, Y_rel = _bin_centre_mesh(Lx, Ly, n_x_bins, n_y_bins)
    return _eval_fourier_surface(X_rel, Y_rel, Lx, Ly, A00, coeffs, derivatives=False)


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
    Fourier height and analytic partial derivatives on bin centres (steps 1-6).

    Same pipeline as :func:`fourier_height_from_atoms`, but also performs step 6 by
    requesting analytic derivatives from :func:`_eval_fourier_surface` with ``derivatives=True``.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atom coordinates in the same length unit as the box; column 2 is height.
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    n_x_bins, n_y_bins : int
        Output grid size at bin centres.
    M, N : int
        Maximum Fourier mode indices; see :func:`fourier_mode_list`.
    rcond : float or None, optional
        Passed to :func:`numpy.linalg.lstsq`.

    Returns
    -------
    Z, fx, fy, fxx, fyy, fxy : ndarray
        Each array has shape ``(n_x_bins, n_y_bins)``.

    Warns
    -----
    UserWarning
        If the least-squares system is rank-deficient or underdetermined.
    """
    Lx, Ly, A00, coeffs = _fourier_fit_from_atoms(
        np.asarray(positions, dtype=np.float64),
        x_range,
        y_range,
        M,
        N,
        rcond=rcond,
    )
    X_rel, Y_rel = _bin_centre_mesh(Lx, Ly, n_x_bins, n_y_bins)
    return _eval_fourier_surface(X_rel, Y_rel, Lx, Ly, A00, coeffs, derivatives=True)
