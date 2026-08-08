r"""

--------------------
Monge gauge formulas
--------------------

In MembraneCurvature, we calculate curvature from a derived surface
:math:`z = f(x,y)` using the Monge gauge formulas. Partial derivatives of the
surface are obtained numerically or analytically, depending on the surface
method.

With the Monge gauge formulas, we map first derivatives :math:`(\partial_x, \partial_y)`, and second
and mixed derivatives :math:`(\partial_{xx}, \partial_{yy}, \partial_{xy})` to the definitions of
curvature:

.. _mean_curvature_formula:

Mean curvature (:math:`H`)
--------------------------

Defined by:

.. math:: H =
    \frac{(1+\partial_x^2)\partial_{yy}+(1+\partial_y^2)\partial_{xx}-2\partial_x\partial_y\partial_{xy}}
    {2(1+\partial_x^2+\partial_y^2)^{3/2}}.

.. _gaussian_curvature_formula:

Gaussian curvature (:math:`K`)
------------------------------

Defined by:

.. math:: K = \frac{\partial_{xx}\partial_{yy}-\partial_{xy}^2}
   {(1+\partial_x^2+\partial_y^2)^2}.

In general, units of mean curvature are [length] :sup:`-1`,
and units of Gaussian curvature are [length] :sup:`-2`.
Mean curvature is the arithmetic mean of the two principal curvatures. The
default units of :math:`H` are Å\ :sup:`-1`.
Gaussian curvature is the product of the two principal curvatures. The
default units of :math:`K` are Å\ :sup:`-2`.

The Monge gauge formulas are implemented separately in the helpers
:func:`mean_curvature_monge` and :func:`gaussian_curvature_monge`, which accept
precomputed partial derivatives as arguments. This approach intends to separate the
curvature algebra from how the derivatives are obtained.

.. important::

    Depending on the :attr:`~membrane_curvature.base.MembraneCurvature.surface_method`
    selected, partial derivatives are obtained in two ways:

    - :ref:`numerical_derivatives` calculate curvature from the discrete surface using
      :func:`numpy.gradient`, then call the Monge gauge formulas.

    - :ref:`analytic_derivatives` calculate curvature from analytic derivatives of the
      fitted Fourier series.

For details on curvature calculation, check the :ref:`algorithm` page.
"""

import numpy as np

from .fourier_surface import (
    _bin_centre_mesh,
    _eval_fourier_surface,
    _fourier_fit_from_atoms,
    _unpack_coefficients,
)
from .padding import clip_padded_grid


def mean_curvature_monge(fx, fy, fxx, fyy, fxy):
    r"""
    Helper to calculate :ref:`mean_curvature_formula` from partial derivatives of :math:`z = f(x, y)`.

    Parameters
    ----------
    fx : np.ndarray.
        Values of :math:`\partial_x` on the grid.
    fy : np.ndarray.
        Values of :math:`\partial_y` on the grid.
    fxx : np.ndarray.
        Values of :math:`\partial_{xx}` on the grid.
    fyy : np.ndarray.
        Values of :math:`\partial_{yy}` on the grid.
    fxy : np.ndarray.
        Values of :math:`\partial_{xy}` on the grid.

    Returns
    -------
    H : np.ndarray.
        Mean curvature, same shape as the input arrays.

    """
    H = (1 + fx**2) * fyy + (1 + fy**2) * fxx - 2 * fx * fy * fxy
    H = H / (2 * (1 + fx**2 + fy**2) ** (1.5))
    return H


def gaussian_curvature_monge(fx, fy, fxx, fyy, fxy):
    r"""
    Helper to calculate :ref:`gaussian_curvature_formula` from partial derivatives of :math:`z = f(x, y)`.

    Parameters
    ----------
    fx : np.ndarray.
        Values of :math:`\partial_x` on the grid.
    fy : np.ndarray.
        Values of :math:`\partial_y` on the grid.
    fxx : np.ndarray.
        Values of :math:`\partial_{xx}` on the grid.
    fyy : np.ndarray.
        Values of :math:`\partial_{yy}` on the grid.
    fxy : np.ndarray.
        Values of :math:`\partial_{xy}` on the grid.

    Returns
    -------
    K : np.ndarray.
        Gaussian curvature, same shape as the input arrays.

    """
    return (fxx * fyy - (fxy**2)) / (1 + (fx**2) + (fy**2)) ** 2


def gaussian_curvature(Z, *varargs):
    """
    Calculate :ref:`gaussian_curvature_formula` from a surface array.

    Partial derivatives are estimated with :func:`numpy.gradient` and passed to
    :func:`gaussian_curvature_monge`.

    Parameters
    ----------
    Z : np.ndarray
        Surface array of shape ``(n_x_bins, n_y_bins)``.
    *varargs : list of scalar or array, optional
        Spacing between values. Default unitary spacing for all dimensions.
        See :func:`numpy.gradient` for more information.

    Returns
    -------
    K : np.ndarray
        Gaussian curvature, shape ``(n_x_bins, n_y_bins)``.

    """

    Zx, Zy = np.gradient(Z, *varargs)
    Zxx, Zxy = np.gradient(Zx, *varargs)
    _, Zyy = np.gradient(Zy, *varargs)

    return gaussian_curvature_monge(Zx, Zy, Zxx, Zyy, Zxy)


def mean_curvature(Z, *varargs):
    """
    Calculate :ref:`mean_curvature_formula` from a surface array.

    Partial derivatives are estimated with :func:`numpy.gradient` and passed to
    :func:`mean_curvature_monge`.

    Parameters
    ----------
    Z : np.ndarray
        Surface array of shape ``(n_x_bins, n_y_bins)``.
    *varargs : list of scalar or array, optional
        Spacing between values. Default unitary spacing for all dimensions.
        See :func:`numpy.gradient` for more information.

    Returns
    -------
    H : np.ndarray
        Mean curvature, shape ``(n_x_bins, n_y_bins)``.

    """

    (
        Zx,
        Zy,
    ) = np.gradient(Z, *varargs)
    Zxx, Zxy = np.gradient(Zx, *varargs)
    _, Zyy = np.gradient(Zy, *varargs)

    return mean_curvature_monge(Zx, Zy, Zxx, Zyy, Zxy)


def fourier_curvature_from_coefficients(
    A00,
    coeffs,
    x_range,
    y_range,
    n_x_bins,
    n_y_bins,
):
    r"""
    Calculate the surface and mean and Gaussian curvature from Fourier
    coefficients.

    Partial derivatives are obtained analytically from the fitted series and
    passed to :func:`mean_curvature_monge` and :func:`gaussian_curvature_monge`.

    Parameters
    ----------
    A00 : float
        Mean term :math:`A_{00}`.
    coeffs : dict
        Mapping ``(m, n) -> (A, B)`` cosine and sine coefficients per mode.
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    n_x_bins, n_y_bins : int
        Output grid size at bin centres.

    Returns
    -------
    z_surface, mean_H, gaussian_K : ndarray
        Each array has shape ``(n_x_bins, n_y_bins)``.
    """
    box_length_x = x_range[1] - x_range[0]
    box_length_y = y_range[1] - y_range[0]
    x_rel, y_rel = _bin_centre_mesh(box_length_x, box_length_y, n_x_bins, n_y_bins)
    z_surface, fx, fy, fxx, fyy, fxy = _eval_fourier_surface(
        x_rel,
        y_rel,
        box_length_x,
        box_length_y,
        A00,
        coeffs,
        derivatives=True,
    )
    mean_H = mean_curvature_monge(fx, fy, fxx, fyy, fxy)
    gaussian_K = gaussian_curvature_monge(fx, fy, fxx, fyy, fxy)
    return z_surface, mean_H, gaussian_K


def fourier_curvature_from_theta(
    theta,
    modes,
    x_range,
    y_range,
    n_x_bins,
    n_y_bins,
):
    r"""
    Calculate the surface and mean and Gaussian curvature from a Fourier
    coefficient vector :math:`\theta`.

    Splits :math:`\theta` with :func:`~membrane_curvature.fourier_surface._unpack_coefficients`
    and calls :func:`fourier_curvature_from_coefficients`.

    Parameters
    ----------
    theta : ndarray, shape (n_coefficients,)
        Least-squares coefficient vector with mean at index 0, then pairs
        :math:`(A_{mn}, B_{mn})` in the same order as ``modes``.
    modes : list of tuple of int
        Mode list from :func:`~membrane_curvature.fourier_surface.fourier_mode_list`.
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    n_x_bins, n_y_bins : int
        Output grid size at bin centres.

    Returns
    -------
    z_surface, mean_H, gaussian_K : ndarray
        Each array has shape ``(n_x_bins, n_y_bins)``.
    """
    A00, coeffs = _unpack_coefficients(theta, modes)
    return fourier_curvature_from_coefficients(
        A00,
        coeffs,
        x_range,
        y_range,
        n_x_bins,
        n_y_bins,
    )


def fourier_curvature(
    positions,
    x_range,
    y_range,
    n_x_bins,
    n_y_bins,
    M,
    N,
    rcond=None,
    return_theta=False,
):
    r"""
    Calculate the surface and mean and Gaussian curvature from a Fourier fit to
    atom positions.

    Fits coefficients with :func:`~membrane_curvature.fourier_surface._fourier_fit_from_atoms`
    and calls :func:`fourier_curvature_from_coefficients`.

    Parameters
    ----------
    positions : ndarray, shape (n_atoms, 3)
        Atom coordinates in the same length unit as the box (typically Å).
    x_range, y_range : tuple of float
        ``(min, max)`` domain extents.
    n_x_bins, n_y_bins : int
        Output grid size at bin centres.
    M, N : int
        Maximum Fourier mode indices; see
        :func:`~membrane_curvature.fourier_surface.fourier_mode_list`.
    rcond : float or None, optional
        Passed to singular-value truncation for the Fourier fit.
        See :func:`~membrane_curvature.fourier_surface._solve_design_least_squares_svd`.
    return_theta : bool, optional
        If ``True``, also return the coefficient vector :math:`\theta`.

    Returns
    -------
    z_surface, mean_H, gaussian_K : ndarray
        Each array has shape ``(n_x_bins, n_y_bins)``.
    theta : ndarray, optional
        Coefficient vector. Returned only when ``return_theta`` is ``True``.

    Warns
    -----
    UserWarning
        If the Fourier least-squares system is rank-deficient or underdetermined.
    """
    _box_length_x, _box_length_y, A00, coeffs, theta = _fourier_fit_from_atoms(
        positions,
        x_range,
        y_range,
        M,
        N,
        rcond=rcond,
    )
    z_surface, mean_H, gaussian_K = fourier_curvature_from_coefficients(
        A00,
        coeffs,
        x_range,
        y_range,
        n_x_bins,
        n_y_bins,
    )
    if return_theta:
        return z_surface, mean_H, gaussian_K, theta
    return z_surface, mean_H, gaussian_K


def curvature_with_edge_pad(z_surface_padded, dx, dy, edge_pad_bins):
    r"""
    Calculate mean and Gaussian curvature from a padded surface and clip the
    result to the primary box.

    Parameters
    ----------
    z_surface_padded : ndarray
        Surface on the expanded grid, shape
        ``(n_x_bins + 2*edge_pad_bins, n_y_bins + 2*edge_pad_bins)``.
    dx : float
        Bin spacing along :math:`x` (Å).
    dy : float
        Bin spacing along :math:`y` (Å).
    edge_pad_bins : int
        Buffer width in bins.

    Returns
    -------
    z_surface : ndarray
        Clipped surface, shape ``(n_x_bins, n_y_bins)``.
    mean : ndarray
        Clipped mean curvature, same shape.
    gaussian : ndarray
        Clipped Gaussian curvature, same shape.
    """
    zx, zy = np.gradient(z_surface_padded, dx, dy)
    zxx, zxy = np.gradient(zx, dx, dy)
    _, zyy = np.gradient(zy, dx, dy)
    mean_padded = mean_curvature_monge(zx, zy, zxx, zyy, zxy)
    gauss_padded = gaussian_curvature_monge(zx, zy, zxx, zyy, zxy)
    return (
        clip_padded_grid(z_surface_padded, edge_pad_bins),
        clip_padded_grid(mean_padded, edge_pad_bins),
        clip_padded_grid(gauss_padded, edge_pad_bins),
    )


def curvature_from_primary_with_edge_pad(z_surface, dx, dy, edge_pad_bins):
    r"""
    Calculate mean and Gaussian curvature from a primary surface with a periodic
    buffer.

    Adds a periodic buffer around the primary surface and calls
    :func:`curvature_with_edge_pad`.

    Parameters
    ----------
    z_surface : ndarray
        Surface in the primary box, shape ``(n_x_bins, n_y_bins)``.
    dx : float
        Bin spacing along :math:`x` (Å).
    dy : float
        Bin spacing along :math:`y` (Å).
    edge_pad_bins : int
        Buffer width in bins on each side.

    Returns
    -------
    z_surface : ndarray
        The original surface, returned after removing the extra padding added for
        calculation.
    mean : ndarray
        Mean curvature on the primary grid.
    gaussian : ndarray
        Gaussian curvature on the primary grid.
    """
    p = int(edge_pad_bins)
    if p < 1:
        raise ValueError(f'edge_pad_bins must be >= 1, got {edge_pad_bins!r}')
    z_padded = np.pad(np.asarray(z_surface, dtype=float), pad_width=p, mode='wrap')
    return curvature_with_edge_pad(z_padded, dx, dy, p)
