r"""

--------------------
Curvature
--------------------

In MembraneCurvature, we calculate curvature from a height field :math:`z = f(x,y)`
using the Monge gauge formulas, where partial derivatives are estimated from the
derived surface.

With the Monge gauge formulas, we map first derivatives :math:`(\partial_x, \partial_y)`, and second
and mixed derivatives :math:`(\partial_{xx}, \partial_{yy}, \partial_{xy})` to the definitions of
curvature.

Mean curvature :math:`H`, defined by:

.. math:: H =
    \frac{(1+\partial_x^2)\partial_{yy}+(1+\partial_y^2)\partial_{xx}-2\partial_x\partial_y\partial_{xy}}
    {2(1+\partial_x^2+\partial_y^2)^{3/2}},

and Gaussian curvature :math:`K`, defined by:

.. math:: K = \frac{\partial_{xx}\partial_{yy}-\partial_{xy}^2}
   {(1+\partial_x^2+\partial_y^2)^2}.


The Monge gauge formulas are implemented separately in the helpers
:func:`mean_curvature_monge` and :func:`gaussian_curvature_monge`, which accept
precomputed partial derivatives as arguments. This approach intends to separate the
curvature algebra from the derivative estimation.

To estimate partial derivatives, we use :func:`mean_curvature` and :func:`gaussian_curvature`
from the discrete height field using :func:`numpy.gradient`, then pass them to the
Monge gauge functions above. Optional spacing arguments (``dx``, ``dy``) are
forwarded to :func:`numpy.gradient` for physical-unit calculations.

Since the mean curvature calculates the arithmetic mean of two
principal curvatures, the default units of :math:`H` are Å\ :sup:`-1`.
On the other hand, Gaussian curvature calculates the geometric mean of the
two principal curvatures. Therefore, the default units of :math:`K` are Å\ :sup:`-2`.
In general, units of mean curvature are [length] :sup:`-1`,
and units of Gaussian curvature are [length] :sup:`-2`.

.. note::

    When spacing is provided to :func:`numpy.gradient` (for example, ``dx`` and
    ``dy`` from the simulation box and grid definition), curvature calculations are
    performed in physical units and are less sensitive to grid bin count than with
    implicit unit spacing. However, results are not exactly bin-independent because
    finite-difference discretization error increases for coarse grids.

.. warning::

    Numpy cannot calculate the gradient for arrays with inner array of
    `length==1` unless `axis=0` is specified. Therefore in the functions here included
    for mean and Gaussian curvature, shape of arrays must be at least (2,2).
    In general, to calculate a numerical gradients shape of arrays must be >=(`edge_order` +
    1).

For a periodic **Fourier** fit to atom heights, use
:func:`~membrane_curvature.fourier_surface.fourier_height_from_atoms` or
:func:`~membrane_curvature.fourier_surface.fourier_height_derivatives_from_atoms`
(surface only in :mod:`~membrane_curvature.fourier_surface`), then
:func:`fourier_curvature` for Monge :math:`H` and :math:`K`
in one call.

For binning with ``padding=True``, curvature is evaluated with
:func:`curvature_with_edge_pad`.


Functions
---------

"""

import numpy as np

from .fourier_surface import fourier_height_derivatives_from_atoms
from .padding import clip_padded_grid


def mean_curvature_monge(fx, fy, fxx, fyy, fxy):
    r"""
    Helper to calculate mean curvature :math:`H` from partial derivatives of :math:`z = f(x, y)`.

    Same expression as in :func:`mean_curvature` once first derivatives
    :math:`(\partial_x, \partial_y)` and second and mixed derivatives
    :math:`(\partial_{xx}, \partial_{yy}, \partial_{xy})` of :math:`z` are known,
    supplied as ``fx``, ``fy``, ``fxx``, ``fyy``, ``fxy`` (see the module equations
    above). This function does not compute gradients.

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
    Helper to calculate Gaussian curvature :math:`K` from partial derivatives of :math:`z = f(x, y)`.

    Same expression as in :func:`gaussian_curvature` once first derivatives
    :math:`(\partial_x, \partial_y)` and second and mixed derivatives
    :math:`(\partial_{xx}, \partial_{yy}, \partial_{xy})` of :math:`z` are known,
    supplied as ``fx``, ``fy``, ``fxx``, ``fyy``, ``fxy`` (see the module equations
    above). This function does not compute gradients.

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


def fourier_curvature(
    positions,
    x_range,
    y_range,
    n_x_bins,
    n_y_bins,
    M,
    N,
    rcond=None,
):
    """
    Bin-centre height and Monge mean / Gaussian curvature from a Fourier fit to atoms.

    Runs :func:`~membrane_curvature.fourier_surface.fourier_height_derivatives_from_atoms`
    once, then :func:`mean_curvature_monge` and :func:`gaussian_curvature_monge`.

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

    Returns
    -------
    z_surface, mean_H, gaussian_K : ndarray
        Each array has shape ``(n_x_bins, n_y_bins)``.

    Warns
    -----
    UserWarning
        If the Fourier least-squares system is rank-deficient or underdetermined.
    """
    Z, fx, fy, fxx, fyy, fxy = fourier_height_derivatives_from_atoms(
        positions,
        x_range,
        y_range,
        n_x_bins,
        n_y_bins,
        M,
        N,
        rcond=rcond,
    )
    H = mean_curvature_monge(fx, fy, fxx, fyy, fxy)
    K = gaussian_curvature_monge(fx, fy, fxx, fyy, fxy)
    return Z, H, K


def gaussian_curvature(Z, *varargs):
    """
    Calculate Gaussian curvature from Z cloud points.

    Uses :func:`numpy.gradient` on `Z`, then :func:`gaussian_curvature_monge`.

    Parameters
    ----------
    Z: np.ndarray.
        Multidimensional array of shape (n,n).
    varargs : list of scalar or array, optional
        Spacing between f values. Default unitary spacing for all dimensions.
        See np.gradient docs for more information.

    Returns
    -------
    K : np.ndarray.
        The result of Gaussian curvature of Z. Returns multidimensional
        array object with values of Gaussian curvature of shape `(n, n)`.

    """

    Zx, Zy = np.gradient(Z, *varargs)
    Zxx, Zxy = np.gradient(Zx, *varargs)
    _, Zyy = np.gradient(Zy, *varargs)

    return gaussian_curvature_monge(Zx, Zy, Zxx, Zyy, Zxy)


def mean_curvature(Z, *varargs):
    """
    Calculates mean curvature from Z cloud points.

    Uses :func:`numpy.gradient` on `Z`, then :func:`mean_curvature_monge`.

    Parameters
    ----------
    Z: np.ndarray.
        Multidimensional array of shape (n,n).
    varargs : list of scalar or array, optional
        Spacing between f values. Default unitary spacing for all dimensions.
        See np.gradient docs for more information.

    Returns
    -------
    H : np.ndarray.
        The result of mean curvature of Z. Returns multidimensional
        array object with values of mean curvature of shape `(n, n)`.

    """

    (
        Zx,
        Zy,
    ) = np.gradient(Z, *varargs)
    Zxx, Zxy = np.gradient(Zx, *varargs)
    _, Zyy = np.gradient(Zy, *varargs)

    return mean_curvature_monge(Zx, Zy, Zxx, Zyy, Zxy)


def curvature_with_edge_pad(z_surface_padded, dx, dy, edge_pad_bins):
    """
    Evaluate Monge curvature on a padded height field and clip to the primary box.

    Parameters
    ----------
    z_surface_padded : ndarray
        Height field on the expanded grid, shape
        ``(n_x_bins + 2*edge_pad_bins, n_y_bins + 2*edge_pad_bins)``.
    dx : float
        Bin spacing along ``x`` (Å).
    dy : float
        Bin spacing along ``y`` (Å).
    edge_pad_bins : int
        Buffer width in bins.

    Returns
    -------
    z_surface : ndarray
        Clipped height, shape ``(n_x_bins, n_y_bins)``.
    mean : ndarray
        Clipped mean curvature, same shape.
    gaussian : ndarray
        Clipped Gaussian curvature, same shape.

    .. note::

        Function available only when running with ``padding=True``. Finite differences
        run on the padded array and then clipped to the primary box. Gradients are
        computed once and shared by mean and Gaussian curvature.
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
    """
    Evaluate Monge curvature on a primary box height field using a wrap-pad buffer.

    Parameters
    ----------
    z_surface : ndarray
        Primary height field, shape ``(n_x_bins, n_y_bins)``.
    dx : float
        Bin spacing along ``x`` (Å).
    dy : float
        Bin spacing along ``y`` (Å).
    edge_pad_bins : int
        Buffer width in bins on each side.

    Returns
    -------
    z_surface : ndarray
        Same primary height (clipped identity of the wrap-padded array).
    mean : ndarray
        Mean curvature on the primary grid.
    gaussian : ndarray
        Gaussian curvature on the primary grid.

    Notes
    -----
    Used for average maps after optional FFT filtering, when atom images are no
    longer available and the height field is already periodic on the primary
    grid. The buffer is built with :func:`numpy.pad` ``mode='wrap'``.
    """
    p = int(edge_pad_bins)
    if p < 1:
        raise ValueError(f'edge_pad_bins must be >= 1, got {edge_pad_bins!r}')
    z_padded = np.pad(np.asarray(z_surface, dtype=float), pad_width=p, mode='wrap')
    return curvature_with_edge_pad(z_padded, dx, dy, p)
