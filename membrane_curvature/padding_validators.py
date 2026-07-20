"""

------------------
Padding Validators
------------------

.. versionadded:: 2.0.0

Helper functions to validate inputs when running
:class:`~membrane_curvature.base.MembraneCurvature` with ``padding=True``
on the binning surface method.

These validators check the simulation box and user-provided buffer width
before padded binning and finite differences run, so invalid inputs fail
fast with clear error messages.

.. warning::

   We provide the API documentation for the validators **for reference only.**

   **Validators are not intended to be called directly by the user**.
   They are called automatically when ``padding=True`` on
   :class:`~membrane_curvature.base.MembraneCurvature`.


Functions
---------

"""

import warnings

import numpy as np


def validate_orthorhombic(dimensions, atol=1e-3):
    """
    Check that the simulation cell is orthorhombic in the membrane plane.

    Parameters
    ----------
    dimensions : array_like
        MDAnalysis universe dimensions ``[lx, ly, lz, alpha, beta, gamma]``.
    atol : float, optional
        Absolute tolerance in degrees for angles to equal 90.

    Returns
    -------
    None
        Returns ``None`` when the box is accepted.

    Raises
    ------
    ValueError
        If ``dimensions`` is too short, or any box angle differs from 90°
        beyond ``atol``.

    Notes
    -----
    Image shifts :math:`\\pm L_x`, :math:`\\pm L_y` for padding are only valid
    for orthorhombic cells. Non-orthorhombic boxes need lattice-vector
    replicas, which are not included in the current implementation.
    """
    dims = np.asarray(dimensions, dtype=float)
    if dims.size < 6:
        raise ValueError('dimensions must include box lengths and angles')
    angles = dims[3:6]
    if not np.allclose(angles, 90.0, atol=atol):
        raise ValueError(f'edge padding requires an orthorhombic box (angles ~90°), got {angles}')


def validate_edge_pad_bins(edge_pad_bins, n_x_bins, n_y_bins):
    """
    Check ``edge_pad_bins`` when ``padding=True``.

    Parameters
    ----------
    edge_pad_bins : int
        Buffer width in bins on each side of the primary grid.
    n_x_bins : int
        Primary bin count along ``x``.
    n_y_bins : int
        Primary bin count along ``y``.

    Returns
    -------
    edge_pad_bins : int
        The validated buffer width.

    Raises
    ------
    ValueError
        If ``edge_pad_bins`` is not an integer ``>= 1``, or if it exceeds
        ``n_x_bins`` or ``n_y_bins``.

    Warns
    -----
    UserWarning
        If ``edge_pad_bins > 4`` (unlikely to change curvature; mainly cost).
    """
    if not isinstance(edge_pad_bins, (int, np.integer)) or edge_pad_bins < 1:
        raise ValueError(f'edge_pad_bins must be an int >= 1, got {edge_pad_bins!r}')
    if edge_pad_bins > 4:
        warnings.warn(
            'edge_pad_bins > 4 is unlikely to change curvature results '
            'and will slow down the calculation. Prefer values near 2.',
            category=UserWarning,
            stacklevel=3,
        )
    if edge_pad_bins > n_x_bins or edge_pad_bins > n_y_bins:
        raise ValueError(f'edge_pad_bins must be <= n_x_bins and n_y_bins, got {edge_pad_bins!r}')
    return int(edge_pad_bins)
