r"""
.. _fft_filtering:

-------------
FFT filtering
-------------

.. versionadded:: 2.0.0

.. important::

    The brick-wall FFT filter is available only when :class:`~membrane_curvature.base.MembraneCurvature`
    runs with ``surface_method='binning'``. Filtering defaults to ``fft_filter=None``.
    Pass ``fft_filter='auto'`` to enable automatic filtering with low pass set as ``(0, 0.5 * q_Nyq)``.

    **Note that per-frame arrays are not FFT-filtered!**

    For average maps, :class:`~membrane_curvature.base.MembraneCurvature` applies the brick-wall
    filter once to the averaged ``z_surface`` and computes mean and Gaussian curvature on that filtered average.


This module implements a brick-wall filter on a binned height field in the reciprocal space.
This implementation uses :func:`numpy.fft.fft2` with physical bin widths :math:`\Delta x` and
:math:`\Delta y` in Å along the x and y dimensions, respectively.

Usage
------

Users can control filtering in three ways via the ``fft_filter`` argument:

- **Disabled (default for binning):** ``fft_filter=None``

  No FFT filtering is applied to the time-averaged surface.

- **Automatic:** ``fft_filter='auto'``

  Low-pass with ``q_high = 0.5 * q_Nyq``, resolved at runtime via
  :func:`resolve_fft_filter`.

- **Manual:** ``fft_filter={'q': (q_low, q_high)}``

  Users explicitly specify the ``(q_low, q_high)`` pair in rad/Å.

.. warning::

    **We recommend** ``fft_filter='auto'`` **when enabling the brick-wall filter**,
    unless there is a specific reason to choose bounds manually.
    With ``fft_filter='auto'``, bounds are ``(0, 0.5 * q_Nyq)`` from the bin grid
    ``dx`` and ``dy``.

    When using the manual mode, the following rules apply:

    - Keep ``q_low=0`` for a normal smoothing filter. If you set ``q_low`` above zero,
      the filter drops the average height of the membrane and behaves more like a
      high-pass (ripples only, no overall level).
    - With ``fft_filter='auto'``, the upper cutoff ``0.5 * q_Nyq`` is a safe default for
      most grids. **Note that it is not a single correct value for every system**. Change ``q_high`` only
      if you know you need sharper or smoother results.

    Filter cutoffs ``q_low`` and ``q_high`` are given in rad/Å. They describe how
    rippled the membrane height can be along the grid: small :math:`|q|` represents a broad,
    gentle undulation, while large :math:`|q|` represents fine, short-wavelength detail.
    The finest detail the bin grid can represent is set by the Nyquist frequency
    :func:`nyquist_q`, here in this module referred to as ``q_Nyq``, and which depends on
    ``dx`` and ``dy``.

The filter is applied in two steps. First, :func:`resolve_fft_filter` picks the pass band.
Then :func:`apply_fft_filter` FFTs the height, zeros modes outside the band, and inverse-FFTs back
to the surface in real space.

Limitations
-----------

- :func:`numpy.fft.fft2` treats the grid as **periodic**. Prefer ``wrap=True`` for
  binning on trajectories with periodic boundaries, or accept edge artifacts when the
  patch is non-periodic.
- The filter is a **brick-wall** mask in :math:`|q|`, which can cause ringing (Gibbs
  oscillations). Curvature uses second derivatives of the height field, so ringing can
  be amplified. ``fft_filter='auto'`` uses a conservative cutoff. Tune manually only
  if you know what you are doing.
- **Empty bins** (``NaN`` after binning): before filtering, empty bins are filled with
  the average height of occupied bins (not zero or neighbour interpolation). Large gaps
  can therefore add artifacts near empty regions in the filtered surface. Finer binning
  or fewer empty bins usually helps. See :func:`_height_for_fft`.
- **Non-square bins** (``dx != dy``): the pass band uses an **isotropic** cutoff on
  :math:`|q| = \sqrt{q_x^2 + q_y^2}`, while :func:`nyquist_q` uses
  :math:`\min(\pi/\Delta x, \pi/\Delta y)`. Modes that are Nyquist-resolvable along the
  finer-spaced axis but have :math:`|q|` above the coarse-axis limit are removed by the
  filter even though the grid could represent them along that axis.


Functions
----------
"""

import warnings

import numpy as np


def _validate_bin_widths(dx, dy):
    """
    Validate bin widths before calculating the Nyquist wavevector magnitude.

    Parameters
    ----------
    dx: float
        Bin width along ``x`` (Å).
    dy: float
        Bin width along ``y`` (Å).

    Raises
    ------
    ValueError
        If ``dx`` or ``dy`` is not strictly positive.
    """
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError(f'Bin widths must be positive nonzero. Got dx={dx}, dy={dy}')


def nyquist_q(dx, dy):
    r"""
    Return a conservative Nyquist wavevector magnitude for a bin grid after validation.
    The Nyquist wavevector magnitude is calculated as
    :math:`q_{\mathrm{Nyq}} = \min(\pi/\Delta x, \pi/\Delta y)` in :math:`\mathrm{rad}/\mathrm{\AA}`.

    Parameters
    ----------
    dx : float
        Bin width along ``x`` (Å).
    dy : float
        Bin width along ``y`` (Å).

    Returns
    -------
    q_nyq : float
        Nyquist wavevector magnitude in rad/Å.

    Raises
    ------
    ValueError
        If ``dx`` or ``dy`` is not strictly positive.
    """
    _validate_bin_widths(dx, dy)
    q_nyq = min(np.pi / dx, np.pi / dy)
    return q_nyq


def resolve_fft_filter(fft_filter, dx, dy):
    """
    Resolve ``fft_filter`` into ``(q_low, q_high)`` for binning analyses.

    Parameters
    ----------
    fft_filter : str or dict
        ``'auto'`` applies a low-pass with ``q_high = 0.5 * q_Nyq``.
        Otherwise ``{'q': (q_low, q_high)}`` in rad/Å.
    dx : float
        Bin width along ``x`` (Å).
    dy : float
        Bin width along ``y`` (Å).

    Returns
    -------
    q_bounds : tuple of (float, float)
        Resolved ``(q_low, q_high)`` in rad/Å.

    Raises
    ------
    ValueError
        If ``fft_filter`` is not ``'auto'`` or a dict with key ``'q'``.
        Invalid ``(q_low, q_high)`` values and bounds above Nyquist are raised or
        warned from :func:`_validate_q_pair`. Non-positive ``dx`` / ``dy`` are
        raised from :func:`nyquist_q` when bounds are resolved.

    Warns
    -----
    UserWarning
        See :func:`_validate_q_pair`.
    """
    if fft_filter == 'auto':
        q_bounds = _auto_q_bounds(dx, dy)
    elif isinstance(fft_filter, dict):
        if set(fft_filter.keys()) != {'q'}:
            raise ValueError('fft_filter must be None, "auto", or a dict of the form {"q": (q_low, q_high)}')
        q_bounds = _validate_q_pair(fft_filter['q'], dx, dy)
    else:
        raise ValueError('fft_filter must be None, "auto", or a dict of the form {"q": (q_low, q_high)}')
    return q_bounds


def _auto_q_bounds(dx, dy):
    """
    Calculates the default tuple ``(q_low, q_high)`` for ``fft_filter='auto'``.
    Keeps modes with :math:`|q| \leq 0.5\, q_{\mathrm{Nyq}}` (low-pass).

    Returns
    -------
    tuple_low_high : tuple of (float, float)
        ``(0.0, 0.5 * nyquist_q(dx, dy))``.
    """
    tuple_low_high = (0.0, 0.5 * nyquist_q(dx, dy))
    return tuple_low_high


def _validate_q_pair(bounds, dx, dy):
    """
    Parse and validate the bounds ``(q_low, q_high)`` for ``fft_filter``.
    If valid, the pass band is defined by ``q_low <= |q| <= q_high`` in rad/Å.

    Parameters
    ----------
    bounds : tuple or list
        Pair ``(q_low, q_high)`` in rad/Å.
    dx : float
        Bin width along ``x`` (Å).
    dy : float
        Bin width along ``y`` (Å).

    Returns
    -------
    (q_low, q_high) : tuple of (float, float)
        Validated ``(q_low, q_high)``.

    Raises
    ------
    ValueError
        If ``bounds`` is invalid, ``q_low == q_high``, or ``q_low`` exceeds Nyquist.

    Warns
    -----
    UserWarning
        If ``q_high`` is above what the bin grid can represent (Nyquist), or if
        ``q_low > 0`` (the filter no longer keeps the average height).
        Use ``q_low=0`` for a usual low-pass).
    """
    if not isinstance(bounds, (tuple, list)) or len(bounds) != 2:
        raise ValueError('fft_filter["q"] must be a pair (q_low, q_high)')
    q_low, q_high = bounds[0], bounds[1]
    if q_low < 0.0 or q_high < 0.0:
        raise ValueError('fft_filter bounds must be non-negative')
    if q_low > q_high:
        raise ValueError(f'fft_filter q_low must be <= q_high, got ({q_low}, {q_high})')
    if q_low == q_high:
        raise ValueError(
            f'fft_filter q_low and q_high must differ. Got ({q_low}, {q_high}). '
            'An equal band removes all but modes at a single |q| (often none on the grid), '
            'which zeros the filtered height. q_low must be <= q_high.'
        )
    if q_low > 0.0:
        warnings.warn(
            f'fft_filter q_low={q_low} removes the mean height. Use q_low=0.0 for a standard low-pass.',
            stacklevel=3,
        )

    q_nyq = nyquist_q(dx, dy)
    if q_low > q_nyq:
        raise ValueError(
            f'fft_filter q_low={q_low} exceeds the Nyquist limit q_Nyq={q_nyq:.6g} rad/Å '
            'for this bin grid. No wavelengths on the grid fall in the pass band. '
            'Use a smaller q_low (typically 0.0) or fft_filter="auto".'
        )
    if q_high > q_nyq:
        warnings.warn(
            f'fft_filter q_high={q_high} exceeds Nyquist q_Nyq={q_nyq:.6g} rad/Å. '
            'Modes above Nyquist are not resolved by the grid.',
            stacklevel=3,
        )

    return q_low, q_high


def apply_fft_filter(height, dx, dy, q_bounds):
    """
    Apply brick-wall FFT filtering to a binned height field.

    Parameters
    ----------
    height : numpy.ndarray
        Binned height field, shape ``(n_x_bins, n_y_bins)``.
    dx : float
        Bin width along ``x`` (Å).
    dy : float
        Bin width along ``y`` (Å).
    q_bounds : tuple of (float, float)
        ``(q_low, q_high)`` in rad/Å from :func:`resolve_fft_filter`.

    Returns
    -------
    filtered_height : numpy.ndarray
        Filtered height field.
    """
    q_low, q_high = q_bounds
    n_x_bins, n_y_bins = height.shape
    low_bound_sq = q_low**2
    high_bound_sq = q_high**2

    # handle empty bins
    nan_mask = np.isnan(height)
    masked_height = _height_for_fft(height, nan_mask)

    # apply the filter
    spectrum = np.fft.fft2(masked_height)
    radius_sq = _squared_radius_grid(n_x_bins, n_y_bins, dx, dy)
    keep_modes = (radius_sq >= low_bound_sq) & (radius_sq <= high_bound_sq)
    spectrum = np.where(keep_modes, spectrum, 0.0)

    # inverse FFT and restore NaN bins
    filtered_height = np.fft.ifft2(spectrum).real
    filtered_height[nan_mask] = np.nan
    return filtered_height


def _squared_radius_grid(n_x_bins, n_y_bins, dx, dy):
    r"""
    Build a per-bin grid of squared angular wavevector magnitude.

    Each entry matches the corresponding bin in the output of
    :func:`numpy.fft.fft2` (via :func:`numpy.fft.fftfreq`).

    Parameters
    ----------
    n_x_bins : int
        Number of bins along ``x``.
    n_y_bins : int
        Number of bins along ``y``.
    dx : float
        Bin width along ``x`` (Å).
    dy : float
        Bin width along ``y`` (Å).

    Returns
    -------
    grid_squared_radius : numpy.ndarray
        Array of shape ``(n_x_bins, n_y_bins)`` with values :math:`|q|^2` in
        (rad/Å)\ :sup:`2`.
    """
    qx = np.fft.fftfreq(n_x_bins, d=dx) * (2.0 * np.pi)
    qy = np.fft.fftfreq(n_y_bins, d=dy) * (2.0 * np.pi)
    qx_grid, qy_grid = np.meshgrid(qx, qy, indexing='ij')
    grid_squared_radius = qx_grid**2 + qy_grid**2
    return grid_squared_radius


def _height_for_fft(height, nan_mask):
    """
    Prepare a height field for :func:`numpy.fft.fft2` when some bins are empty.

    Parameters
    ----------
    height : numpy.ndarray
        Binned height field, shape ``(n_x_bins, n_y_bins)``.
    nan_mask : numpy.ndarray
        Boolean mask of shape ``(n_x_bins, n_y_bins)`` with ``True`` for empty bins.

    Returns
    -------
    filled_height : numpy.ndarray
        Height field with empty bins filled to the average height of occupied bins.

    Notes
    -----

    Empty bins are stored as ``NaN``. Before the forward FFT, each empty bin is
    replaced with the mean height of occupied bins. The fill value is uniform
    across empty bins and it is neither zero nor obtained by spatial interpolation.

    This is a simple placeholder to run the filter. It can make the surface look
    too smooth or wavy near large gaps, especially when many bins are empty. After
    the inverse FFT, empty bins are set back to ``NaN``. If all bins are empty, the
    fill value is ``0.0``.
    """
    if not nan_mask.any():
        return height
    fill_value = np.nanmean(height)
    if not np.isfinite(fill_value):
        fill_value = 0.0
    filled_height = height.copy()
    filled_height[nan_mask] = fill_value
    return filled_height
