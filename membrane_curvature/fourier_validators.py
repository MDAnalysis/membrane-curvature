"""

.. versionadded:: 2.0.0

Input checks used before running a Fourier-surface calculation.

.. warning::

    These helpers check the atom coordinates, grid sizes, and related inputs
    that are passed into the Fourier surface API and raise a clear error if
    something is off. They are called automatically by the high-level entry
    points in :mod:`membrane_curvature.fourier_surface`
    (:func:`~membrane_curvature.fourier_surface.fourier_height_from_atoms`
    and :func:`~membrane_curvature.fourier_surface.fourier_height_derivatives_from_atoms`),
    so you do not normally have to call them yourself.

.. note::

    If an input is invalid, you will see a :class:`ValueError` message
    explaining what is wrong. For example ``positions must contain at least
    one row`` or ``n_x_bins and n_y_bins must be positive``. For coordinate
    input, the original NumPy message is kept attached to the error (chained
    via ``__cause__``) so the full traceback still shows where the problem started.

"""

import numpy as np


def _coerce_positions(positions) -> np.ndarray:
    """
    Convert atom coordinates into a NumPy array of floats.

    Internal helper used by the public Fourier entry points to accept
    list/tuple/array-like input and produce a contiguous
    :class:`numpy.ndarray` of ``float64`` values.

    Shape checks (2D, at least three columns, non-empty) are handled
    separately by :func:`validate_positions`.

    Parameters
    ----------
    positions : array-like
        Atom coordinates. Any data type that can be turned into a numeric
        NumPy array works: nested lists/tuples, a NumPy array, etc.

    Returns
    -------
    positions : numpy.ndarray, dtype float64
        The input as a float NumPy array.

    Raises
    ------
    ValueError
        If the input cannot be converted to a numeric array. Common causes
        are non-numeric entries (for example strings) or rows of different
        lengths (a "ragged" list of lists). The original NumPy error is
        chained via ``__cause__`` so the traceback shows what NumPy
        actually complained about.
    """
    try:
        return np.asarray(positions, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            'positions could not be converted to a float64 NumPy array '
            '(probably due to non-numeric values or ragged sequences)'
        ) from exc


def validate_positions(positions: np.ndarray) -> None:
    """
    Check that an atom coordinate array has the expected shape.

    A Fourier fit needs a 2D array of shape ``(n_atoms, 3)`` (or more
    columns; only the first three are used as :math:`x`, :math:`y`,
    :math:`z`). This helper raises a clear error if that is not the case.
    Typically called after coercing the input to a NumPy float array.

    Parameters
    ----------
    positions : numpy.ndarray
        Atom coordinates. Usually the output of :func:`_coerce_positions`.

    Raises
    ------
    ValueError
        Raised when the array does not look like a list of 3D points. The
        message states which expectation failed:

        - not 2D (for example, a flat list of numbers),
        - fewer than three columns (missing :math:`x`, :math:`y`, or
          :math:`z`),
        - no rows at all (the selection is empty).
    """
    if positions.ndim != 2:
        raise ValueError(f'positions must be a 2D array; got ndim={positions.ndim}')
    if positions.shape[1] < 3:
        raise ValueError(f'positions must have at least 3 columns (got shape {positions.shape})')
    if positions.shape[0] == 0:
        raise ValueError('positions must contain at least one row')


def validate_positive_domain_widths(
    x_range: tuple[float, float],
    y_range: tuple[float, float],
) -> None:
    r"""
    Check that the periodic domain has strictly positive width on each axis.

    The Fourier fit treats :math:`L_x = x_{\max} - x_{\min}` and
    :math:`L_y = y_{\max} - y_{\min}` as the period lengths. A zero or
    negative width would make wavevectors and the bin spacing wrong.
    This helper catches that mistake early.

    Parameters
    ----------
    x_range, y_range : tuple of (float, float)
        ``(min, max)`` extents of the periodic domain along :math:`x` and
        :math:`y`.

    Raises
    ------
    ValueError
        If ``x_range[1] - x_range[0]`` or ``y_range[1] - y_range[0]`` is
        not strictly positive (``max`` must be greater than ``min``).
    """
    Lx = float(x_range[1] - x_range[0])
    Ly = float(y_range[1] - y_range[0])
    if Lx <= 0 or Ly <= 0:
        raise ValueError(
            'x_range and y_range must have positive width; '
            f'got x_range={tuple(x_range)} (Lx={Lx}), '
            f'y_range={tuple(y_range)} (Ly={Ly})'
        )


def validate_positive_bin_counts(n_x_bins, n_y_bins) -> None:
    """
    Validate that grid bin counts are strictly positive integers.

    Catches invalid grid sizes before they propagate to operations that
    would otherwise produce an empty mesh or divide by zero (for example
    ``Lx / n_x_bins`` and ``Ly / n_y_bins`` in
    :func:`~membrane_curvature.fourier_surface._bin_centre_mesh`).

    Parameters
    ----------
    n_x_bins : int
        Number of bins along :math:`x`. Python :class:`int` and
        :class:`numpy.integer` are accepted; floats and :class:`bool` are
        rejected even though :class:`bool` is an :class:`int` subclass.
    n_y_bins : int
        Number of bins along :math:`y`. Same constraint as ``n_x_bins``.

    Raises
    ------
    ValueError
        Raised if either input is not a positive integer. Three cases:

        - non-integer type, that is, not ``int`` / :class:`numpy.integer`
          (for example ``1.5``, ``"10"``, or ``None``),
        - :class:`bool` (``True`` / ``False``), rejected explicitly even
          though ``bool`` is an :class:`int` subclass,
        - integer that is zero or negative.

        The message names the offending argument, its value, and its type
        so the mistake is easy to spot.
    """
    type_problems = []
    for name, value in (('n_x_bins', n_x_bins), ('n_y_bins', n_y_bins)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            type_problems.append(f'{name}={value!r} (type {type(value).__name__})')
    if type_problems:
        raise ValueError('n_x_bins and n_y_bins must be integers; got ' + ', '.join(type_problems))
    if n_x_bins <= 0 or n_y_bins <= 0:
        raise ValueError(f'n_x_bins and n_y_bins must be positive; got {n_x_bins=}, {n_y_bins=}')
