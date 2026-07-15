r"""
--------------------
Edge padding
--------------------

.. versionadded:: 2.0.0

.. warning::

    This feature is optional and only available for orthorhombic boxes (angles 90°)
    and when ``surface_method='binning'``.

Module to apply periodic edge padding. A periodic buffer around the simulation
box is created for the binning surface method. Mean and Gaussian curvature on
the padded height field are evaluated in :mod:`~membrane_curvature.curvature`,
then the border is clipped so results match the primary ``n_x_bins x n_y_bins``
grid.

.. note::

    Typically, a buffer of 2 to 3 bins is sufficient to fix finite difference
    artifacts at edges and corners from ``numpy.gradient`` for second derivatives
    in mean and Gaussian curvature.

"""

import numpy as np

from .binning_surface import normalized_grid


def _expand_inclusive_ranges(lo, hi):
    """Expand inclusive ``[lo, hi]`` ranges into flat item and value index arrays."""
    lo = np.asarray(lo, dtype=np.int64)
    hi = np.asarray(hi, dtype=np.int64)
    counts = np.maximum(hi - lo + 1, 0)
    total = int(counts.sum())
    if total == 0:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty

    item = np.repeat(np.arange(lo.size, dtype=np.int64), counts)
    starts = np.empty(lo.size, dtype=np.int64)
    starts[0] = 0
    if lo.size > 1:
        np.cumsum(counts[:-1], out=starts[1:])
    offsets = np.arange(total, dtype=np.int64) - np.repeat(starts, counts)
    return item, lo[item] + offsets


def tile_xy_buffer(positions, box_length_x, box_length_y, pad_x, pad_y):
    r"""
    Build atom coordinates that populate the primary box plus a buffer of width $\Delta$.

    Parameters
    ----------
    positions : ndarray
        Primary box coordinates, shape ``(n_atoms, 3)``.
    box_length_x : float
        Box length along ``x`` (Å).
    box_length_y : float
        Box length along ``y`` (Å).
    pad_x : float
        Buffer width $\Delta_x$ in Å beyond each ``x`` edge.
    pad_y : float
        Buffer width $\Delta_y$ in Å beyond each ``y`` edge.

    Returns
    -------
    tiled : ndarray
        Coordinates of every periodic image
        ``(x + i L_x, y + j L_y, z)`` that falls inside the half-open domain
        ``[-pad_x, box_length_x + pad_x) x [-pad_y, box_length_y + pad_y)``.
        Shape ``(n_tiled, 3)``.
    """
    pos = np.asarray(positions, dtype=float)
    if pos.shape[0] == 0:
        return np.empty((0, 3), dtype=float)

    x_min = -float(pad_x)
    x_max = float(box_length_x) + float(pad_x)
    y_min = -float(pad_y)
    y_max = float(box_length_y) + float(pad_y)

    xs = pos[:, 0]
    ys = pos[:, 1]
    zs = pos[:, 2]

    # Half-open: x_min <= x + i*Lx < x_max
    i_lo = np.ceil((x_min - xs) / box_length_x).astype(np.int64)
    i_hi = np.ceil((x_max - xs) / box_length_x).astype(np.int64) - 1
    j_lo = np.ceil((y_min - ys) / box_length_y).astype(np.int64)
    j_hi = np.ceil((y_max - ys) / box_length_y).astype(np.int64) - 1

    atom_i, i_vals = _expand_inclusive_ranges(i_lo, i_hi)
    if atom_i.size == 0:
        return np.empty((0, 3), dtype=float)

    # Each (atom, i) image is paired with all j images of that atom.
    j_counts = np.maximum(j_hi[atom_i] - j_lo[atom_i] + 1, 0)
    total = int(j_counts.sum())
    if total == 0:
        return np.empty((0, 3), dtype=float)

    atom_ij = np.repeat(atom_i, j_counts)
    i_ij = np.repeat(i_vals, j_counts)

    starts = np.empty(atom_i.size, dtype=np.int64)
    starts[0] = 0
    if atom_i.size > 1:
        np.cumsum(j_counts[:-1], out=starts[1:])
    j_offsets = np.arange(total, dtype=np.int64) - np.repeat(starts, j_counts)
    j_vals = j_lo[atom_ij] + j_offsets

    out = np.empty((total, 3), dtype=float)
    out[:, 0] = xs[atom_ij] + i_ij * box_length_x
    out[:, 1] = ys[atom_ij] + j_vals * box_length_y
    out[:, 2] = zs[atom_ij]
    return out


def _mean_z_grid(coordinates, n_x_bins, n_y_bins, x_min, y_min, width_x, width_y):
    """Bin mean ``z`` on a half-open rectangular grid in the xy-plane."""
    if coordinates.size == 0:
        return np.full((n_x_bins, n_y_bins), np.nan)

    x_factor = n_x_bins / width_x
    y_factor = n_y_bins / width_y
    cell_x = np.floor((coordinates[:, 0] - x_min) * x_factor).astype(np.intp)
    cell_y = np.floor((coordinates[:, 1] - y_min) * y_factor).astype(np.intp)
    valid = (cell_x >= 0) & (cell_x < n_x_bins) & (cell_y >= 0) & (cell_y < n_y_bins)
    cell_x = cell_x[valid]
    cell_y = cell_y[valid]
    z_coords = coordinates[valid, 2]

    grid_z = np.zeros((n_x_bins, n_y_bins), dtype=float)
    grid_count = np.zeros((n_x_bins, n_y_bins), dtype=float)
    np.add.at(grid_z, (cell_x, cell_y), z_coords)
    np.add.at(grid_count, (cell_x, cell_y), 1.0)
    return normalized_grid(grid_z, grid_count)


def padded_grid_spec(n_x_bins, n_y_bins, x_range, y_range, edge_pad_bins):
    """
    Expand the binning grid by ``edge_pad_bins`` on each side.

    Parameters
    ----------
    n_x_bins : int
        Number of bins along ``x`` on the primary grid.
    n_y_bins : int
        Number of bins along ``y`` on the primary grid.
    x_range : tuple of float
        Primary ``(xmin, xmax)``, usually ``(0, Lx)``.
    y_range : tuple of float
        Primary ``(ymin, ymax)``, usually ``(0, Ly)``.
    edge_pad_bins : int
        Number of bins of buffer on each side. Must be ``>= 1``.

    Returns
    -------
    spec : dict
        Mapping with keys:

        - ``n_x_bins_pad``, ``n_y_bins_pad`` : int
            Expanded bin counts, ``n + 2 * edge_pad_bins``.
        - ``x_range_pad``, ``y_range_pad`` : tuple of float
            Expanded coordinate ranges
            ``(xmin - Δx, xmax + Δx)`` and ``(ymin - Δy, ymax + Δy)``.
        - ``dx``, ``dy`` : float
            Bin widths (unchanged by padding).
        - ``pad_x``, ``pad_y`` : float
            Physical buffer widths ``edge_pad_bins * dx`` and
            ``edge_pad_bins * dy``.
        - ``crop`` : tuple of int
            ``(edge_pad_bins, edge_pad_bins)`` used to slice the primary block
            from padded arrays as
            ``arr[p:-p, p:-p]`` with ``p = edge_pad_bins``.
    """
    dx = (x_range[1] - x_range[0]) / n_x_bins
    dy = (y_range[1] - y_range[0]) / n_y_bins
    pad_x = edge_pad_bins * dx
    pad_y = edge_pad_bins * dy

    return {
        'n_x_bins_pad': n_x_bins + 2 * edge_pad_bins,
        'n_y_bins_pad': n_y_bins + 2 * edge_pad_bins,
        'x_range_pad': (x_range[0] - pad_x, x_range[1] + pad_x),
        'y_range_pad': (y_range[0] - pad_y, y_range[1] + pad_y),
        'dx': dx,
        'dy': dy,
        'pad_x': pad_x,
        'pad_y': pad_y,
        'crop': (edge_pad_bins, edge_pad_bins),
    }


def clip_padded_grid(array, edge_pad_bins):
    """
    Remove the buffer border and return the primary box region.

    Parameters
    ----------
    array : ndarray
        Padded field of shape
        ``(n_x_bins + 2*edge_pad_bins, n_y_bins + 2*edge_pad_bins)``.
    edge_pad_bins : int
        Buffer width in bins used when the array was built.

    Returns
    -------
    clipped : ndarray
        Primary block of shape ``(n_x_bins, n_y_bins)``.
    """
    p = int(edge_pad_bins)
    return array[p:-p, p:-p].copy()


def get_z_surface_with_edge_pad(
    coordinates,
    n_x_bins,
    n_y_bins,
    x_range,
    y_range,
    edge_pad_bins,
    box_length_x=None,
    box_length_y=None,
):
    """
    Derive a primary box height field using a periodic buffer during binning.

    Parameters
    ----------
    coordinates : ndarray
        Atom coordinates, shape ``(n_atoms, 3)``. Need not be wrapped into the
        primary cell. Periodic images are selected by :func:`tile_xy_buffer`.
    n_x_bins : int
        Primary bin count along ``x``.
    n_y_bins : int
        Primary bin count along ``y``.
    x_range : tuple of float
        Primary ``x`` range.
    y_range : tuple of float
        Primary ``y`` range.
    edge_pad_bins : int
        Buffer width in bins on each side.
    box_length_x : float, optional
        Box length along ``x``. Defaults to ``x_range[1] - x_range[0]``.
    box_length_y : float, optional
        Box length along ``y``. Defaults to ``y_range[1] - y_range[0]``.

    Returns
    -------
    z_surface : ndarray
        Mean ``z`` per primary bin, shape ``(n_x_bins, n_y_bins)``. Empty bins
        remain ``nan``, same convention as
        :func:`~membrane_curvature.binning_surface.get_z_surface`.

    """
    z_padded = get_z_surface_padded(
        coordinates,
        n_x_bins,
        n_y_bins,
        x_range,
        y_range,
        edge_pad_bins,
        box_length_x=box_length_x,
        box_length_y=box_length_y,
    )
    return clip_padded_grid(z_padded, edge_pad_bins)


def get_z_surface_padded(
    coordinates,
    n_x_bins,
    n_y_bins,
    x_range,
    y_range,
    edge_pad_bins,
    box_length_x=None,
    box_length_y=None,
):
    r"""
    Bin height on the expanded ``box + $\Delta$`` grid (no clipping).

    Positions are tiled into the buffer, then binned with a vectorized
    mean `z` accumulator.
    """
    spec = padded_grid_spec(n_x_bins, n_y_bins, x_range, y_range, edge_pad_bins)
    lx = box_length_x if box_length_x is not None else (x_range[1] - x_range[0])
    ly = box_length_y if box_length_y is not None else (y_range[1] - y_range[0])

    tiled = tile_xy_buffer(coordinates, lx, ly, spec['pad_x'], spec['pad_y'])
    x0, x1 = spec['x_range_pad']
    y0, y1 = spec['y_range_pad']
    return _mean_z_grid(
        tiled,
        spec['n_x_bins_pad'],
        spec['n_y_bins_pad'],
        x0,
        y0,
        x1 - x0,
        y1 - y0,
    )
