r"""

----------------------------
Binning Nearest Surface
----------------------------

.. versionadded:: 2.0.0

With ``surface_method='binning_nearest'``, :class:`~membrane_curvature.base.MembraneCurvature`
derives a surface by assigning each grid point the ``z`` of the nearest lipid in
the ``xy`` plane following the minimum-image convention. Grid points sit at bin corners
(``left + i * bin_size``), following the nearest-lipid assignment in g_lomepro
[CAG2009]_.

This height field is then used to calculate mean and Gaussian curvature using
finite differences.

.. [CAG2009] Vytautas Gapsys, Bert L. de Groot, Rodolfo Briones,
    *Computational analysis of local membrane properties*,
    Journal of Computer-Aided Molecular Design (2013), doi: `doi:10.1007/s10822-013-9684-0`_.

.. _`10.1007/s10822-013-9684-0`: https://doi.org/10.1007/s10822-013-9684-0

.. note::

    The binning nearest method supports a brick-wall filter to remove
    high-frequency noise from the height field. The filter is applied to the
    height field before calculating curvature. For more details, check the API
    documentation for the :ref:`fft_filtering` module.


.. warning::

    ``binning`` and ``fourier`` evaluate the surface on bin centres, so for the
    same bin counts and box extent their maps are offset by half a bin relative
    to ``binning_nearest``.

Functions
---------
"""

import MDAnalysis.lib.distances as mda_distances
import numpy as np
from MDAnalysis.exceptions import NoDataError


def lipid_center_positions(atomgroup):
    """
    One reference point per lipid: the center of the **selected** atoms in each
    residue.

    .. note::

        Only atoms in ``atomgroup`` contribute. Same selection as the surface method,
        not the full residue. For a one-atom reference such as ``name PO4`` the point is
        that bead's position. This is the same reference used by the ``binning`` and
        ``fourier`` surface methods. When several atoms per residue are selected they
        are reduced to their mass weighted center. Geometric center when masses are
        unavailable or all zero.

    Parameters
    ----------
    atomgroup : MDAnalysis.core.groups.AtomGroup
        AtomGroup of reference selection (one leaflet) that defines the surface.

    Returns
    -------
    centers : numpy.ndarray
        One ``(x, y, z)`` reference point per residue, shape ``(n_residues, 3)``.
    """
    try:
        masses = atomgroup.masses
    except (NoDataError, AttributeError):
        masses = None
    if masses is not None and np.any(masses != 0):
        return atomgroup.center_of_mass(compound='residues')
    return atomgroup.center_of_geometry(compound='residues')


def _grid_extent_from_box(box, grid_origin, lipid_centers):
    """Return (left_x, left_y, right_x, right_y) for the bin grid."""
    left_x, left_y = 0.0, 0.0
    right_x, right_y = float(box[0]), float(box[1])
    if grid_origin == 'lipid_bbox':
        left_x = float(np.min(lipid_centers[:, 0]))
        left_y = float(np.min(lipid_centers[:, 1]))
        right_x = float(np.max(lipid_centers[:, 0]))
        right_y = float(np.max(lipid_centers[:, 1]))
    return left_x, left_y, right_x, right_y


def _bin_corner_grid(left_x, left_y, bin_sizex, bin_sizey, n_x_bins, n_y_bins):
    """Bin reference points at lower-left corners."""
    x_corners = left_x + np.arange(n_x_bins, dtype=float) * bin_sizex
    y_corners = left_y + np.arange(n_y_bins, dtype=float) * bin_sizey
    x_grid, y_grid = np.meshgrid(x_corners, y_corners, indexing='ij')
    refs_xy = np.column_stack([x_grid.ravel(), y_grid.ravel(), np.zeros(x_grid.size)])
    return refs_xy, x_grid.shape


def _nearest_grid_axes(box, grid_origin, lipid_centers, n_x_bins, n_y_bins, x_range, y_range):
    """Return grid origin and bin widths for the nearest-neighbor surface."""
    left_x, left_y, right_x, right_y = _grid_extent_from_box(box, grid_origin, lipid_centers)
    if x_range is not None:
        left_x, right_x = x_range
    if y_range is not None:
        left_y, right_y = y_range
    bin_sizex = (right_x - left_x) / n_x_bins
    bin_sizey = (right_y - left_y) / n_y_bins
    return left_x, left_y, bin_sizex, bin_sizey


def _z_surface_from_nearest(refs_xy, lipid_centers, box, grid_shape):
    """Assign each grid point the ``z`` of the nearest lipid in the ``xy`` plane."""
    lipid_xy = lipid_centers.copy()
    lipid_xy[:, 2] = 0.0
    distances = mda_distances.distance_array(refs_xy, lipid_xy, box=box)
    nearest_index = np.argmin(distances, axis=1)
    return lipid_centers[nearest_index, 2].reshape(grid_shape)


def get_z_surface_nearest(
    atomgroup,
    n_x_bins=10,
    n_y_bins=10,
    box=None,
    grid_origin='box',
    x_range=None,
    y_range=None,
):
    """
    Derive surface by assigning each grid corner the ``z`` of the nearest lipid in ``xy``.

    The atomgroup should contain reference atoms for a single leaflet. The height
    field is sampled at bin corners (``left + i * bin_size``), not bin centres.
    Minimum-image distances use :func:`~MDAnalysis.lib.distances.distance_array`.

    Parameters
    ----------
    atomgroup : MDAnalysis.core.groups.AtomGroup
        Reference atoms for one leaflet (typically one atom per lipid).
    n_x_bins : int
        Number of bins along ``x``.
    n_y_bins : int
        Number of bins along ``y``.
    box : array-like, optional
        Unit cell in MDAnalysis format ``[lx, ly, lz, alpha, beta, gamma]``.
        Required for periodic minimum-image distances.
    grid_origin : {'box', 'lipid_bbox'}
        How to set the xy grid extent. ``box`` uses ``[0, lx]`` and ``[0, ly]``.
    x_range : tuple of (float, float), optional
        Override grid extent along ``x`` as ``(left_x, right_x)``.
    y_range : tuple of (float, float), optional
        Override grid extent along ``y`` as ``(left_y, right_y)``.

    Returns
    -------
    z_surface : numpy.ndarray
        Height field of shape ``(n_x_bins, n_y_bins)``.
    grid_axes : dict
        Per-frame grid metadata: ``left_x``, ``left_y``, ``dx``, ``dy``.
    """
    if box is None:
        raise ValueError('box is required for minimum-image distances')

    lipid_centers = lipid_center_positions(atomgroup)
    left_x, left_y, bin_sizex, bin_sizey = _nearest_grid_axes(
        box, grid_origin, lipid_centers, n_x_bins, n_y_bins, x_range, y_range
    )
    refs_xy, grid_shape = _bin_corner_grid(left_x, left_y, bin_sizex, bin_sizey, n_x_bins, n_y_bins)
    z_surface = _z_surface_from_nearest(refs_xy, lipid_centers, box, grid_shape)

    grid_axes = {
        'left_x': left_x,
        'left_y': left_y,
        'dx': bin_sizex,
        'dy': bin_sizey,
    }
    return z_surface, grid_axes
