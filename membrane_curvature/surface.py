r"""
.. role:: raw-math(raw) :format: latex html

--------------------
Surface
--------------------

Calculation of curvature requires a surface of reference. In MembraneCurvature,
the surface of reference is defined by the `z` position of the `atoms` in `AtomGroup`.


Functions
---------


"""

import numpy as np
import warnings


def derive_surface(atoms, n_cells_x, n_cells_y, max_width_x, max_width_y):
    """
    Derive surface from `atom` positions in `AtomGroup`.

    Parameters
    ----------
    atoms: AtomGroup.
        AtomGroup of reference selection to define the surface
        of the membrane.
    n_cells_x: int.
        number of cells in the grid of size `max_width_x`, `x` axis.
    n_cells_y: int.
        number of cells in the grid of size `max_width_y`, `y` axis.
    max_width_x: float.
        Maximum width of simulation box in x axis. (Determined by simulation box dimensions)
    max_width_y: float.
        Maximum width of simulation box in y axis. (Determined by simulation box dimensions)

    Returns
    -------
    z_coordinates: numpy.ndarray
        Average z-coordinate values. Return Numpy array of floats of
        shape `(n_cells_x, n_cells_y)`.

    """
    coordinates = atoms.positions
    return get_z_surface(coordinates, n_x_bins=n_cells_x, n_y_bins=n_cells_y,
                         x_range=(0, max_width_x), y_range=(0, max_width_y))


def get_z_surface(coordinates, n_x_bins=10, n_y_bins=10, x_range=(0, 100), y_range=(0, 100)):
    """
    Derive surface from distribution of z coordinates in grid.

    Parameters
    ----------
    coordinates : numpy.ndarray 
        Coordinates of AtomGroup. Numpy array of shape=(n_atoms, 3).
    n_x_bins : int.
        Number of bins in grid in the `x` dimension. 
    n_y_bins : int.
        Number of bins in grid in the `y` dimension. 
    x_range : tuple of (float, float)
        Range of coordinates (min, max) in the `x` dimension with shape=(2,).
    y_range : tuple of (float, float)
        Range of coordinates (min, max) in the `y` dimension with shape=(2,). 

    Returns
    -------
    z_surface: np.ndarray
        Surface derived from set of coordinates in grid of `x_range, y_range` dimensions.
        Returns Numpy array of floats of shape (`n_x_bins`, `n_y_bins`)

    """

    grid_z_coordinates = np.zeros((n_x_bins, n_y_bins))
    grid_norm_unit = np.zeros((n_x_bins, n_y_bins))

    x_factor = n_x_bins / (x_range[1] - x_range[0])
    y_factor = n_y_bins / (y_range[1] - y_range[0])

    x_coords, y_coords, z_coords = coordinates.T

    cell_x_floor = np.floor(x_coords * x_factor).astype(int)
    cell_y_floor = np.floor(y_coords * y_factor).astype(int)

    for l, m, z in zip(cell_x_floor, cell_y_floor, z_coords):

        try:
            if l < 0 or m < 0:
                msg = ("Atom outside grid boundaries. Skipping atom.")
                warnings.warn(msg)
                continue

            grid_z_coordinates[l, m] += z
            grid_norm_unit[l, m] += 1

        except IndexError:
            msg = ("Atom outside grid boundaries. Skipping atom.")
            warnings.warn(msg)

    z_surface = normalized_grid(grid_z_coordinates, grid_norm_unit)

    return z_surface


def normalized_grid(grid_z_coordinates, grid_norm_unit):
    """
    Calculates average `z` coordinates in unit cell.

    Parameters
    ----------

    z_ref: np.array
        Empty array of `(l,m)`
    grid_z_coordinates: np.array
        Array of size `(l,m)` with `z` coordinates stored in unit cell.
    grid_norm_unit: np.array
        Array of size `(l,m)` with number of atoms in unit cell.

    Returns
    -------
    z_surface: np.ndarray
        Normalized `z` coordinates in grid.
        Returns Numpy array of floats of shape (`n_x_bins`, `n_y_bins`)

    """

    grid_norm_unit = np.where(grid_norm_unit > 0, grid_norm_unit, np.nan)
    z_normalized = grid_z_coordinates / grid_norm_unit

    return z_normalized
