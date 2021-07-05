import numpy as np


def derive_surface(n_cells_x, n_cells_y, selection, max_width_x, max_width_y):
    """
    Derive surface from AtomGroup positions.

    Parameters
    ----------
    n_cells_x : int.
        number of cells in the grid of size `max_width_x`, `x` axis.
    n_cells_y : int.
        number of cells in the grid of size `max_width_y`, `y` axis.
    selection : AtomGroup.
        AtomGroup of reference selection to define the surface
        of the membrane.
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
    coordinates = selection.positions
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
        Range of indexes in grid in the `x` dimension with shape=(2,).
    y_range : tuple of (float, float)
        Range of indexes in grid in the `y` dimension with shape=(2,).


    Returns
    -------
    z_surface: numpy.ndarray 
        Surface derived from set of coordinates in grid of `x_range, y_range` dimensions.
        Returns Numpy array of floats of shape (`n_x_bins`, `n_y_bins`)

    """

    z_ref = np.zeros((n_x_bins, n_y_bins))
    grid_z_coordinates = np.zeros((n_x_bins, n_y_bins))
    grid_norm_unit = np.zeros((n_x_bins, n_y_bins))

    x_factor = n_x_bins / (x_range[1] - x_range[0])
    y_factor = n_y_bins / (y_range[1] - y_range[0])

    x_coords, y_coords, z_coords = coordinates.T

    cell_x_floor = np.floor(x_coords * x_factor).astype(int)
    cell_y_floor = np.floor(y_coords * y_factor).astype(int)

    for l, m, z in zip(cell_x_floor, cell_y_floor, z_coords):

        try:
            grid_z_coordinates[l, m] += z
            grid_norm_unit[l, m] += 1

        except IndexError:
            print("Atom outside grid boundaries. Skipping atom.")

    z_surface = avg_unit_cell(z_ref, grid_z_coordinates, grid_norm_unit)

    return z_surface


def avg_unit_cell(z_ref, grid_z_coordinates, grid_norm_unit):
    """
    Calculates average z coordinate in unit cell

    Parameters
    ----------

    z_ref: np.array
        Empty array of (l,m) 
    grid_z_coordinates: np.array
        Array of size (l,m) with z coordinates stored in unit cell.
    grid_norm_unit: np.array
        Array of size (l,m) with number of atoms in unit cell.

    Returns
    -------
    Returns nornalized set of z coordinates in grid.

    """

    grid_norm_unit = np.where(grid_norm_unit > 0, grid_norm_unit, np.nan)
    z_ref = grid_z_coordinates / grid_norm_unit

    return z_ref


def gaussian_curvature(Z):
    """
    Calculate gaussian curvature from Z cloud points.


    Parameters
    ----------
    Z : Numpy array, cloud of points.


    Returns
    -------
    K : 2d-array
        Returns 2-dimensional array object with values of mean curvature.

    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    K = (Zxx * Zyy - (Zxy ** 2)) / (1 + (Zx ** 2) + (Zy ** 2)) ** 2

    return K


def mean_curvature(Z):
    """
    Calculates mean curvature from Z cloud points.


    Parameters
    ----------
    Z : Numpy array, cloud of points.


    Returns
    -------
    H : 2d-array
        Returns 2-dimensional array object with values of gaussian curvature.

    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    H = (Zx**2 + 1) * Zyy - 2 * Zx * Zy * Zxy + (Zy**2 + 1) * Zxx
    H = -H / (2 * (Zx**2 + Zy**2 + 1)**(1.5))

    return H
