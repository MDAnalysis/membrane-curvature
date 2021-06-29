import itertools as it
import numpy as np


def get_positions(index):
    """ 
    Get positions of bead by index 

    Parameters
    ----------
    index: Indexed assigned to bead in AtomGroup

    """
    x, y, z = index
    return x, y, z


def grid_map(coords, factor):
    """ 
    Maps (x,y) coordinates to unit cell in grid.

    Parameters
    ----------
    (x, y):  tuple
        Value of (x, y) coordinates
    factor:  float
        Mapping factor to assign grid.

    Returns
    -------
    Returns l, m  with l,m as int.

    """

    index_grid_l = int(abs(coords[0]) * factor)
    index_grid_m = int(abs(coords[1]) * factor)

    return index_grid_l, index_grid_m


def derive_surface(n_cells, selection, max_width):
    """
    Derive surface from distribution of z coordinates in grid.

    Parameters
    ----------
    n_cells : int.
        number of cells in the grid of size `max_width`.
    selection : AtomGroup
        AtomGroup of reference selected to define the surface
        of the membrane.
    max_width : int.
        Maximum width of simulation box.

    Returns
    -------
    Returns set of z coordinates in grid.

    """
    z_ref = np.zeros([n_cells, n_cells])
    grid_z_coordinates = np.zeros([n_cells, n_cells])
    grid_norm_unit = np.zeros([n_cells, n_cells])

    factor = np.float32(n_cells / max_width)

    for bead in selection:

        x, y, z = get_positions(bead)

        l, m = grid_map((x, y), factor)

        try:
            grid_z_coordinates[l, m] += z
            grid_norm_unit[l, m] += 1

        except ValueError:
            print("Atom outside grid boundaries. Skipping atom {}".format(bead))

    surface = avg_unit_cell(z_ref, n_cells, grid_z_coordinates, grid_norm_unit)

    return surface


def avg_unit_cell(z_ref, n_cells, grid_z_coordinates, grid_norm_unit):

    for i, j in it.product(range(n_cells), range(n_cells)):
        if grid_norm_unit[i, j] > 0:
            z_ref[i, j] += grid_z_coordinates[i, j] / grid_norm_unit[i, j]
        else:
            z_ref[i, j] = np.nan

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
