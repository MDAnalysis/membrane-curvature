from ast import NameConstant
import numpy as np


def derive_surface(n_cells, selection, max_width):
    """
    Derive surface from distribution of z coordinates in grid.

    Parameters
    ----------
    n_cells : int.
        number of cells in the grid of size `max_width`.
    selection : AtomGroup
        AtomGroup of reference selection to define the surface
        of the membrane.
    max_width : int.
        Maximum width of simulation box.

    Returns
    -------
    Returns set of z coordinates in grid.

    """
    NM_TO_ANGSTROM = 10  # Factor needed temporarily until we can make this code unit-aware or unitless
    z_ref = np.zeros((n_cells, n_cells))
    grid_z_coordinates = np.zeros((n_cells, n_cells))
    grid_norm_unit = np.zeros([n_cells, n_cells])

    # max_width *= NM_TO_ANGSTROM # scaling the max width of the box. Will remove
    factor = np.float32(n_cells / max_width)

    cell_xy_floor = np.int32(selection.positions[:, :2] * factor)
    z_coordinate = selection.positions[:, 2]

    for (l, m), z in zip(cell_xy_floor, z_coordinate):

        try:
            grid_z_coordinates[l, m] += z
            grid_norm_unit[l, m] += 1

        except IndexError:
            print("Atom outside grid boundaries. Skipping atom.")

    surface = avg_unit_cell(z_ref, grid_z_coordinates, grid_norm_unit)

    return surface


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

    """

    normed = grid_norm_unit > 0
    z_ref = grid_z_coordinates / grid_norm_unit
    z_ref[~normed] = np.nan

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
