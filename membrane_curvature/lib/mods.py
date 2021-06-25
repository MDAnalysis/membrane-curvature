import itertools as it
import numpy as np


def grid_map(coords, factor):
    """ Maps (x,y) coordinates to unit cell in grid.

    Parameters
    ----------
    (x, y):  tuple
        Value of (x, y) coordinates
    factor:  float
        Mappign factor to assign grid.

        Returns
        -------
        Returns l, m  with l,m as int.

    """

    index_grid_l = int(abs(coords[0]) * factor)
    index_grid_m = int(abs(coords[1]) * factor)

    return index_grid_l, index_grid_m


def core_fast_leaflet(universe, z_Ref, n_cells, selection, max_width):
    """
    Runs core_fast_leaflet for selected surface

    Parameters
    ----------
    universe: MDA Universe
        Provides coordinates and trajectory of the coordinates.
    z_Ref : dict.
        Dictionary to store values in every iteration over frames.
    n_cells : int.
        number of cells in the grid of size `max_width`.
    selection : AtomGroup
        AtomGroup of reference selection to define the surface
        of the membrane.
    max_width : int.
        Maximum width of simulation box.

    Returns
    -------

    Returns a dictionary with tuples of [i,j] as keys, where 0<i<`max_width`,
    and 0<j<`max_width`. The resulting values for each key are a list of
    `len=(n_frames)` containing the extracted z_coordinate a given [i,j] cell
    in each frame.

    """

    grid_count_frames = np.zeros([n_cells, n_cells])

    for ts in universe.trajectory:

        grid_z_coordinates = np.zeros([n_cells, n_cells])
        grid_norm_unit = np.zeros([n_cells, n_cells])

        factor = np.float32(n_cells / max_width)

        for bead in selection:
            x, y, z = bead.position/10

            l, m = grid_map((x, y), factor)

            try:
                grid_z_coordinates[l, m] += z
                grid_norm_unit[l, m] += 1

            except ValueError:
                print("Atom outside grid boundaries. Skipping atom {}".format(bead))

        for i, j in it.product(range(n_cells), range(n_cells)):
            if grid_norm_unit[i, j] > 0:
                z_Ref[i, j] += grid_z_coordinates[i, j] / grid_norm_unit[i, j]
                grid_count_frames[i, j] += 1

    for i, j in it.product(range(n_cells), range(n_cells)):
        if grid_count_frames[i, j] > 0:
            z_Ref[i, j] /= grid_count_frames[i, j]

        else:
            z_Ref[i, j] = np.nan


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
