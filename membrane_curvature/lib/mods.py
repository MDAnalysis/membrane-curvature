import itertools as it
import numpy as np

def grid_map(coords, factor):
    """ Maps coordinates to grid.

    Parameters
    ----------
    x: float
        Value of x coordinate
    y: float
        Value of y coordinate

        Returns
        -------
        Return a tuple [l,m] with l,m as int.

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

    grid_count = np.zeros([n_cells, n_cells])

    for ts in universe.trajectory:

        grid_1 = np.zeros([n_cells, n_cells])
        grid_2 = np.zeros([n_cells, n_cells])

        factor = np.float32(n_cells / max_width)

        for bead in selection:
            x, y, z = bead.position/10

            index_l, index_m = grid_map( (x, y), factor)

            try:
                grid_1[index_l, index_m] += z
                grid_2[index_l, index_m] += 1

            except ValueError:
                print("Atom outside grid boundaries. Skipping atom {}".format(bead))

        for i, j in it.product(range(n_cells), range(n_cells)):
            if grid_2[i, j] > 0:
                z_Ref[i, j] += grid_1[i, j] / grid_2[i, j]
                grid_count[i, j] += 1

    for i, j in it.product(range(n_cells), range(n_cells)):
        if grid_count[i, j] > 0:
            z_Ref[i, j] /= grid_count[i, j]

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


def curvature(dict_reference, leaflets, n_cells):
    """
    Calculates mean curvature from Z cloud points.
    Parameters
    ----------
    dict_reference : dict { [index]:[z_coords]}.
        Dictionary with indexes [i,j] as keys and z_coordinates as values.
    leaflets : str. Default ["lower", "upper"]
        Leaflets of bilayer.
    n_cells : int. default 20.
        Number of cells in grid.
    Returns
    -------
    K, H : 2d-array, 2d-array
        Returns 2-dimensional array objects for mean
        and gaussian curvature, respectively.
    """

    H, K = [{key2: [] for key2 in leaflets} for i in range(2)]

    reference_avg = {key2:
                     {key3: [] for key3 in np.ndindex(n_cells, n_cells)}
                     for key2 in leaflets}

    for leaflet in leaflets:
        reference_avg[leaflet] = np.rot90(np.fliplr(dict_reference[leaflet]))

        K[leaflet] = gaussian_curvature(reference_avg[leaflet])
        H[leaflet] = mean_curvature(reference_avg[leaflet])

    else:
        print('No interpolation performed. Plot may display empty values.')
        return H, K
