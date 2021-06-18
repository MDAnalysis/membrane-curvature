import itertools as it
import numpy as np
import pickle
import math


def def_all_beads(lipid_types, leaflets, head_list, topology):
    """
    Select reference elements to derive membrane surface.


    Parameters
    ----------
    lipid_types : list.
        List of lipid types in simulation box.
    leaflets : str. Default ["lower", "upper"]
        Leaflets of bilayer.
    head_list : list
        List of indexes by leaflet. Provided by user in -ii and -io.
    topology: obj
        Topology from grofile.
    Returns
    -------
    [ element_1, ..., element_n] : list
        List with elements used as a reference to derive membrane surface.

    """

    dic_all_beads = {lf: {lt: [] for lt in lipid_types} for lf in leaflets}
    print('==== Lipid types in membrane ==== ')
    for lt in lipid_types:
        print('====>', lt)
        dic_all_beads['upper'][lt] = np.concatenate((topology.select('resname ' +
                                                                     lt +
                                                                     ' and index ' +
                                                                     str(head_list[0]) +
                                                                     ' to ' +
                                                                     str(head_list[1]) +
                                                                     ' and name PO4'), topology.select('resname ' +
                                                                                                       lt +
                                                                                                       ' and index ' +
                                                                                                       str(head_list[0]) +
                                                                                                       ' to ' +
                                                                                                       str(head_list[1]) +
                                                                                                       ' and name GM1'))).astype(int).tolist()
        dic_all_beads['lower'][lt] = np.concatenate((topology.select(
            'resname ' + lt + ' and index ' + str(head_list[1] + 1) + ' to ' + str(head_list[2]) + ' and name PO4'),
            topology.select(
            'resname ' + lt + ' and index ' + str(head_list[1] + 1) + ' to ' + str(head_list[2]) + ' and name GM1'))).astype(int).tolist()

        print("upper", len(dic_all_beads['upper'][lt]))
        print("lower", len(dic_all_beads['lower'][lt]))

    return dic_all_beads


def dict2pickle(name, dict_):
    """
    Exports values stored in a dictionary to a pickle file.

    Parameters
    ----------
    name : str,
        name of pickle file.
    dict_ : dict like {index -> value}
        dictionary to export.

    Returns
    -------
        Returns a pickled dictionary.

    """

    with open(name + ".pickle", 'wb') as pk:
        pickle.dump(dict_, pk, protocol=pickle.HIGHEST_PROTOCOL)


def core_fast(traj, jump, n_cells, leaflets, lipid_types, lipid_ref, max_width, prefix):
    """
    Runs core_fast_leaflet for each leaflet

    Parameters
    ----------
    traj: trajectory.
        MD trajectory to analyze.
    jump: int. Default 1.
        Skip <jump> number of frames in calculation.
    n_cells : int.
        number of cells in the grid of size `max_width`.
    leaflets : str. Default ["lower", "upper"]
        Leaflets of bilayer.
    lipid_types : list.
        List of lipid types in simulation box.
    lipid_Ref : list
        List of elements used as reference to define membrane surface.
    box_size : float
        Size of box (x dimension).
    max_width : int.
        Maximum width of simulation box.
    prefix : str. Default `None`
        Name of pickle
    Returns
    -------

    Returns a dictionary with name of leaflet as keys. The resulting values
    are nested dictionaries with tuples of [i,j] as keys, where
    0<i<`max_width`, and 0<j<`max_width`. Each key [i,j] has a list of
    `len=(n_frames)` containing the extracted z_coordinate a given [i,j]
    cell in each frame.

    """
    z_ref = {key1: np.zeros([n_cells, n_cells]) for key1 in leaflets}

    for leaflet in leaflets:
        core_fast_leaflet(z_ref[leaflet], leaflet, traj, jump, n_cells,
                          lipid_types, lipid_ref, max_width)

    dict2pickle(prefix, z_ref)

    return z_ref


def core_fast_leaflet(z_Ref, leaflet, traj, jump, n_cells, lipid_types,
                      lipid_ref, max_width):
    """
    Runs core_fast_leaflet for each leaflet

    Parameters
    ----------
    z_Ref : dict.
        Dictionary to store values in every iteration over frames.
    leaflet : str. {"lower", "upper"}
        Leaflets of bilayer.
    traj: trajectory.
        MD trajectory to analyze.
    jump: int. Default 1.
        Skip <jump> number of frames in calculation.
    n_cells : int.
        number of cells in the grid of size `max_width`.
    lipid_types : list.
        List of lipid types in simulation box.
    lipid_ref : list
        List of reference selection for atoms that define the surface
        of the membrane.
    box_size : float
        Size of box (x dimension).
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

    for frame in range(0, traj.n_frames, jump):
        grid_1 = np.zeros([n_cells, n_cells])
        grid_2 = np.zeros([n_cells, n_cells])

        factor = np.float32(n_cells / max_width)

        for lipid_type in lipid_types:

            for bead in lipid_ref[leaflet][lipid_type]:

                x, y, z = traj.xyz[frame, bead, :]

                l = int(abs(x) * factor)
                m = int(abs(y) * factor)

                try:
                    grid_1[l, m] += z
                    grid_2[l, m] += 1

                except BaseException:
                    pass

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
