Algorithm
=========================================================

The membrane curvature MDAnalysis module calculates the mean and Gaussian curvature 
from surfaces derived by a group of reference. By default, the group of reference 
implemented are lipid headgroups.

The algorithm locates all the atoms in the group of reference and maps each element 
of the group to a grid of `m x n` dimensions. The dimensions of the grid are 
determined by the length of the simulation box in the first frame of the 
trajectory.

The grid is then divided in `m x n` elements, with unit cells of width `uw` = 2nm 
(default value). In this way, for every atom in the group of reference an indexed cell 
in the grid is assigned according to their respective `x` and `y` coordinates.
i.e. `(x, y) ↦ [i, j]`. 

|grid|

Once the grid is populated with the atoms in the group of reference, the associated
`z` coordinate of each atom is stored in an array assigned to each [i, j] unit cell.
The length of these arrays will depend on the number of atoms of reference that
are mapped to a given `[i, j]` unit cell. Subsequently, once all values of `z` are
stored, the mean value of `z` stored in that array is calculated. Since this operation 
is performed for every frame, the length of the array containing the mean `z` position
will be determined by the number of frames iterated in the trajectory.

.. code-block:: python

    len(reference_atoms[i,j]) == len(u.trajectory)

After iterating over every frame of the trajectory, the mean value of the stored `z`
coordinate are calculated. These mean values of `z` are the ones used to derive the
surface.

With the defined surface, then values of mean (`H`) and Gaussian (`K`) curvature
are calculated. This analysis returns a 2-dimensional array of `m x n` with the values of mean
mean and Gaussian curvature. The length of this array is `m x n`.

More information on how to visualize the results of the MDAnalysis Membrane 
Curvature toolkit can be found in the Visualization_ page .

.. |grid| image:: ../_static/gridmap.png
  :width: 400
  :alt: GridMap

.. _Visualization: