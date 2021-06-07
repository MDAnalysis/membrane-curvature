Algorithm
=========================================================

Overall, this algorithm calculates mean and gaussian curvature from surfaces 
formed by a group of reference. By default, the group of reference to derive
surfaces are lipid headgroups.

First, the coordinates of each element of the group of references are 
calculated and then mapped to a grid of $m x n$ dimensions. To define the 
dimensions of the grid, we divide the length of the simulation box along 
the $x$ and $y$ axis over the unit cell width (default value 2 nm). 

Once the dimensions of the grid are defined, for every atom in the group of 
reference, the mapping function assigns an indexed cell in the grid. i.e. 
(x,y)↦[i,j]. For each indexed element in the grid, the respective $z$ 
coordinate is stored in a separate array. This operation is performed for 
every frame in the trajectory. After iterating over every frame of the 
trajectory, the mean value of the $z$ coordinate are computed. 

This analysis returns a 2-dimensional array which will be used to calculate
mean and gaussian curvature. 