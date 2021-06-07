Visualization
=========================================================

The results of the membrane curvature analysis can be visualized via 
Matplotlib_ and Pandas_. Since DataFrame structures are 2-dimensional tabular
data, back-transformation to Numpy arrays is inexpensive. 

Two different approaches are suggested.

#. Via `imshow`, from the Matplotlib package.
#. Via `contourf` or `contour`, from the Matplotlib package.

With the first, each element of the array is plotted as a square in a matrix 
of $m x n$ elements. The color of each square is determined by the value of 
the corresponding array element and the colormap of preference. 

With the second, contour lines and filled contours of the obtained two-dimensional
data are plotted. A contour line connects points with the same curvature values.

A divergent colormap is suggested. 

.. _Matplotlib: https://matplotlib.org
.. _Pandas: https://pandas.pydata.org
