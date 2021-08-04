Visualization
=========================================================

The results of the membrane curvature analysis can be visualized via 
Matplotlib_.

Two different approaches are suggested.

1. Via imshow_.

A simple plot using `imshow` can be obtained by::


        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.imshow(avg_mean_curvature, cmap='bwr', interpolation='gaussian', origin='lower')
        ax.set_title('Mean Curvature')
        plt.show()

With imshow_, each element of the array is plotted as a square in a matrix 
of `m x n` elements. The color of each square is determined by the value of 
the corresponding array element and the colormap of preference. 


2. Via `contourf` or `contour`.

You can use contour plots using `contourf`::


        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.contourf(avg_mean_curvature, cmap='bwr, origin='lower')
        ax.set_title('Mean Curvature')
        plt.show()


With this approach, contour lines and filled contours of the obtained two-dimensional
data are plotted. A contour line connects points with the same curvature values.

A divergent colormap is suggested. 

.. _Matplotlib: https://matplotlib.org
.. _imshow: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html
