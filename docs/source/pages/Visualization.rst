.. _visualization:

Visualization
=========================================================

This page provides examples of how to visualize membrane curvature via 
Matplotlib_. Two different approaches are suggested.

:ref:`plt.imshow`

:ref:`plt.contourf`

.. warning::
      Please note that when plotting via :func:`~matplotlib.pyplot.imshow` or
      :func:`~matplotlib.pyplot.contourf`, the orientation of the plot in the
      final rendering is determined by the ``origin`` argument. By default,
      in these Matplotlib_ functions, the origin is set to the ``(left, top)``
      corner. Additionally, for an array of shape ``(M, N)``, the first index
      ``M`` runs along the vertical axis, and the second index ``N`` runs along
      the horizontal axis.
      
      For these reasons, the returned arrays of surface, mean, and Gaussian
      curvature should be transposed, and ``origin=lower`` should be passed to
      :func:`~matplotlib.pyplot.imshow` or :func:`~matplotlib.pyplot.contourf`
      to generate the correct plots. This is a critical step when visualizing
      results obtained from MembraneCurvature analyses to avoid generating
      plots of membrane surface and curvature with the wrong orientation.

.. _plt.imshow:

1. imshow
----------------

A simple plot using :func:`~matplotlib.pyplot.imshow` can be obtained like so::


        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.imshow(avg_mean_curvature.T, cmap='bwr', interpolation='gaussian', origin='lower')
        ax.set_title('Mean Curvature')
        plt.show()

With :func:`~matplotlib.pyplot.imshow`, each element of the array is plotted as
a square in a matrix of `M x N` elements. Since we are interested in plotting
the array wwhere the first index runs along the horizontal axis, while the
second index runs along the vertical, we should transpose the array obtained
from the MembraneCurvature analysis. In the generated plot, the color of each
square is determined by the value of the corresponding array element and the
colormap of preference. 

For example, to visualize the results obtained in :ref:`membrane-only`, we can run:

.. ipython::
   :okwarning:
   
   In [0]: import MDAnalysis as mda
      ...: from membrane_curvature.base import MembraneCurvature
      ...: from MDAnalysis.tests.datafiles import Martini_membrane_gro
      ...: import matplotlib.pyplot as plt
   
   In [1]: universe = mda.Universe(Martini_membrane_gro)

   In [2]: curvature_upper_leaflet = MembraneCurvature(universe,
      ...:                                             select='resid 1-225 and name PO4',
      ...:                                             n_x_bins=8, 
      ...:                                             n_y_bins=8, 
      ...:                                             wrap=True).run()

   In [3]: mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

   In [4]: curvature_lower_leaflet = MembraneCurvature(universe,
      ...:                                             select='resid 226-450 and name PO4',
      ...:                                             n_x_bins=8, 
      ...:                                             n_y_bins=8, 
      ...:                                             wrap=True).run()

   In [5]: mean_lower_leaflet = curvature_lower_leaflet.results.average_mean
   
   In [6]: leaflets = ['Lower', 'Upper']

   In [7]: curvatures = [mean_lower_leaflet, mean_upper_leaflet]
   
   @savefig mycurvature.png width=8in
   In [8]: fig, [ax1, ax2] = plt.subplots(ncols=2, figsize=(6,3), dpi=200)
      ...: for ax, mc, lf in zip((ax1, ax2), curvatures, leaflets):
      ...:     ax.imshow(mc.T, origin='lower', interpolation='gaussian', cmap='seismic')
      ...:     ax.set_aspect('equal')
      ...:     ax.set_title('{} Leaflet'.format(lf))
      ...:     ax.axis('off')


.. _plt.contourf:

2. contourf
-------------------------------

You can use contour plots using :func:`~matplotlib.pyplot.contourf`. With this
approach, contour lines and filled contours of the obtained two-dimensional data
are plotted. A contour line connects points with the same curvature values.

When plotting using :func:`~matplotlib.pyplot.contourf`, an extra step is required
to perform an interpolation. We suggest using
:func:`scipy.ndimage.gaussian_filter` as in:

.. ipython::
   :okwarning:
   
   In [0]: from scipy import ndimage

   In [1]: leaflets = ['Lower', 'Upper']

   @savefig mycontours.png width=8in
   In [2]: fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(5,3))
      ...: for ax, mc, lf in zip((ax1, ax2), curvatures, leaflets):
      ...:     arr_ = ndimage.gaussian_filter(mc, sigma=1, order=0, mode='reflect')
      ...:     ax.contourf(arr_.T, 
      ...:                 cmap='bwr',
      ...:                 levels=30)
      ...:     ax.set_aspect('equal')
      ...:     ax.set_title('{} Leaflet'.format(lf))
      ...:     ax.axis('off')

.. _Matplotlib: https://matplotlib.org