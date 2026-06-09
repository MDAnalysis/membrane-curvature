**************************************
MembraneCurvature's API Documentation
**************************************

The API of MembraneCurvature centers on the :class:`~membrane_curvature.base.MembraneCurvature` class. It loads topology and coordinate
data and provides routines to derive a surface from atom positions, then used to compute mean and Gaussian curvature either
per frame or averaged over frames.

To run MembraneCurvature, typically users provide an MDAnalysis :class:`~MDAnalysis.core.universe.Universe`, an :class:`~MDAnalysis.core.groups.AtomGroup`,
and  a method to derive the surface according to the :attr:`~membrane_curvature.base.MembraneCurvature.surface_method`.
The selected ``AtomGroup`` are used to extract the z-coordinates that are used to reconstruct the surface according to the chosen method. Once the surface
is derived, it is used to calculate mean and Gaussian curvature using the Monge gauge formulas.

This page provides the API documentation for the ``MembraneCurvature`` class, the two surface derivation methods with its validators, and the curvature functions.
For usage examples, refer to the :ref:`usage` page and :ref:`tutorials` page.

MembraneCurvature class
=======================

The ``MembraneCurvature`` class ties the workflow together to calculate mean and Gaussian curvature from a surface derived from a reference
atom selection.

It accepts an :class:`~MDAnalysis.core.universe.Universe` and an :class:`~MDAnalysis.core.groups.AtomGroup` as input, and from these,
inputs it builds a surface from the selected reference atoms, then stores the resulting surface,and calculates mean curvature, and Gaussian curvature
for each frame as well as trajectory averages.

.. toctree::
   :maxdepth: 1

   api/base/base


Surface Methods
===============

Using the selected reference atoms, :class:`~membrane_curvature.base.MembraneCurvature` uses two methods to derive the surface
used to calculate mean and Gaussian curvature:

.. toctree::
   :maxdepth: 1

   api/surface_methods/fourier_surface
   api/surface_methods/binning_surface

**The Fourier method** fits a periodic 2D Fourier sum to atom heights by linear least squares at each frame,
evaluates the fitted surface, and obtains partial derivatives analytically from that sum.

**The binning method** assigns atoms to a regular grid and stores the mean :math:`\bar{z}` per cell, and estimates partial derivatives
with numpy.gradient using the physical bin spacing with finite differences.

.. attention::

   **The calculation of mean and Gaussian curvature differs between the two methods:**

   - With Fourier, **curvature is calculated analytically** from the fitted surface.
   - With binning, **curvature is calculated numerically** from the discrete height field using finite differences.
   For details on the two methods, see the API documentation in :mod:`~membrane_curvature.fourier_surface` and :mod:`~membrane_curvature.binning_surface`.

Curvature
=========

The functions in :mod:`membrane_curvature.curvature` include both analytical and numerical methods to estimate surface derivatives from a height field.

For the binning method (``surface_method='binning'``), ``MembraneCurvature`` estimates surface derivatives from a height field using finite differences with
:func:`numpy.gradient` while for the Fourier method (set by default, ``surface_method='fourier'``), the surface derivatives are estimated from the fitted
surface using analytical methods.

In both cases, the derivatives are then used to calculate mean and Gaussian curvature using the Monge-gauge formulas.

.. toctree::
   :maxdepth: 1

   api/curvature/curvature


Validators
==========

Validators are helper functions used to verify that input parameters provided by the user are valid and compatible with the
:class:`~membrane_curvature.base.MembraneCurvature` class and the surface derivation methods.

.. toctree::
   :maxdepth: 1

   api/validators/fourier_validators

.. warning::

   We provide the API documentation for the validators **for reference only.**

   Do not call them directly. **Validators are not intended to be called directly by the user**.
   They are called automatically by the high-level entry points in the surface derivation methods:
   :mod:`~membrane_curvature.fourier_surface` and :mod:`~membrane_curvature.binning_surface`.

