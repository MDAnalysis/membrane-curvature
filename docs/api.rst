**************************************
MembraneCurvature's API Documentation
**************************************

The MembraneCurvature API is built around the :class:`~membrane_curvature.base.MembraneCurvature` class, which loads trajectory and coordinate
files, and provides routines to derive a surface from atom positions. The resulting surface is then used to compute mean and Gaussian curvature,
either for individual frames or averaged across multiple frames of an MD trajectory.

To use :class:`~membrane_curvature.base.MembraneCurvature`, users must provide at least two inputs: an MDAnalysis
:class:`~MDAnalysis.core.universe.Universe`, and an :class:`~MDAnalysis.core.groups.AtomGroup` used as the reference to derive the surface and
calculate curvature.

:class:`~membrane_curvature.base.MembraneCurvature` supports three different methods to derive surfaces from reference atoms:

- :mod:`~membrane_curvature.fourier_surface` — the default method.
- :mod:`~membrane_curvature.binning_surface` — selected explicitly by ``surface_method='binning'``.
- :mod:`~membrane_curvature.binning_nearest_surface` - selected explicitly by ``surface_method='binning_nearest'``.

The height coordinates from the selected ``AtomGroup`` are extracted, and then used to reconstruct the surface using the specified surface derivation method.
Once the surface has been generated, mean and Gaussian curvature are calculated.

This page provides the API documentation for the ``MembraneCurvature`` class, the available surface derivation methods and their validators, and the curvature
calculation functions.

For practical examples, refer to the :ref:`usage` and :ref:`tutorials` page.

MembraneCurvature class
=======================

The :class:`~membrane_curvature.base.MembraneCurvature` class provides the main entrypoint to calculate mean and Gaussian curvature from a
surface derived from a reference atom selection.

.. toctree::
   :maxdepth: 1
   :hidden:

   api/base/base

Given a ``Universe`` and an ``AtomGroup``, :class:`~membrane_curvature.base.MembraneCurvature` reconstructs a surface from the heights
of the atoms in the selected ``AtomGroup`` of reference. 

.. note::

   The method selected by the user in the :attr:`~membrane_curvature.base.MembraneCurvature.surface_method` parameter defines:

   - The set of required parameters to run ``MembraneCurvature``.
   - The specific operations to derive the surface.
   
   See the API documentation in :class:`~membrane_curvature.base.MembraneCurvature` for more details.


Surface Methods
===============

:class:`~membrane_curvature.base.MembraneCurvature` supports three different methods to derive a surface from the selected
reference atoms:

.. toctree::
   :maxdepth: 1
   :hidden:

   api/surface_methods/fourier_surface
   api/surface_methods/binning_surface
   api/surface_methods/binning_nearest_surface
   api/surface_methods/padding
   api/surface_methods/fft_filtering

- :mod:`~membrane_curvature.fourier_surface` fits a periodic 2D Fourier sum to atom heights by linear least squares fit at each frame,
  evaluates the fitted surface, and obtains partial derivatives analytically from that sum.

- :mod:`~membrane_curvature.binning_surface` assigns atoms to a regular grid, stores the mean height per cell, and estimates partial derivatives
  using the physical bin spacing with finite differences.

- :mod:`~membrane_curvature.binning_nearest_surface` assigns each grid corner the ``z`` of the nearest lipid in the ``xy`` plane
  (minimum-image distances) and estimates partial derivatives with finite differences using the physical bin spacing.

Binning methods support a :mod:`~membrane_curvature.padding` option for orthorhombic boxes, which reduces finite difference artifacts, and a
a brick-wall filter to remove high-frequency noise from the height field. The filter is applied to the height field before calculating curvature.
:mod:`~membrane_curvature.fft_filtering` provides the functions to apply the filter and to resolve the pass-band limits.

Curvature
=========

The functions in :mod:`~membrane_curvature.curvature` include both analytical and numerical methods to estimate curvature from surface derivatives.

.. note::

   The calculation of mean and Gaussian curvature differs between the two methods:

   - With :mod:`~membrane_curvature.fourier_surface`, **curvature is calculated analytically** from the fitted surface.
   - With :mod:`~membrane_curvature.binning_surface` or :mod:`~membrane_curvature.binning_nearest_surface`, **curvature is calculated numerically**
     from the discrete height field using finite differences.

In both cases, the derivatives are then used to calculate mean and Gaussian curvature using the Monge-gauge formulas.

.. toctree::
   :maxdepth: 1
   :hidden:

   api/curvature/curvature

   
Validators
==========

The :mod:`~membrane_curvature.fourier_validators` and
:mod:`~membrane_curvature.padding_validators` are helper functions used to verify
that input parameters provided by the user are valid and compatible with the
``MembraneCurvature`` class and the surface derivation methods.
They are called automatically by the high-level entry points when
``MembraneCurvature`` runs with ``surface_method="fourier"`` or with
``padding=True`` on binning.

.. toctree::
   :maxdepth: 1
   :hidden:

   api/validators/fourier_validators
   api/validators/padding_validators

.. warning::

   **Validators are not considered a public API, should not be called directly, and are subject to change.**
   API documentation is provided for reference only.

