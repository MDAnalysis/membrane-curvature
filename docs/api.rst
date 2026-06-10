**************************************
MembraneCurvature's API Documentation
**************************************

The MembraneCurvature API is built around the :class:`~membrane_curvature.base.MembraneCurvature` class, which loads trajectory and coordinate
data and provides routines to derive a surface from atom positions. The resulting surface is then used to compute mean and Gaussian curvature,
either for individual frames or averaged across multiple frames of an MD trajectory.

To use :class:`~membrane_curvature.base.MembraneCurvature`, users typically provide an MDAnalysis :class:`~MDAnalysis.core.universe.Universe`,
an :class:`~MDAnalysis.core.groups.AtomGroup`, and a method to derive the surface via the :attr:`~membrane_curvature.base.MembraneCurvature.surface_method`
parameter.

Currently, :class:`~membrane_curvature.base.MembraneCurvature` supports two surface derivation methods:

- :mod:`~membrane_curvature.fourier_surface` — the default method.
- :mod:`~membrane_curvature.binning_surface` — set explicitly by ``surface_method='binning'``.

The z-coordinates extracted from the selected ``AtomGroup`` are used to reconstruct the surface using the specified surface derivation method.
Once the surface has been generated, mean and Gaussian curvature are calculated using the Monge gauge formulas.

This page provides the API documentation for the ``MembraneCurvature`` class, the available surface derivation methods and their validators, and the curvature
calculation functions. For practical examples, refer to the :ref:`usage` and :ref:`tutorials` page.

MembraneCurvature class
=======================

The ``MembraneCurvature`` class provides the main entrypoint to calculate mean and Gaussian curvature from a surface derived from a reference
atom selection.

.. toctree::
   :maxdepth: 1

   api/base/base

Given a :class:`~MDAnalysis.core.universe.Universe` and an :class:`~MDAnalysis.core.groups.AtomGroup`,
:class:`~membrane_curvature.base.MembraneCurvature` reconstructs a surface from the heights of the atoms in the selected
:class:`~MDAnalysis.core.groups.AtomGroup`.
The specific operations used to derive the surface depend on the method selected by the user
in the :attr:`~membrane_curvature.base.MembraneCurvature.surface_method` parameter. 

.. warning::

   **The set of required parameters to run** :class:`~membrane_curvature.base.MembraneCurvature`
   **varies depending on the selected** ``surface_method`` **:**
   
   See the API documentation in :mod:`~membrane_curvature.base` for more details.

The reconstructed surfaces are then stored and used to calculate mean and Gaussian curvature. This calculation is available per-frame,
as well as an avergaed over the trajectory.


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

