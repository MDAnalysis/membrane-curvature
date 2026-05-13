Algorithm
=========================================================

Overview
---------


MembraneCurvature calculates mean and Gaussian curvature of surfaces derived
from atoms of reference in 4 steps:

:ref:`select-atoms`

:ref:`choose-surface-method`

:ref:`derive-surface-curvature` 

:ref:`iterate`

A summary of the algorithm used in MembraneCurvature is shown in the following
diagram:

|diagram|

.. _select-atoms:

1. Select atoms of reference
-----------------------------

The first step in the algorithm consists of selecting atoms that will be used as
a reference to derive a surface. This selection will be contained in an
:class:`~MDAnalysis.core.groups.AtomGroup`. Typically in biological membranes,
lipid headgroups are the most common elements to use as an AtomGroup of
reference. 

|atoms|

.. _choose-surface-method:

2. Choose surface method
------------------------

Two surface-derivation methods are available, selected via the
``surface_method`` argument of
:class:`~membrane_curvature.base.MembraneCurvature`.

- :ref:`binning_method` (``surface_method='binning'``, default)

  Atoms are assigned to bins on a regular grid and the height of each bin
  is set to the mean :math:`z`-coordinate of its atoms. Partial derivatives
  are estimated numerically from the discrete height field using
  :func:`numpy.gradient` with the physical bin spacing.
  See :mod:`~membrane_curvature.binning_surface` for details.

- :ref:`fourier_method` (``surface_method='fourier'``)

  A truncated 2D Fourier series is fitted to atom heights by linear
  least squares each frame. Partial derivatives are evaluated analytically
  from the fitted series on the same bin-centre grid, avoiding
  finite-difference discretization error.
  See :mod:`~membrane_curvature.fourier_surface` for details.

.. _binning_method:

2.1. Binning method
-------------------

The binning method derives a surface by partitioning the simulation box into
a regular ``n_x_bins`` x ``n_y_bins`` grid along the :math:`x` and
:math:`y` directions. For each bin we compute the average :math:`z` coordinate
of the atoms that fall inside the bin to form the discrete height field.

The resulting discrete height field is then differentiated numerically
using :func:`numpy.gradient` to obtain the partial derivatives required for
curvature calculation. This is the default method
(``surface_method='binning'``) and requires no additional parameters
beyond the grid dimensions.

In the next section, we describe the details of the binning method.

2.1.1. Set grid
^^^^^^^^^^^^^^^
The dimensions of the grid are determined by the size of the simulation box
contained in the :class:`~MDAnalysis.core.universe.Universe`.
The grid covers a rectangular domain in the membrane plane. By default that
domain matches the `x` and `y` edges of the simulation box from MDAnalysis'
:class:`~MDAnalysis.core.universe.Universe`.The grid comprises ``n_x_bins`` x 
``n_y_bins`` bins along the `x` and `y` directions.

For every atom in the :class:`~MDAnalysis.core.groups.AtomGroup` of reference,
MembraneCurvature assigns it to a grid cell based on its `x` and `y` coordinates.
In practice, each coordinate pair ``(x, y)`` is mapped to a grid location
``[l, m]`` corresponding to a bin in the discretized membrane surface.

|grid|

Here, :math:`L_x` and :math:`L_y` denote the size of the simulation box in the 
`x` and `y` directions, respectively, while the size of the region covered by
the grid is represented by ``x_range`` and ``y_range``, along the `x` and `y`
directions, respectively. The spacing between grid points in each direction is
then determined by dividing these lengths by the number of bins.

.. note::
  Unless the user provides a different input, MembraneCurvature will determine
  the dimensions of the grid based on the size of the box on the first frame via
  :attr:`~MDAnalysis.core.universe.Universe.dimensions`.

  .. code-block:: python

      grid_dimension_x = (0, universe.dimensions[0])
      grid_dimension_y = (0, universe.dimensions[1])

2.1.2. Populate grid
^^^^^^^^^^^^^^^^^^^^^
Once the grid has been populated, the `z` coordinates of atoms assigned to each
cell are collected to form a height field over the grid.

Coordinates are converted to integer bin indices via scale factors

.. math::
   x\_factor = \frac{n\_x\_bins}{x_{max} - x_{min}}, \qquad
   y\_factor = \frac{n\_y\_bins}{y_{max} - y_{min}}

and then floored ``(index = floor(x * x_factor))``. Atoms with negative computed
indices are skipped and trigger a one-time detailed warning; atoms whose
indices exceed the grid bounds are also skipped and generate a warning. The
module uses a WarnOnce helper so the full diagnostic is shown on first
occurrence and a short message on subsequent occurrences.

Empty bins (zero samples) are represented as ``NaN`` in the returned
``(n_x_bins, n_y_bins)`` array: the implementation replaces zero counts
with ``NaN`` and divides summed z-values by the per-bin counts. As a result,
trajectory averages use ``numpy.nanmean`` and therefore ignore empty bins.

.. note::
  
  The binning routine itself does not apply periodic wrapping;
  the caller (for example :class:`~membrane_curvature.base.MembraneCurvature`)
  applies ``AtomGroup.wrap()`` when ``wrap=True`` before building the surface.
  Therefore whether atoms are wrapped into the primary box (and thus their
  bin assignment) depends on the wrapping option supplied to the analysis.

.. _fourier_method:

2.2 Fourier method
--------------------
The Fourier method fits a truncated, doubly periodic 2D Fourier series to
the atom height values in each frame. The basis is periodic on the
simulation box with periods :math:`L_x` and :math:`L_y`, so the fitted
surface is consistent with periodic boundary conditions in :math:`x` and
:math:`y`.

The full workflow comprises six steps, where the first four steps build and
solve a linear model that reconstructs the height field from plane-wave basis
functions.

2.2.1 Choose Fourier modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The first step is to choose a non-redundant set of Fourier modes with
:func:`~membrane_curvature.fourier_surface.fourier_mode_list` and determine the
total parameter count with :func:`~membrane_curvature.fourier_surface.n_fourier_parameters`.
This removes conjugate redundancy for real-valued surfaces and isolates the
mean term.

|fourier_modes|

2.2.2 Compute wavevectors
^^^^^^^^^^^^^^^^^^^^^^^^^^
Compute the wavevector components (:math:`k_x`, :math:`k_y`) for each mode
using :func:`~membrane_curvature.fourier_surface._compute_wavevector`. These
set the phase :math:`\phi = k_x x + k_y y` that appears in each cosine/sine
basis function.

.. warning::
  
  Because the derivatives are analytic, the Fourier method is not subject
  to finite-difference discretization error from the bin grid. Curvature
  still depends on the Fourier series truncation: use
  ``fourier_m = fourier_n = 2`` unless you need shorter wavelengths, and
  increase them only while curvature improves systematically rather than
  starts adding noise.

2.2.4 Build design matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^
Build the design matrix with :func:`~membrane_curvature.fourier_surface._build_fourier_matrix`.
Each row corresponds to an atom position and columns are the constant offset
followed by ``cos(k·r)``, ``sin(k·r)`` pairs for every retained mode; the
design matrix encodes the linear relation between Fourier coefficients and
observed heights.

|fourier_matrix|

2.2.5 Solve least-squares system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Solve the least-squares system for the Fourier coefficients using
:func:`~membrane_curvature.fourier_surface._fourier_fit_from_atoms` (which
calls :func:`numpy.linalg.lstsq`). Split the solution into the mean term and
the per-mode amplitudes with :func:`~membrane_curvature.fourier_surface._unpack_coefficients`.
The returned coefficients are then used to evaluate the fitted height and its
analytic derivatives on the bin-centre mesh.

The Fourier coefficients are obtained by linear least squares

.. math::

   \hat{\boldsymbol{\theta}} =
   \operatorname{arg\,min}_{\boldsymbol{\theta}}
   \lVert \mathbf{\Phi}\,\boldsymbol{\theta} - \mathbf{z} \rVert_2^2

via :func:`numpy.linalg.lstsq`. Because the model is linear in
:math:`\boldsymbol{\theta}`, no nonlinear optimisation is required.
If the effective rank of :math:`\mathbf{\Phi}` is smaller than :math:`P`,
a ``UserWarning`` is emitted; see :ref:`practical-considerations` for
guidance on choosing ``fourier_m`` and ``fourier_n``.

.. _derive-surface-curvature:

3. Derive surface and calculate surface derivatives
---------------------------------------------------

For every frame of the trajectory, the surface derived from the 
:class:`~MDAnalysis.core.groups.AtomGroup` is
calculated and stored in :attr:`~MembraneCurvature.results.z_surface`.
Similarly, the calculation of mean and Gaussian curvature is performed in every
frame and stored in :attr:`MembraneCurvature.results.mean_curvature` and
:attr:`MembraneCurvature.results.gaussian_curvature`, respectively.

The following sections describe the details of the two methods used
to derive the surface and calculate its derivatives.

.. _derive-surface:

3.1. Derive surface
^^^^^^^^^^^^^^^^^^^^

The surface is derived from atom positions using the selected method
and stored in :attr:`~MembraneCurvature.results.z_surface` for each frame.

With ``surface_method='binning'``, the height field is the discrete
:math:`N_x \times N_y` array of per-cell mean :math:`z` coordinates
assembled in :ref:`set-grid`. Bins containing no atoms are set to ``NaN``
and excluded from trajectory averages.

With ``surface_method='fourier'``, the height field is the fitted Fourier
series evaluated at bin centres after the least-squares fit described in
:ref:`fourier_method`. Because the representation is continuous and periodic
on the simulation box, every bin centre receives a value and no ``NaN``
entries arise.

|derive_surfaces_comparison|

.. _calculate-derivatives:

3.2. Calculate derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Five partial derivatives of the height field are required by the
curvature formulas: the first derivatives :math:`\partial_x` and
:math:`\partial_y`, the second derivatives :math:`\partial_{xx}`
and :math:`\partial_{yy}`, and the mixed derivative :math:`\partial_{xy}`.

|surf_fourier|

3.2a. Binning + finite differences (``surface_method='binning'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``surface_method='binning'``, these are estimated numerically from
the discrete height field using :func:`numpy.gradient` with the physical
spacings :math:`\Delta x` and :math:`\Delta y`, so that curvature values
are in physical units. Because this step uses finite differences, very
coarse grids may introduce discretization error.

.. note::

  These derivatives are evaluated using the actual grid spacing (``dx``, ``dy``),
  so that changes   in height are measured per unit distance in real space rather
  than per grid index. This makes curvature values physically meaningful and
  reduces their sensitivity to the chosen grid resolution.

  Curvature is computed from surface derivatives evaluated using the grid spacing
  (``dx``, ``dy``), ensuring results are expressed in physical units and are less
  sensitive to grid resolution. Because the derivatives are computed numerically,
  very coarse grids may still affect curvature estimates due to finite-difference
  discretization error.

For details on the binning method, see API documentation in
:mod:`~membrane_curvature.binning_surface` that describes every
associated functions.

3.2b. Fourier fit + analytic derivatives (``surface_method='fourier'``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``surface_method='fourier'``, the partial derivatives are evaluated
analytically from the fitted Fourier series (see :ref:`fourier_method`),
so they are not subject to finite-difference discretization error.
Curvature accuracy is instead governed by the Fourier series truncation.

.. warning::
  Use ``fourier_m = fourier_n = 2`` unless shorter wavelengths are needed, and
  increase the mode indices only while curvature improves systematically.

.. note::
  Because no finite-difference step is involved, the analytic derivatives
  are exact with respect to the fitted surface, regardless of grid resolution.

For details on the Fourier method, see API documentation in
:mod:`~membrane_curvature.fourier_surface` that describes every step with its
associated functions.

.. _mean-curvature:

3.3. Mean curvature
^^^^^^^^^^^^^^^^^^^^

Mean curvature :math:`H` is calculated from the five partial derivative
arrays using the Monge-gauge formula

.. math::

   H = \frac{(1+f_x^2)\,f_{yy} + (1+f_y^2)\,f_{xx}
             - 2\,f_x f_y f_{xy}}
            {2\,(1+f_x^2+f_y^2)^{3/2}}

via :func:`~membrane_curvature.curvature.mean_curvature_monge`.

The result has units Å :sup:`-1` and is stored in
:attr:`MembraneCurvature.results.mean_curvature` for each frame.

.. _gaussian-curvature:

3.4. Gaussian curvature
^^^^^^^^^^^^^^^^^^^^^^^^

Gaussian curvature :math:`K` is calculated from the same derivative
arrays using

.. math::

   K = \frac{f_{xx}\,f_{yy} - f_{xy}^2}
            {(1+f_x^2+f_y^2)^{2}}

via :func:`~membrane_curvature.curvature.gaussian_curvature_monge`.


Similarly, the calculation of mean and Gaussian curvature is performed
in every frame and the result has units Å :sup:`-2` and is stored in
:attr:`MembraneCurvature.results.gaussian_curvature` for each frame.


.. warning::

   The Monge-gauge formulas in steps 3.3 and 3.4 are exact. There is no
   small-gradient approximation applied. Both methods feed the same
   five derivative arrays into the same formulas; the only difference is
   how those derivatives were obtained in step :ref:`calculate-derivatives`.

.. _iterate:

4. Average over frames
-----------------------------------

The attributes :attr:`MembraneCurvature.results.average_mean` and
:attr:`MembraneCurvature.results.average_gaussian` contain the computed
values of mean and Gaussian curvature averaged over all the 
:attr:`~n_frames` in the trajectory. 

After performing the average over frames, the information of average surface,
mean, and Gaussian curvature are stored in the 
:attr:`MembraneCurvature.results.average_z_surface<membrane_curvature.base.MembraneCurvature.results.average_z_surface>`,
:attr:`MembraneCurvature.results.average_mean`, and
:attr:`MembraneCurvature.results.average_gaussian` arrays, respectively.
Each array has shape ``(n_x_bins, n_y_bins)``.

|avg_frames|

.. |diagram| image:: ../_static/Algorithm_v200.png
  :width: 800
  :alt: MembraneCurvature_diagram

.. |atoms| image:: ../_static/AtomsReference.png
  :width: 600
  :alt: atoms_ref

.. |grid| image:: ../_static/grid.png
  :width: 600
  :alt: Grid

.. |fourier_modes| image:: ../_static/Fourier_modes.png
  :width: 800
  :alt: FourierModes

.. |fourier_matrix| image:: ../_static/Fourier_matrix.png
  :width: 800
  :alt: FourierMatrix

.. |derive_surfaces_comparison| image:: ../_static/derive_surfaces_comparison.png
  :width: 800
  :alt: DeriveSurfacesComparison

.. |surf_fourier| image:: ../_static/DeriveSurfCurv_Fourier.png
  :width: 800
  :alt: CurvDiaagram

.. |avg_frames| image:: ../_static/AvgFrames.png
  :width: 800
  :alt: avgFrames
