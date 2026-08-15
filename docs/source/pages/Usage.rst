.. _usage:

======
Usage
======

This page includes a practical guide to using MembraneCurvature with the surface
derivation method of choice (Fourier, binning, or binning nearest). The available
methods are introduced using :ref:`membrane-only systems <membrane-only-systems>`,
followed by examples of membrane-protein systems
:ref:`with position restraints <membrane-protein-pr>`,
and :ref:`with no position restraints <membrane-protein-no-pr>`.

.. note::

   Examples in this page use data files from `MDAnalysisTests`_ and files from
   ``membrane_curvature.tests.datafiles``. To run these examples, `MDAnalysisTests`_
   and :ref:`MembraneCurvature <install_membrane_curvature>` must be installed.

.. _membrane-only-systems:

1. Membrane-only systems
========================

The examples in this section use a system that comprises a lipid bilayer of
DPPC:CHOL using the Martini force field. Since we have a bilayer, we select
atoms of phospholipid head groups in the upper leaflet only using the
:attr:`~select` parameter. Coordinate wrapping is applied only where required
by the selected surface method.

Some points to keep in mind when calculating membrane curvature in
membrane-only systems are addressed in this `blog post`_.

.. warning::

  MembraneCurvature does not identify leaflets. It relies on the user to pass
  the correct :class:`~MDAnalysis.core.groups.AtomGroup` with ``select``. Be aware 
  of the selection you provide: if both leaflets are included in the selection,
  they are treated as a single surface, hence producing incorrect results.

1.1 Fourier surface method (default)
------------------------------------

We can calculate membrane curvature using the Fourier surface method by either
setting ``surface_method='fourier'`` with default values ``fourier_m=2``,
``fourier_n=2``, or omitting ``surface_method``, ``fourier_m``, and ``fourier_n``
to rely on the defaults:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    # run with all default values
    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

.. warning::

  Use ``fourier_m=2`` and ``fourier_n=2``, the default values, unless you need shorter
  wavelengths. Increasing ``fourier_m`` and ``fourier_n`` may introduce noise in the fit
  rather than improving curvature results.

.. tip::

  When using the Fourier method, do not set ``wrap=True``. Because periodic boundary
  conditions are handled by the Fourier fit, passing ``wrap=True`` raises a :exc:`ValueError`.

The code above is equivalent to running with explicit parameters for ``surface_method``,
``fourier_m``, ``fourier_n``, and ``curvature_on`` with the default values:

.. code-block:: python

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='fourier',
                                                fourier_m=2,
                                                fourier_n=2,
                                                curvature_on='per_frame',
                                                ).run()

Note that by default, MembraneCurvature runs with ``curvature_on='per_frame'`` when
the parameter ``curvature_on`` is omitted. See section :ref:`curvature-on` for the
difference between ``'per_frame'`` and ``'average_surface'``.

.. important::

  The optional parameters ``padding`` and ``fft_filter`` are not available with
  ``surface_method='fourier'``.

Advanced: Tune the Fourier least-squares cutoff (``fourier_rcond``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Fourier surface is fit by solving a linear least-squares system with singular-value truncation `SVD`_.
The optional cutoff ``fourier_rcond`` controls which singular values are treated as "effectively zero" and
therefore removed from the solve. Smaller values keep more directions, and potentially noisier if
the system is underdetermined. Larger values regularize more aggressively.

In :class:`~membrane_curvature.base.MembraneCurvature`, pass this cutoff as ``fourier_rcond``:

.. code-block:: python

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='fourier',
                                                fourier_m=2,
                                                fourier_n=2,
                                                fourier_rcond=1e-12
                                                ).run()


If the effective rank of the design matrix is smaller than the number of
fitted parameters, a :class:`UserWarning` is emitted. The solver still returns a
valid least squares solution of smallest Euclidean norm, but the coefficients are
not uniquely determined by the data.

.. warning::

  ``fourier_rcond`` controls how aggressively we ignore poorly constrained combinations of
  Fourier coefficients. **We strongly recommend using** ``fourier_rcond`` **with its
  default value** ``None``.
  Larger values keep fewer singular-value directions (more stable / more
  regularized). Smaller values keep more directions (closer fit but potentially noisier).

  Rough intuition:

  - ``fourier_rcond=None``: sensible default. Uses NumPy's heuristic cutoff.
  - ``fourier_rcond=0``: truncate only *exactly* zero singular values.
  - ``fourier_rcond=1e-12`` or ``1e-10``: more aggressive truncation. It can reduce noise when the
    fit is underdetermined.

.. note::

  All Fourier least-squares steps (design matrix, SVD, coefficients) use
  64-bit floating point :func:`numpy.float64`. The cutoff defined by ``fourier_rcond`` is
  a **relative threshold**: singular values with :math:`s \le rcond \cdot s_{\max}` are dropped.
  The meaningful scale is therefore relative to the largest singular value :math:`s_{\max}`, not
  absolute coordinates or heights.

  With ``fourier_rcond=None``, the cutoff scales with the size of the least-squares problem, typically
  the number of atoms in the selection, or the number of fitted coefficients if that is larger.
  Passing a value much smaller than :math:`\sim 10^{-16}`, or much smaller than that automatic
  cutoff, usually has no visible effect. To smooth an underdetermined fit on purpose, use larger
  values such as ``1e-12`` or ``1e-10`` (see warning above).

1.2 Binning methods
-------------------

Binning methods calculate curvature using finite differences on the derived surfaces.
Two binning methods are available: ``'binning'`` and ``'binning_nearest'``.

For these methods, the optional parameter ``padding`` is available to reduce edge artifacts at the
simulation box boundaries, particularly for Gaussian curvature given that second order derivatives
amplify edge effects. Section :ref:`padding_usage` includes a usage example.

The optional parameter ``fft_filter`` is also available for binning methods, with a different
purpose. It applies a brick-wall filter in reciprocal space to remove high-frequency noise from
the averaged height field, and it is not a fix for edge artifacts. See section
:ref:`fft_filtering_usage` for a practical example.

1.2.1 Binning
^^^^^^^^^^^^^
Set ``surface_method='binning'`` to run membrane curvature with raw binning by
providing the values for the grid ``n_x_bins`` and ``n_y_bins``.

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    # run with binning and coordinate wrapping
    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='binning',
                                                n_x_bins=8,
                                                n_y_bins=8,
                                                wrap=True,
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

.. tip::

    When passing a raw trajectory with the ``surface_method='binning'``, set
    ``wrap=True``. With this setting, the :math:`x` and :math:`y` coordinates of the
    AtomGroup are wrapped into the primary unit cell. Atoms outside the primary cell
    are otherwise skipped during binning, which can leave bins empty or build the surface
    from fewer lipids than expected.

.. warning::

    With ``surface_method='binning'``, empty bins are stored as
    :data:`numpy.nan` and can propagate through the finite difference
    calculation of curvature. Edge and corner grid points also have fewer
    neighbours, which can introduce artifacts that are strongest in Gaussian
    curvature because it uses second derivatives.

    If Gaussian curvature is important for your analysis, consider running with
    ``padding=True`` to reduce edge artifacts at the simulation box boundaries.
    **Padding does not fill bins that have no reference atoms.**
    Empty bins remain as :data:`numpy.nan` when the grid is too coarse or the
    sampling is too sparse. See :ref:`Periodic Edge Padding<padding_usage>` for a
    practical example.

1.2.2 Binning nearest
^^^^^^^^^^^^^^^^^^^^^

Set ``surface_method='binning_nearest'`` to run membrane curvature using a binning method
that assigns each grid point the :math:`z` coordinate of the nearest lipid in the :math:`xy`
plane. The nearest lipid is searched over the whole selection, and the reference point of each
lipid is the center of its selected atoms.

.. tip::

    When running with ``surface_method='binning_nearest'``, omit the ``wrap`` parameter,
    or set it to ``False``. Coordinate wrapping is not valid with this method, because
    distances between grid points and lipids are calculated to the closest periodic image.
    Passing ``wrap=True`` raises a :exc:`ValueError`.

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    # run with binning nearest - no coordinate wrapping
    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='binning_nearest',
                                                n_x_bins=8,
                                                n_y_bins=8,
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

With this method, grid points sit at the bin corners, and not at the bin centers as in the
``binning`` and ``fourier`` methods. Therefore, for the same values of ``n_x_bins`` and
``n_y_bins``, the maps derived with ``binning_nearest`` are offset by half a bin with respect to
the maps derived with ``binning`` and ``fourier``. Because every grid point is assigned the
:math:`z` of a lipid, this method does not produce empty bins.

The details of the difference between the two binning methods are described in the
:ref:`binning_nearest_method` section of the Algorithm page.

.. important::

  Unlike raw binning, ``binning_nearest`` does not produce empty bins. Every grid point
  is assigned the :math:`z` of the nearest lipid in the selection, so the generated surface
  has a value at every bin.

  However, edge artifacts are still present in Gaussian curvature because second derivatives
  are calculated via finite differences. To reduce edge artifacts, consider running with
  ``padding=True``.

.. _padding_usage:

1.2.3 Periodic edge padding
^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the parameter ``padding``, MembraneCurvature adds a buffer around the
primary grid, calculates curvature on the expanded grid, and then clips the
results back to the initial ``n_x_bins`` x ``n_y_bins`` size.

.. warning::

  Padding is supported only for ``surface_method='binning'`` and
  ``surface_method='binning_nearest'``, and only when the simulation box
  is orthorhombic.

Padding is optional and is set to ``False`` by default. To apply edge periodic padding,
run with ``padding=True``:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    # run with binning and padding
    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='binning', # or 'binning_nearest'
                                                n_x_bins=8,
                                                n_y_bins=8,
                                                wrap=False,
                                                padding=True,
                                                ).run()

In this example, we use ``surface_method='binning'`` with ``wrap=False``,
because the periodic images added in the padded buffer provide the lateral
periodic boundary conditions. For this reason, MembraneCurvature does not
raise the usual ``wrap == False`` warning when padding is enabled.

.. tip::

  **Running with** ``padding=True`` **is particularly useful for analysis of Gaussian curvature,
  where edge artifacts can be significant.**

  Padding runs with ``edge_pad_bins=2`` by default. A buffer of 2 bins on each side
  is usually enough to fix edge artifacts in the calculation of second derivatives of
  Gaussian curvature. Values above ``edge_pad_bins=4`` are unlikely to change the curvature
  values and mainly increase runtime. In that case MembraneCurvature emits a
  :class:`UserWarning`. If ``edge_pad_bins`` is larger than ``n_x_bins`` or
  ``n_y_bins``, a :exc:`ValueError` is raised.

Alternatively, use ``edge_pad_bins`` if you need to change the width of the default buffer
regions of 2 bins on each side.

For example, if you need to change the width of the buffer regions to 3 bins on each side,
you can do:

.. code-block:: python

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='binning',
                                                n_x_bins=8,
                                                n_y_bins=8,
                                                wrap=False,
                                                padding=True,
                                                edge_pad_bins=3,
                                                ).run()

See :ref:`Periodic Edge Padding <binning-edge-padding>` in the Algorithm page for more details.


.. _fft_filtering_usage:

1.2.4 FFT filtering
^^^^^^^^^^^^^^^^^^^

FFT filtering is available for the binning methods ``surface_method='binning'`` and
``surface_method='binning_nearest'`` via the parameter ``fft_filter``. This option
applies a brick-wall passband to the averaged surface array.

Three filtering modes are available:

- ``fft_filter=None`` (default): no FFT filtering applied.
- ``fft_filter='auto'``: low-pass with pass band ``(0, 0.5 * q_Nyq)``, where
  :math:`q_{\mathrm{Nyq}} = \min(\pi/\Delta x,\ \pi/\Delta y)` from bin widths ``dx`` and ``dy``.
- ``fft_filter={'q': (q_low, q_high)}``: custom pass band in the reciprocal space ``(q_low, q_high)`` range in rad/Å (advanced).

An example of running an analysis with the auto filtering mode (``fft_filter='auto'``) looks like this:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)
    fft_auto_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='binning',
                                                n_x_bins=8,
                                                n_y_bins=8,
                                                wrap=True,
                                                fft_filter='auto').run()

.. tip::

  You can combine ``fft_filter`` with ``padding=True``. Enable ``fft_filter``
  only if you want to filter the averaged height field in reciprocal space to remove
  high-frequency noise.

.. important::

  When ``fft_filter`` is enabled and the parameter ``curvature_on`` is omitted,
  MembraneCurvature selects ``curvature_on='average_surface'``. Therefore,
  :attr:`~membrane_curvature.base.MembraneCurvature.results.average_mean` and
  :attr:`~membrane_curvature.base.MembraneCurvature.results.average_gaussian`
  are calculated from the filtered average surface.

  With an explicit ``curvature_on='per_frame'``, the filtered surface is still
  stored in
  :attr:`~membrane_curvature.base.MembraneCurvature.results.average_z_surface`,
  but the average curvature maps remain the averages of the unfiltered
  per-frame curvature arrays.


Advanced: Tuning the brick-wall filter
""""""""""""""""""""""""""""""""""""""

The brick-wall filter can be configured with custom cutoff values by passing a tuple
``(q_low, q_high)`` in rad/Å as ``fft_filter={'q': (q_low, q_high)}`` to define a custom
filter band. These values define the filter passband directly in reciprocal space.

.. warning::

  Use custom ``fft_filter={'q': (q_low, q_high)}`` values with caution. The
  pass band you choose directly affects the filtered average surface and the
  curvature maps derived from it. Poorly chosen bounds can produce inaccurate
  curvature values.

  We recommend ``fft_filter='auto'`` by default. The manual mode is available
  for flexibility, but you should only use it when you can justify the chosen
  bounds for your system.

To set the pass band relative to the Nyquist magnitude :math:`q_{\mathrm{Nyq}}`, compute
it with :func:`~membrane_curvature.fft_filtering.nyquist_q` from the bin widths ``dx`` and ``dy``,
then scale to set the pass band. For example, to set the passband to ``(0, 0.25 * q_Nyq)``:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature import nyquist_q
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    # set the number of bins
    n_x_bins = 8
    n_y_bins = 8

    # calculate the bin widths
    dx = universe.dimensions[0] / n_x_bins
    dy = universe.dimensions[1] / n_y_bins

    # calculate the Nyquist magnitude
    q_Nyq = nyquist_q(dx, dy)

    # run with binning and manual filter band
    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='binning',
                                                n_x_bins=n_x_bins,
                                                n_y_bins=n_y_bins,
                                                wrap=True,
                                                fft_filter={'q': (0, 0.25 * q_Nyq)},
                                                ).run()

The bin widths are also stored in the analysis instance as ``curvature_upper_leaflet.dx`` and
``curvature_upper_leaflet.dy``. You can use them to calculate the Nyquist magnitude of an
existing analysis with
``nyquist_q(curvature_upper_leaflet.dx, curvature_upper_leaflet.dy)``.

.. _curvature-on:

1.3 Average curvature mode (``curvature_on``)
---------------------------------------------

The parameter ``curvature_on`` controls how
:attr:`~membrane_curvature.base.MembraneCurvature.results.average_mean` and
:attr:`~membrane_curvature.base.MembraneCurvature.results.average_gaussian`
are obtained. It is available for all surface methods.

Two modes are available:

- ``curvature_on='per_frame'`` (default): mean and Gaussian curvature are
  calculated for every frame, and the average maps are the time averages of
  those per-frame arrays. Passing ``curvature_on=None`` is equivalent to
  ``'per_frame'``.
- ``curvature_on='average_surface'``: the average surface is built first, and
  mean and Gaussian curvature are calculated once from that average surface.

The scientific difference between the two paths is described in the Algorithm
page, section :ref:`calculate-curvature`.

To run with the default per-frame averaging mode, you can omit ``curvature_on``,
or set it explicitly:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                curvature_on='per_frame',
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

To calculate curvature from the average surface, set the parameter ``curvature_on='average_surface'``
explicitly:

.. code-block:: python

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                curvature_on='average_surface',
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian


.. important::

  In both examples above, the per-frame arrays
  :attr:`~membrane_curvature.base.MembraneCurvature.results.mean` and
  :attr:`~membrane_curvature.base.MembraneCurvature.results.gaussian` are still
  stored. The choice of ``curvature_on`` changes how the average maps are
  obtained, not whether per-frame results are stored.
  **The per-frame arrays are always stored in every run.**

.. _membrane-protein:

2. Membrane-protein systems
===========================

.. tip::

  When passing raw trajectories in membrane-protein systems:

  - **With position restraints:** set ``wrap=True`` when running with
    ``surface_method='binning'``.
  - **With no position restraints:** set ``wrap=False`` and preprocess the
    trajectory with a rotational and translational fit. The preprocessing steps
    are described in section :ref:`membrane-protein-no-pr`.

  Some points to keep in mind when calculating membrane curvature in
  membrane-protein systems with position restraints are addressed in this
  `blog post`_.

.. _membrane-protein-pr:

2.1 Protein with position restraints
------------------------------------

In this example, we have a simulation box comprising a copy of the Yiip
transporter, embedded in a lipid bilayer of POPE:POPG. Similar to the example
for membrane-only, we select the atoms for the upper leaflet to run the analysis.

2.1.1 Fourier surface method (default)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can calculate membrane curvature using the default values with:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import XTC_MEMPROT, GRO_MEMPROT

    universe = mda.Universe(GRO_MEMPROT, XTC_MEMPROT)
    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 297-517 and name P'
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

2.1.2 Binning surface method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The same trajectory and selection can use the binning method by setting
``surface_method='binning'`` with the values for ``n_x_bins`` and ``n_y_bins``
and apply coordinate wrapping with ``wrap=True``.

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import XTC_MEMPROT, GRO_MEMPROT

    universe = mda.Universe(GRO_MEMPROT, XTC_MEMPROT)
    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 297-517 and name P',
                                                surface_method='binning',
                                                n_x_bins=8,
                                                n_y_bins=8,
                                                wrap=True
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

For edge padding and for ``surface_method='binning_nearest'``, the same options
described in sections :ref:`padding_usage` and :ref:`fft_filtering_usage` apply.


.. _membrane-protein-no-pr:

2.2 Protein without position restraints
---------------------------------------

For membrane-protein systems where the simulation setup has no position
restraints on the protein, a trajectory preprocessing by the user is required.
If the goal is to assess membrane curvature induced by the protein, the
preprocessed trajectory should have the protein centered in the simulation box
with translational and rotational fit.

For example, in `Gromacs`_, the trajectory would be preprocessed with:

.. code-block:: bash

        gmx trjconv -pbc whole -ur compact -c
        gmx trjconv -fit rot+trans
        gmx trjconv -fit transxy

After you have preprocessed the trajectory, use ``wrap=False`` (the trajectory is
already fitted and centered).

2.2.1 Fourier surface method (default)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the default ``surface_method='fourier'`` and ``fourier_m=2``, ``fourier_n=2``:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from membrane_curvature.tests.datafiles import XTC_MEMBPROT_FIT, GRO_MEMBPROT_FIT

    universe = mda.Universe(GRO_MEMBPROT_FIT, XTC_MEMBPROT_FIT)

    curvature_lower_leaflet = MembraneCurvature(universe,
                                                select='resid 2583-3042 and name PO4'
                                                ).run()

    mean_lower_leaflet = curvature_lower_leaflet.results.average_mean

    gaussian_lower_leaflet = curvature_lower_leaflet.results.average_gaussian

2.2.2 Binning nearest surface method with padding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set ``surface_method='binning_nearest'`` with ``padding=True`` and the values
for ``n_x_bins`` and ``n_y_bins``. Keep ``wrap=False``, because coordinate
wrapping is not valid with ``surface_method='binning_nearest'``:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from membrane_curvature.tests.datafiles import XTC_MEMBPROT_FIT, GRO_MEMBPROT_FIT

    universe = mda.Universe(GRO_MEMBPROT_FIT, XTC_MEMBPROT_FIT)

    curvature_lower_leaflet = MembraneCurvature(universe,
                                                select='resid 2583-3042 and name PO4',
                                                surface_method='binning_nearest',
                                                n_x_bins=10,
                                                n_y_bins=10,
                                                wrap=False,
                                                padding=True,
                                                ).run()

    mean_lower_leaflet = curvature_lower_leaflet.results.average_mean

    gaussian_lower_leaflet = curvature_lower_leaflet.results.average_gaussian

.. note::

  With ``surface_method='binning_nearest'``, distances between grid corners and
  lipid reference points are calculated using the periodic simulation box.
  Therefore, a lipid across the opposite box boundary can be selected as the
  nearest lipid. Padding then adds periodic neighbors around the primary grid
  before finite differences are calculated.

More information on how to visualize the results of the MDAnalysis Membrane
Curvature tool can be found in the :ref:`visualization` page.

.. _`blog post`: https://ojeda-e.com/blog/Considerations-curvature-MD-simulations-PartI

.. _`MDAnalysisTests`: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests

.. _`Gromacs`: https://www.gromacs.org/

.. _`SVD`: https://en.wikipedia.org/wiki/Singular_value_decomposition
