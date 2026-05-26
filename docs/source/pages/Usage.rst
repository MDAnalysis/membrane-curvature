.. _usage:

Usage
=========================================================

This page includes a practical guide to using MembraneCurvature with the surface
derivation method of choice (binning or Fourier).
It also includes examples showing how curvature can be calculated for three different
simulation systems.

:ref:`surface-methods`

:ref:`examples-usage`

:ref:`membrane-only`

:ref:`membrane-protein`

        :ref:`membrane-protein-pr`

        :ref:`membrane-protein-no-pr`

.. warning::
   Examples included in this page show how to use MembraneCurvature
   using data files from `MDAnalysisTests`_. To run these examples,
   `MDAnalysisTests`_ must be installed.

.. _surface-methods:

1. Surface derivation methods
-----------------------------

:class:`~membrane_curvature.base.MembraneCurvature` uses an AtomGroup as a reference,
user-defined via the ``select`` parameter, to derive a surface and calculate mean and
Gaussian curvature.

There are two methods available to derive the surface:

- **Fourier** (``surface_method='fourier'``, default method) fits a truncated periodic 2D Fourier sum
  to atom heights by linear least squares at each frame, evaluates the fitted surface,
  and obtains partial derivatives analytically from that sum (no finite-difference 
  on the grid). Optional arguments ``fourier_m``, ``fourier_n``,
  tune the truncation for the Fourier fit and the least-squares solve.

- **Binning** (``surface_method='binning'``) assigns atoms to a
  regular ``n_x_bins`` x ``n_y_bins`` grid, stores the mean :math:`z` per cell,
  and estimates partial derivatives with :func:`numpy.gradient` using the
  physical bin spacing.

.. warning::

  Use ``fourier_m = fourier_n = 2`` (the constructor with default values) unless you need
  shorter-wavelength structure; increase ``fourier_m`` and ``fourier_n`` only while curvature
  improves systematically, rather than becoming noisier.

For copy-paste examples of both methods see:

- :ref:`membrane-only`: Fourier (default), then binning on the Martini bilayer.
- :ref:`membrane-protein-pr`: Fourier (default), then binning on the membrane-protein trajectory.


.. _examples-usage:

2. Examples of how to use MembraneCurvature to derive curvature profiles
------------------------------------------------------------------------

The following sections offer examples of how to use MembraneCurvature to derive curvature profiles in three types of systems:

- :ref:`membrane-only`
- :ref:`membrane-protein-pr`
- :ref:`membrane-protein-no-pr`

.. _membrane-only:

2.1. Membrane-only systems
------------------------------

In this example, we show a basic usage of MembraneCurvature in a system that
comprises a lipid bilayer of DPPC:CHOL using the Martini force field. Since we
have a bilayer, we select atoms of phospholipid head groups in the upper
leaflet only using the :attr:`~select` parameter and apply coordinate wrapping.

**Fourier surface method (default)**

We can calculate membrane curvature using the Fourier surface method by either
setting ``surface_method='fourier'`` with ``fourier_m=2``, ``fourier_n=2``, or
omitting ``surface_method``, ``fourier_m``, and ``fourier_n`` to rely on the defaults
(default ``fourier_m=2``, ``fourier_n=2``):

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                ).run()

    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

Note the code to calculate curvature for the upper leaflet with Fourier (default method)
is equivalent to:

.. code-block:: python

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='fourier',
                                                fourier_m=2,
                                                fourier_n=2,
                                                ).run()

.. tip::

  When using the Fourier method, ``wrap`` is not required. Periodic boundary conditions are handled
  inside the Fourier fit. 🙂
  Use ``wrap=True`` only with ``surface_method='binning'`` to pack atoms into the primary unit cell
  on raw trajectories.

**Binning surface method**

Alternatively, set ``surface_method='binning'`` and provide the values for
the binning grid ``n_x_bins`` and ``n_y_bins``. Note that you need to apply
coordinate wrapping with ``wrap=True``:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from MDAnalysis.tests.datafiles import Martini_membrane_gro

    universe = mda.Universe(Martini_membrane_gro)

    curvature_upper_leaflet = MembraneCurvature(universe,
                                                select='resid 1-225 and name PO4',
                                                surface_method='binning',
                                                n_x_bins=8,
                                                n_y_bins=8,
                                                wrap=True
                                                ).run()
  
    mean_upper_leaflet = curvature_upper_leaflet.results.average_mean

    gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian


You can find more detailed examples in the notebooks available in the :ref:`tutorials` page.


.. _membrane-protein:

2.2 Membrane-protein systems
----------------------------

.. tip::
        To improve sampling when passing raw trajectories:

        - In systems of membrane-only or membrane-protein with position restraints,
          set ``wrap=True`` to translate the atoms of the AtomGroup back in the unit cell.
        - In membrane-protein systems with no position restraints, set ``wrap=False`` and
          preprocess the trajectory with rotational/translational fit.

        Some points to keep in mind when calculating membrane curvature in membrane-only and
        membrane-protein systems with position restraints are addressed in this `blog post`_. 

.. _membrane-protein-pr:

2.2.1 Membrane-protein systems, protein with position restraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we have a simulation box comprising a copy of the Yiip
transporter, embedded in a lipid bilayer of POPE:POPG. Similar to the example
for membrane-only, we select the atoms for the upper leaflet to run the analysis.

**Fourier surface method (default)**

We can calculate membrane curvature using the default values with:

.. code-block:: python

   import MDAnalysis as mda
   from membrane_curvature.base import MembraneCurvature
   from MDAnalysis.tests.datafiles import XTC_MEMPROT, GRO_MEMPROT

   universe = mda.Universe(GRO_MEMPROT, XTC_MEMPROT)
   curvature_upper_leaflet = MembraneCurvature(universe,
                                               select='resid 297-517 and name P'
                                               ).run()

   avg_mean_curvature_upper_leaflet = curvature_upper_leaflet.results.average_mean

   avg_gaussian_curvature_upper_leaflet = curvature_upper_leaflet.results.average_gaussian

**Binning surface method**

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
                                               n_x_bins=2,
                                               n_y_bins=2,
                                               wrap=True
                                               ).run()

   avg_mean_curvature_upper_leaflet = curvature_upper_leaflet.results.average_mean

   avg_gaussian_curvature_upper_leaflet = curvature_upper_leaflet.results.average_gaussian


.. _membrane-protein-no-pr:

2.2.2. Membrane-protein systems, protein with no position restraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For membrane-protein systems where the simulation setup has no position
restraints on the protein, a trajectory preprocessing by the user is required.
If the goal is to assess membrane curvature induced by the protein, the 
preprocessed trajectory should have the protein centered in the simulation box 
with translational and rotational fit.

for example, in `Gromacs`_, the trajectory would be preprocessed with:

.. code-block:: bash

        gmx trjconv -pbc whole -ur compact -c
        gmx trjconv -fit rot+trans
        gmx trjconv -fit transxy

After you have preprocessed the trajectory, use ``wrap=False`` (the trajectory is
already fitted and centered).

**Fourier surface method (default)**

With the default ``surface_method='fourier'`` and ``fourier_m=2``, ``fourier_n=2``.
Omit ``surface_method``, ``fourier_m``, and ``fourier_n`` unless you need shorter wavelengths:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from membrane_curvature.tests.datafiles import XTC_MEMBPROT_FIT, GRO_MEMBPROT_FIT

    universe = mda.Universe(GRO_MEMBPROT_FIT, XTC_MEMBPROT_FIT)

    curvature_lower_leaflet = MembraneCurvature(universe,
                                                select='resid 2583-3042'
                                                ).run()

    avg_mean_curvature = curvature_lower_leaflet.results.average_mean

    avg_gaussian_curvature = curvature_lower_leaflet.results.average_gaussian

**Binning surface method**

Set ``surface_method='binning'`` with the values for ``n_x_bins`` and ``n_y_bins``, and in this
case set ``wrap=False`` to avoid the warning message since the trajectory is already fitted and centered:

.. code-block:: python

    import MDAnalysis as mda
    from membrane_curvature.base import MembraneCurvature
    from membrane_curvature.tests.datafiles import XTC_MEMBPROT_FIT, GRO_MEMBPROT_FIT

    universe = mda.Universe(GRO_MEMBPROT_FIT, XTC_MEMBPROT_FIT)

    curvature_lower_leaflet = MembraneCurvature(universe,
                                                select='resid 2583-3042',
                                                surface_method='binning',
                                                n_x_bins=10,
                                                n_y_bins=10,
                                                wrap=False
                                                ).run()

    avg_mean_curvature = curvature_lower_leaflet.results.average_mean

    avg_gaussian_curvature = curvature_lower_leaflet.results.average_gaussian

.. note::

        Since you are providing a preprocessed trajectory with translation/rotational fit 
        you can ignore the warning message: 
        ``WARNING   `wrap == False` may result in inaccurate calculation of membrane curvature.`` 


More information on how to visualize the results of the MDAnalysis Membrane 
Curvature tool can be found in the :ref:`visualization` page.

.. _`blog post`: https://ojeda-e.com/blog/Considerations-curvature-MD-simulations-PartI

.. _`MDAnalysisTests`: https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests

.. _`Gromacs`: https://www.gromacs.org/