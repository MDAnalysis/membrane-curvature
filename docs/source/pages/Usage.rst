Usage
=========================================================

MembraneCurvature can be used to derive membrane curvature profiles in three types of 
systems:

:ref:`membrane-only`.

:ref:`membrane-protein`.

        :ref:`membrane-protein-pr`.

        :ref:`membrane-protein-no-pr`.

.. _membrane-only:

1. Membrane-only systems
-----------------------------

A typical usage of MembraneCurvature to calculate the averaged mean curvature
over the trajectory of a system containing a membrane only is provided in the
following example::

        import MDAnalysis as mda
        from membrane_curvature.base import MembraneCurvature
        from MDAnalysis.tests.datafiles import TPR_MEMB, GRO_MEMB

        universe = mda.Universe(TPR_MEMB, GRO_MEMB)
        
        membrane_curvature = MembraneCurvature(universe, 
                                               select='name PO4', 
                                               wrap=True,
                                               n_x_bins=10,
                                               n_y_bins=10)

        membrane_curvature.run()

        avg_mean_curvature  = membrane_curvature.results.mean_curvature

.. _membrane-protein:

2.1 Membrane-protein systems
------------------------------


.. _membrane-protein-pr:

2.1.1 Membrane-protein systems, protein with position restraints
------------------------------------------------------------------
Similarly, in membrane-protein systems, when position restraints are applied on the protein,
we can calculate membrane curvature as::

        import MDAnalysis as mda
        from membrane_curvature.base import MembraneCurvature
        from MDAnalysis.tests.datafiles import TPR_MEMB_PROT, GRO_MEMB_PROT

        universe = mda.Universe(TPR_MEMB,_PROT GRO_MEMB_PROT)
        
        membrane_curvature = MembraneCurvature(universe, 
                                               select='name PO4', 
                                               wrap=True,
                                               n_x_bins=10,
                                               n_y_bins=10)

        membrane_curvature.run()

        avg_mean_curvature  = membrane_curvature.results.mean_curvature

.. note::
        Usage of MembraneCurvature in systems of :ref:`membrane-only` and 
        :ref:`membrane-protein-pr` is the same. 

Some points to keep in mind when calculating membrane curvature in :ref:`membrane-only`
and :ref:`membrane-protein-pr` are addressed in this `blog post`_. 

.. _membrane-protein-no-pr:

2.1.2. Membrane-protein systems, protein with no position restraints
---------------------------------------------------------------------

For membrane-protein systems where the simulation setup has no position
restraints on the protein, a trajectory preprocessing by the user is required.

After you have preprocessed the trajectory, a typical usage of membrane curvature is::

        import MDAnalysis as mda
        from membrane_curvature.base import MembraneCurvature
        from MDAnalysis.tests.datafiles import TPR_MEMB_PROT_FIT, GRO_MEMB_PROT_FIT

        universe = mda.Universe(TPR_MEMB_PROT_FIT, GRO_MEMB_PROT_FIT)
        
        membrane_curvature = MembraneCurvature(universe, 
                                               select='name PO4', 
                                               wrap=False,
                                               n_x_bins=10,
                                               n_y_bins=10)

        membrane_curvature.run()

        avg_mean_curvature  = membrane_curvature.results.mean_curvature


More information on how to visualize the results of the MDAnalysis Membrane 
Curvature tool can be found in the Visualization_ page.

.. _Visualization: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Visualization.html

.. _`blog post`: https://ojeda-e.github.io/blog/Considerations-curvature-MD-simulations-PartI/