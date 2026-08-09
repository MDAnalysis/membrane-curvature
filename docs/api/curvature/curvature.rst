.. automodule:: membrane_curvature.curvature
   :no-members:

.. _numerical_derivatives:

Numerical derivatives
^^^^^^^^^^^^^^^^^^^^^

Functions called to calculate curvature when running the binning methods: ``'binning'``,
``'binning_nearest'``.

.. autofunction:: mean_curvature
.. autofunction:: gaussian_curvature

With optional padding (``padding=True``), curvature is calculated numerically on the padded
surface with:

.. autofunction:: curvature_from_primary_with_edge_pad
.. autofunction:: curvature_with_edge_pad

.. _analytic_derivatives:

Analytic derivatives
^^^^^^^^^^^^^^^^^^^^

Functions called to calculate curvature when running the Fourier surface method.

.. autofunction:: fourier_curvature
.. autofunction:: fourier_curvature_from_coefficients
.. autofunction:: fourier_curvature_from_theta

Helpers
^^^^^^^

Monge gauge formulas from precomputed partial derivatives.

.. autofunction:: mean_curvature_monge
.. autofunction:: gaussian_curvature_monge
