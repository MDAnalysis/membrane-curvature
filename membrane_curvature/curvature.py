r"""
.. role:: raw-math(raw) :format: latex html

--------------------
Curvature
--------------------

In MembraneCurvature, we calculate Gaussian and mean curvature from a cloud of points.

Gaussian curvature is defined by

.. math:: K = \frac{\partial_{xx}\partial_{yy}-\partial_{xy}^2}
   {(1+\partial_x^2+\partial_y^2)^2}.

Mean curvature is defined by

.. math:: H =
    \frac{(1+\partial_x^2)\partial_{yy}+(1+\partial_y^2)\partial_{xx}-2\partial_x\partial_y\partial_{xy}}
    {2(1+\partial_x^2+\partial_y^2)^{3/2}}.


Notes
---------
Numpy cannot calculate the gradient for arrays with inner array of
`length==1` unless `axis=0` is specified. Therefore in the functions here included
for mean and Gaussian curvature, shape of arrays must be at least (2,2).
In general, to calculate a numerical gradients shape of arrays must be >=(`edge_order` +
1).


Functions
---------

"""

import numpy as np


def gaussian_curvature(Z):
    """
    Calculate gaussian curvature from Z cloud points.


    Parameters
    ----------
    Z: np.ndarray.
        Multidimensional array of shape (n,n).


    Returns
    -------
    K : np.ndarray.
        The result of gaussian curvature of Z. Returns multidimensional
        array object with values of gaussian curvature of shape `(n, n)`.

    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    K = (Zxx * Zyy - (Zxy ** 2)) / (1 + (Zx ** 2) + (Zy ** 2)) ** 2

    return K


def mean_curvature(Z):
    """
    Calculates mean curvature from Z cloud points.


    Parameters
    ----------
    Z: np.ndarray.
        Multidimensional array of shape (n,n).


    Returns
    -------
    H : np.ndarray.
        The result of gaussian curvature of Z. Returns multidimensional
        array object with values of gaussian curvature of shape `(n, n)`.

    """

    Zy, Zx = np.gradient(Z)
    Zxy, Zxx = np.gradient(Zx)
    Zyy, _ = np.gradient(Zy)

    H = (Zx**2 + 1) * Zyy - 2 * Zx * Zy * Zxy + (Zy**2 + 1) * Zxx
    H = -H / (2 * (Zx**2 + Zy**2 + 1)**(1.5))

    return H
