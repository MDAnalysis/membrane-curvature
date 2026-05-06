"""
MembraneCurvature
=======================================

:Author: Estefania Barreto-Ojeda
:Year: 2021
:Copyright: GNU Public License v3

MembraneCurvature calculates the mean and Gaussian curvature of
surfaces derived from a selection of reference.

Mean curvature is calculated in units of Å :sup:`-1` and Gaussian curvature
in units of Å :sup:`-2`.
"""

import numpy as np
import warnings
from .binning_surface import get_z_surface
from .curvature import (
    mean_curvature,
    gaussian_curvature,
    fourier_curvature,
)
from .fourier_surface import n_fourier_parameters

import MDAnalysis
from MDAnalysis.analysis.base import AnalysisBase

import logging

MDAnalysis.start_logging()

logger = logging.getLogger('MDAnalysis.MDAKit.membrane_curvature')


class MembraneCurvature(AnalysisBase):
    r"""
    MembraneCurvature is a tool to calculate membrane curvature.

    Parameters
    ----------
    universe : Universe or AtomGroup
        An MDAnalysis Universe object.
    select : str or iterable of str, optional.
        The selection string of an atom selection to use as a
        reference to derive a surface.
    wrap : bool, optional
        Apply coordinate wrapping to pack atoms into the primary unit cell.
    n_x_bins : int, optional, default: '100'
        Number of bins in grid in the x dimension.
    n_y_bins : int, optional, default: '100'
        Number of bins in grid in the y dimension.
    x_range : tuple of (float, float), optional, default: (0, `universe.dimensions[0]`)
        Range of coordinates (min, max) in the x dimension.
    y_range : tuple of (float, float), optional, default: (0, `universe.dimensions[1]`)
        Range of coordinates (min, max) in the y dimension.
    surface_method : {'binning', 'fourier'}, optional
        ``binning`` (default) bins atoms and uses :func:`numpy.gradient` for
        derivatives. ``fourier`` fits a periodic Fourier sum to atom positions
        each frame and evaluates Monge-gauge curvature from analytic derivatives
        on the same bin grid (bin centers); see :mod:`membrane_curvature.fourier_surface`.
    fourier_m : int, optional
        Maximum Fourier mode index in ``x`` when ``surface_method='fourier'``.
        Default ``2``.
    fourier_n : int, optional
        Maximum Fourier mode index in ``y`` when ``surface_method='fourier'``.
        Default ``2``.
    fourier_lstsq_rcond : float, optional
        ``rcond`` passed to :func:`numpy.linalg.lstsq` when ``surface_method='fourier'``.

    Attributes
    ----------
    results.z_surface : ndarray
        Surface derived from atom selection in every frame.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.mean : ndarray
        Mean curvature associated to the surface.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.gaussian : ndarray
        Gaussian curvature associated to the surface.
        Arrays of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.average_z_surface : ndarray
        Average of the array elements in `z_surface`.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_mean : ndarray
        Average of the array elements in `mean_curvature`.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_gaussian: ndarray
        Average of the array elements in `gaussian_curvature`.
        Each array has shape (`n_x_bins`, `n_y_bins`)


    See also
    --------
    :class:`~MDAnalysis.transformations.wrap.wrap`
        Wrap/unwrap the atoms of a given AtomGroup in the unit cell.

    Notes
    -----
    Use `wrap=True` to translates the atoms of your `mda.Universe` back
    in the unit cell. Use `wrap=False` for processed trajectories where
    rotational/translational fit is performed.

    For more details on when to use `wrap=True`, check the :ref:`usage`
    page.

    Fourier mode count: with ``surface_method='fourier'``, the number of fitted
    scalar parameters is :func:`~membrane_curvature.fourier_surface.n_fourier_parameters`
    ``(fourier_m, fourier_n)``. The selection must contain at least that many
    atoms each frame. Use ``fourier_m = fourier_n = 2`` unless you need shorter
    wavelengths; keep the atom count well above the parameter count, and raise
    ``fourier_m`` / ``fourier_n`` only while curvature improves systematically,
    not when it begins to track noise (see :doc:`api/fourier_surface`).

    The derived surface and calculated curvatures are available in the
    :attr:`results` attributes.

    The attribute :attr:`~MembraneCurvature.results.average_z_surface` contains
    the derived surface averaged over the `n_frames` of the trajectory.

    The attributes :attr:`~MembraneCurvature.results.average_mean_curvature` and
    :attr:`~MembraneCurvature.results.average_gaussian_curvature` contain the
    computed values of mean and Gaussian curvature averaged over the `n_frames`
    of the trajectory.

    Example
    -----------
    You can pass a universe containing your selection of reference::

        import MDAnalysis as mda
        from membrane_curvature.base import MembraneCurvature

        u = mda.Universe(coordinates, trajectory)
        mc = MembraneCurvature(u).run()

        surface =  mc.results.average_z_surface
        mean_curvature =  mc.results.average_mean
        gaussian_curvature = mc.results.average_gaussian

    The respective 2D curvature plots can be obtained using the `matplotlib`
    package for data visualization via :func:`~matplotlib.pyplot.contourf` or
    :func:`~matplotlib.pyplot.imshow`.

    For specific examples visit the :ref:`usage` page.
    Check the :ref:`visualization` page for more examples to plot
    MembraneCurvature results using :func:`~matplotlib.pyplot.contourf`
    and :func:`~matplotlib.pyplot.imshow`.

    """

    def __init__(
        self,
        universe,
        select='all',
        n_x_bins=100,
        n_y_bins=100,
        x_range=None,
        y_range=None,
        wrap=True,
        surface_method='binning',
        fourier_m=2,
        fourier_n=2,
        fourier_lstsq_rcond=None,
        **kwargs,
    ):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.ag = universe.select_atoms(select)
        self.wrap = wrap
        self.n_x_bins = n_x_bins
        self.n_y_bins = n_y_bins
        self.x_range = x_range if x_range else (0, universe.dimensions[0])
        self.y_range = y_range if y_range else (0, universe.dimensions[1])
        self.dx = (self.x_range[1] - self.x_range[0]) / n_x_bins
        self.dy = (self.y_range[1] - self.y_range[0]) / n_y_bins

        valid_methods = ('binning', 'fourier')
        if surface_method not in valid_methods:
            raise ValueError(f'surface_method must be one of {valid_methods}, got {surface_method!r}')
        self.surface_method = surface_method
        self.fourier_m = int(fourier_m)
        self.fourier_n = int(fourier_n)
        self.fourier_lstsq_rcond = fourier_lstsq_rcond
        if self.surface_method == 'fourier':
            if self.fourier_m < 0 or self.fourier_n < 0:
                raise ValueError('fourier_m and fourier_n must be non-negative')
            n_param = n_fourier_parameters(self.fourier_m, self.fourier_n)
            n_atom = len(self.ag)
            if n_atom < n_param:
                max_square_mode = max(int((np.sqrt(n_atom) - 1) // 2), 0)
                # If too few atoms for the passed modes fourier_m / fourier_n, raise error.
                # The message suggests the maximum fourier_m / fourier_n given the number
                # of atoms in the selection n_atom = len(self.ag)
                raise ValueError(
                    f"surface_method='fourier' needs at least {n_param} atoms in the selection "
                    f'(fourier_m={self.fourier_m}, fourier_n={self.fourier_n}), got {n_atom}. '
                    f'Suggested max modes for {n_atom} atoms is '
                    f'fourier_m = fourier_n = {max_square_mode}.'
                )

        # Raise if selection doesn't exist
        if len(self.ag) == 0:
            raise ValueError('Invalid selection. Empty AtomGroup.')

        # Only checks the first frame. NPT simulations not properly covered here.
        # Warning message if range doesn't cover entire dimensions of simulation box
        for dim_string, dim_range, num in [('x', self.x_range, 0), ('y', self.y_range, 1)]:
            if dim_range[1] < universe.dimensions[num]:
                msg = (
                    f'Grid range in {dim_string} does not cover entire '
                    'dimensions of simulation box.\n Minimum dimensions '
                    'must be equal to simulation box.'
                )
                warnings.warn(msg)
                logger.warning(msg)

        # Apply wrapping coordinates
        if not self.wrap:
            # Warning
            msg = (
                ' `wrap == False` may result in inaccurate calculation '
                'of membrane curvature. Surfaces will be derived from '
                'a reduced number of atoms. \n '
                ' Ignore this warning if your trajectory has '
                ' rotational/translational fit rotations! '
            )
            warnings.warn(msg)
            logger.warning(msg)

    def _prepare(self):
        # Initialize empty np.array with results
        self.results.z_surface = np.full((self.n_frames, self.n_x_bins, self.n_y_bins), np.nan)
        self.results.mean = np.full((self.n_frames, self.n_x_bins, self.n_y_bins), np.nan)
        self.results.gaussian = np.full((self.n_frames, self.n_x_bins, self.n_y_bins), np.nan)

    def _single_frame(self):
        # Apply wrapping coordinates
        if self.wrap:
            self.ag.wrap()
        if self.surface_method == 'binning':
            self.results.z_surface[self._frame_index] = get_z_surface(
                self.ag.positions,
                n_x_bins=self.n_x_bins,
                n_y_bins=self.n_y_bins,
                x_range=self.x_range,
                y_range=self.y_range,
            )
            self.results.mean[self._frame_index] = mean_curvature(
                self.results.z_surface[self._frame_index], self.dx, self.dy
            )
            self.results.gaussian[self._frame_index] = gaussian_curvature(
                self.results.z_surface[self._frame_index], self.dx, self.dy
            )
        else:
            z_f, mean_f, gauss_f = fourier_curvature(
                self.ag.positions,
                self.x_range,
                self.y_range,
                self.n_x_bins,
                self.n_y_bins,
                self.fourier_m,
                self.fourier_n,
                rcond=self.fourier_lstsq_rcond,
            )
            self.results.z_surface[self._frame_index] = z_f
            self.results.mean[self._frame_index] = mean_f
            self.results.gaussian[self._frame_index] = gauss_f

    def _conclude(self):
        self.results.average_z_surface = np.nanmean(self.results.z_surface, axis=0)
        self.results.average_mean = np.nanmean(self.results.mean, axis=0)
        self.results.average_gaussian = np.nanmean(self.results.gaussian, axis=0)
