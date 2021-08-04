"""
MembraneCurvature
=======================================

:Author: Estefania Barreto-Ojeda
:Year: 2021
:Copyright: GNU Public License v3

MembraneCurvature calculates the mean and Gaussian curvature of
surfaces derived from a selection of reference.
"""

import numpy as np
import warnings
from .surface import get_z_surface
from .curvature import mean_curvature, gaussian_curvature

import MDAnalysis
from MDAnalysis.analysis.base import AnalysisBase

import logging
MDAnalysis.start_logging()

logger = logging.getLogger("MDAnalysis.MDAKit.membrane_curvature")


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

    Attributes
    ----------
    results.z_surface : ndarray
        Surface derived from atom selection in every frame.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.mean_curvature : ndarray
        Mean curvature associated to the surface.
        Array of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.gaussian_curvature : ndarray
        Gaussian curvature associated to the surface.
        Arrays of shape (`n_frames`, `n_x_bins`, `n_y_bins`)
    results.average_z_surface : ndarray 
        Average of the array elements in `z_surface`. 
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_mean_curvature : ndarray 
        Average of the array elements in `mean_curvature`.
        Each array has shape (`n_x_bins`, `n_y_bins`)
    results.average_gaussian_curvature: ndarray 
        Average of the array elements in `gaussian_curvature`.
        Each array has shape (`n_x_bins`, `n_y_bins`)


    See also
    --------
    `MDAnalysis.transformations.wrap
    <https://docs.mdanalysis.org/1.0.0/documentation_pages/transformations/wrap.html>`_

    Notes
    -----
    Use `wrap=True` to translates the atoms of your `mda.Universe` back
    in the unit cell. Use `wrap=False` for processed trajectories where
    rotational/translational fit is performed.

    For more details on when to use `wrap=True`, check the `Usage
    <https://membrane-curvature.readthedocs.io/en/latest/source/pages/Usage.html>`_
    page.


    The derived surface and calculated curvatures are available in the
    :attr:`results` attributes.

    The attribute :attr:`~MembraneCurvature.results.average_z_surface` contains the
    derived surface averaged over the `n_frames` of the trajectory.

    The attributes :attr:`~MembraneCurvature.results.average_mean_curvature` and
    :attr:`~MembraneCurvature.results.average_gaussian_curvature` contain the computed
    values of mean and Gaussian curvature averaged over the `n_frames` of the
    trajectory.

    Example
    -----------
    You can pass a universe containing your selection of reference::

        import MDAnalysis as mda
        from MDAnalysis.analysis import MembraneCurvature

        u = mda.Universe(coordinates, trajectory)
        mc = MembraneCurvature(u).run()

        surface =  mc.results.average_z_surface
        mean_curvature =  mc.results.average_mean_curvature
        gaussian_curvature = mc.results.average_gaussian_curvature

    The respective 2D curvature plots can be obtained using the `matplotlib`
    package for data visualization via `imshow`. We recommend using the
    `gaussian` interpolation.


    """

    def __init__(self, universe, select='all',
                 n_x_bins=100, n_y_bins=100,
                 x_range=None,
                 y_range=None,
                 wrap=True, **kwargs):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.ag = universe.select_atoms(select)
        self.wrap = wrap
        self.n_x_bins = n_x_bins
        self.n_y_bins = n_y_bins
        self.x_range = x_range if x_range else (0, universe.dimensions[0])
        self.y_range = y_range if y_range else (0, universe.dimensions[1])

        # Raise if selection doesn't exist
        if len(self.ag) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup.")

        # Only checks the first frame. NPT simulations not properly covered here.
        # Warning message if range doesn't cover entire dimensions of simulation box
        for dim_string, dim_range, num in [('x', self.x_range, 0), ('y', self.y_range, 1)]:
            if (dim_range[1] < universe.dimensions[num]):
                msg = (f"Grid range in {dim_string} does not cover entire "
                       "dimensions of simulation box.\n Minimum dimensions "
                       "must be equal to simulation box.")
                warnings.warn(msg)
                logger.warn(msg)

        # Apply wrapping coordinates
        if not self.wrap:
            # Warning
            msg = (" `wrap == False` may result in inaccurate calculation "
                   "of membrane curvature. Surfaces will be derived from "
                   "a reduced number of atoms. \n "
                   " Ignore this warning if your trajectory has "
                   " rotational/translational fit rotations! ")
            warnings.warn(msg)
            logger.warn(msg)

    def _prepare(self):
        # Initialize empty np.array with results
        self.results.z_surface = np.full((self.n_frames,
                                          self.n_x_bins,
                                          self.n_y_bins), np.nan)
        self.results.mean = np.full((self.n_frames,
                                     self.n_x_bins,
                                     self.n_y_bins), np.nan)
        self.results.gaussian = np.full((self.n_frames,
                                         self.n_x_bins,
                                         self.n_y_bins), np.nan)

    def _single_frame(self):
        # Apply wrapping coordinates
        if self.wrap:
            self.ag.wrap()
        # Populate a slice with np.arrays of surface, mean, and gaussian per frame
        self.results.z_surface[self._frame_index] = get_z_surface(self.ag.positions,
                                                                  n_x_bins=self.n_x_bins,
                                                                  n_y_bins=self.n_y_bins,
                                                                  x_range=self.x_range,
                                                                  y_range=self.y_range)
        self.results.mean[self._frame_index] = mean_curvature(self.results.z_surface[self._frame_index])
        self.results.gaussian[self._frame_index] = gaussian_curvature(self.results.z_surface[self._frame_index])

    def _conclude(self):
        self.results.average_z_surface = np.nanmean(self.results.z_surface, axis=0)
        self.results.average_mean = np.nanmean(self.results.mean, axis=0)
        self.results.average_gaussian = np.nanmean(self.results.gaussian, axis=0)
