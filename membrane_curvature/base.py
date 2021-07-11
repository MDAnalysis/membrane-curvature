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
from membrane_curvature.surface import get_z_surface
from membrane_curvature.curvature import mean_curvature, gaussian_curvature


from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda


class MembraneCurvature(AnalysisBase):

    r"""
    MembraneCurvature is a tool to calculate membrane curvature.

    Parameters
    ----------
    universe : Universe or AtomGroup
        An MDAnalysis Universe object.
    select : str or iterable of str, optional. Default: 'all' 
        The selection string of an atom selection to use as a 
        reference to derive a surface. 
    pbc : bool, optional. Default: 'True'
        Apply periodic boundary conditions.
    n_x_bins : int, optional. Deafult: '100'
        Number of bins in grid in the `x` dimension. 
    n_y_bins : int, optional. Default: '100'
        Number of bins in grid in the `y` dimension. 
    x_range : tuple of (float, float), optional. Deafult: (0, `universe.dimensions[0]`)
        Range of coordinates (min, max) in the `x` dimension.
    y_range : tuple of (float, float), optional. Deafult: (0, `universe.dimensions[1]`)
        Range of coordinates (min, max) in the `y` dimension. 

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

    Returns
    -------
    average_z_surface : ndarray
        Average of the array elements in `z_surface`
    average_mean_curvature :
        Average of the array elements in `mean_curvature`
    average_gaussian_curvature:
        Average of the array elements in `gaussian_curvature`

    """

    def __init__(self, universe, select='all',
                 n_x_bins=100, n_y_bins=100,
                 x_range=None,
                 y_range=None,
                 pbc=True, **kwargs):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.ag = universe.select_atoms(select)
        self.pbc = pbc
        self.n_x_bins = n_x_bins
        self.n_y_bins = n_y_bins
        self.x_range = x_range if x_range else (0, universe.dimensions[0])
        self.y_range = y_range if y_range else (0, universe.dimensions[1])

        # Raise if selection doesn't exist
        if len(self.ag) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup.")

        # Warning message if range doesn't cover entire dimensions of simulation box
        for dim_string, dim_range, num in [('x', self.x_range, 0), ('y', self.y_range, 1)]:
            if (dim_range[1] < universe.dimensions[num]):
                warnings.warn(f"Grid range in {dim_string} do not cover entire "
                              "dimensions of simulation box.\n Minimum dimensions "
                              "must be equal to simulation box.")

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
        # Populate a slice with np.arrays of surface, mean, and gaussian per frame
        self.results.z_surface[self._frame_index] = get_z_surface(self.ag.positions,
                                                                  n_x_bins=self.n_x_bins,
                                                                  n_y_bins=self.n_y_bins,
                                                                  x_range=self.x_range,
                                                                  y_range=self.y_range)

        self.results.mean[self._frame_index] = mean_curvature(self.results.z_surface[self._frame_index])

        self.results.gaussian[self._frame_index] = gaussian_curvature(self.results.z_surface[self._frame_index])

    def _conclude(self):
        # Calculate mean of np.array of surface, mean, gaussian
        self.results.average_z_surface = np.nanmean(self.results.z_surface, axis=0)

        self.results.average_mean = np.nanmean(self.results.mean, axis=0)

        self.results.average_gaussian = np.nanmean(self.results.gaussian, axis=0)