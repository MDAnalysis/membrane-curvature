import numpy as np

from membrane_curvature.surface import get_z_surface
from membrane_curvature.curvature import mean_curvature, gaussian_curvature


from MDAnalysis.analysis.base import AnalysisBase
import MDAnalysis as mda


class MembraneCurvature(AnalysisBase):

    def __init__(self, universe, select='all',
                 n_x_bins=100, n_y_bins=100,
                 x_range=None,
                 y_range=None,
                 pbc=True, **kwargs):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.selection = universe.select_atoms(select)
        self.pbc = pbc
        self.n_x_bins = n_x_bins
        self.n_y_bins = n_y_bins
        self.x_range = x_range if x_range else (0, universe.dimensions[0])
        self.y_range = y_range if y_range else (0, universe.dimensions[1])

        # Raise if selection doesn't exist
        if (self.selection.n_residues) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup.")

        # Raise if range doesn't cover entire dimensions of simulation box
        if (self.x_range[1] < universe.dimensions[0]):
            raise ValueError("Grid range do not cover entire dimensions of "
                             "simulation box.")

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
        self.results.z_surface[self._frame_index] = get_z_surface(self.selection.positions,
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
