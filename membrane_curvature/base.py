import numpy as np

from membrane_curvature.lib.mods import derive_surface, mean_curvature, gaussian_curvature
import numpy as np


def __init__(self, universe, select,
             x_bins=100, y_bins=100,
             x_range=None,
             y_range=None,
             pbc=True):
    self.universe = universe
    self.selection = select
    self.pbc = pbc
    self.y_bins = y_bins
    self.x_bins = x_bins
    self.x_range = x_range if x_range else (0, universe.dimensions[0])
    self.y_range = y_range if y_range else (0, universe.dimensions[1])

    # Warns if range doesn't cover entire dimensions of simulation box
    if (x_range < universe.dimensions[0]) or (y_range < universe.dimensions[1]):
        raise ValueError("Minimum range must be ({}, {})").format(self.universe.dimensions[0],
                                                                  self.universe.dimenions[1])


def _prepare(self):

    # Initialize empty np.array(cumulative, count) of results
    cumulative = np.zeros((self.x_bins, self.y_bins))
    count = np.zeros((self.x_bins, self.y_bins))

    self.results.surface = (cumulative, count)
    self.results.mean = (cumulative, count)
    self.results.gaussian = (cumulative, count)


def _output(self, value_per_frame, cumulative_output, counter):
    """
    Adds np.arrays of mean curvature, gaussian curvature and surface to results.
    """

    nans_in_grid = np.isnan(value_per_frame)

    # Calculate cumulative per frame when no nans.
    cumulative_output += np.where(nans_in_grid, 0, value_per_frame)

    # Count when no nans.
    counter += np.where(nans_in_grid, 0, 1)


def _single_frame(self):

    surface_ = derive_surface(n_x_bins=self.n_cells_x, n_y_bins=self.n_cells_y,
                              x_range=(0, self.max_width_x), y_range=(0, self.max_width_y))

    mean_ = mean_curvature(surface_)
    gaussian_ = gaussian_curvature(surface_)

    # Save the results in NumPy array of shape (`n_x_bins`, `n_y_bins`).
    self._output(surface_, self.results.surface[0], self.results.surface[1])
    self._output(mean_, self.results.mean[0], self.results.mean[1])
    self._output(gaussian_, self.results.gaussian[0], self.results.gaussian[1])


def _conclude(self):

    # Calculate average of np.array over trajectory normalized over counter
    self.results.average_surface = self.results.surface[0] / self.results.surface[1]

    self.results.average_mean = self.results.mean[0] / self.results.mean[1]

    self.results.average_gaussian = self.results.surface[0] / self.results.surface[1]
