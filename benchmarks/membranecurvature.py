import MDAnalysis as mda
import numpy as np

from membrane_curvature.base import MembraneCurvature
from membrane_curvature.tests.datafiles import GRO_PO4, GRO_PO4_SMALL, XTC_PO4

PO4_SELECT = 'name PO4'


def po4_universe(gro=GRO_PO4):
    """Universe with first frame loaded"""
    universe = mda.Universe(gro)
    universe.trajectory[0]
    return universe


def po4_frame0_context(gro=GRO_PO4):
    """Positions and box ranges for one PO4 frame"""
    universe = po4_universe(gro)
    atomgroup = universe.select_atoms(PO4_SELECT)
    return {
        'positions': atomgroup.positions.copy(),
        'x_range': (0.0, float(universe.dimensions[0])),
        'y_range': (0.0, float(universe.dimensions[1])),
    }


def dummy_height_grid(n_bins=100, seed=0):
    """Smooth height field for gradient / Monge curvature micro-benchmarks."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 2 * np.pi, n_bins)
    y = np.linspace(0, 2 * np.pi, n_bins)
    xx, yy = np.meshgrid(x, y)
    z = np.sin(xx) * np.cos(yy) + 0.1 * rng.standard_normal((n_bins, n_bins))
    return z, 0.1, 0.1


class MembraneCurvatureSmallBenchmark:
    def setup(self):
        self.universe = po4_universe(GRO_PO4_SMALL)

    def time_run_fourier_low_modes(self):
        # Nine PO4 atoms support at most fourier_m = fourier_n = 1.
        MembraneCurvature(self.universe, select=PO4_SELECT, fourier_m=1, fourier_n=1).run()

    def time_run_binning_default_wrap(self):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            surface_method='binning',
            n_x_bins=3,
            n_y_bins=3,
        ).run()

    def time_run_binning_no_wrap(self):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            surface_method='binning',
            n_x_bins=3,
            n_y_bins=3,
            wrap=False,
        ).run()


class MembraneCurvatureBenchmark:
    """Benchmark for PO4-only membrane (~900 atoms), regular grid sizes."""

    params = ([25, 50, 100],)
    param_names = ['n_bins']

    def setup(self, n_bins):
        self.universe = po4_universe(GRO_PO4)

    def time_run_fourier(self, n_bins):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            n_x_bins=n_bins,
            n_y_bins=n_bins,
        ).run()

    def time_run_binning_wrap(self, n_bins):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            surface_method='binning',
            n_x_bins=n_bins,
            n_y_bins=n_bins,
        ).run()

    def time_run_binning_no_wrap(self, n_bins):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            surface_method='binning',
            n_x_bins=n_bins,
            n_y_bins=n_bins,
            wrap=False,
        ).run()


class MembraneCurvatureFourierModesBenchmark:
    """Benchmarks for Fourier modes with PO4 selection."""

    params = ([(2, 2), (3, 3), (5, 5)],)
    param_names = ['modes']

    def setup(self, modes):
        self.universe = po4_universe(GRO_PO4)
        self.fourier_m, self.fourier_n = modes

    def time_run_fourier(self, modes):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            fourier_m=self.fourier_m,
            fourier_n=self.fourier_n,
        ).run()


class MembraneCurvatureTrajectoryBenchmark:
    """Multi-frame trajectory"""

    def setup(self):
        self.universe = mda.Universe(GRO_PO4, XTC_PO4)

    def time_run_fourier_trajectory(self):
        MembraneCurvature(self.universe, select=PO4_SELECT).run()

    def time_run_binning_trajectory(self):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            surface_method='binning',
            n_x_bins=3,
            n_y_bins=3,
        ).run()

    def peakmem_run_fourier_trajectory(self):
        MembraneCurvature(self.universe, select=PO4_SELECT).run()

    def peakmem_run_binning_trajectory(self):
        MembraneCurvature(
            self.universe,
            select=PO4_SELECT,
            surface_method='binning',
            n_x_bins=3,
            n_y_bins=3,
        ).run()
