import MDAnalysis as mda
from membrane_curvature.tests.datafiles import GRO_PO4_SMALL
from membrane_curvature.base import MembraneCurvature


class MembraneCurvatureBenchmark():
    """
    Benchmark for MembraneCurvature class
    """

    def setup(self):
        self.u = mda.Universe(GRO_PO4_SMALL)
        self.sel = mda.select_atoms('name PO4')

    def time_surface(self):
        MembraneCurvature(self.u, select=self.sel).run()
