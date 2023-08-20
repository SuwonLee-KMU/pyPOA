import pyPOA.Optimizer.SurfaceOptimizer as SurfaceOptimizer
import pyPOA.POA.Curves.SpacePHCurve as SpacePHCurve
import pyPOA.POA.Curves.FrenetFrame as FrenetFrame
import numpy as np

class CurvatureConstraints():
    def __init__(self, max_curvature, qs):
        self.max_curvature = max_curvature
        self.planar_PHcurves = qs # shape = nwps

    def __call__(self, X, verbose=False):
        Cmax = self.eval_maximum_curvature(X)
        if verbose:
            print(f"> Computed maximum Curvature: \t{Cmax}, limit = {self.max_curvature}")
        return Cmax

    def eval_maximum_curvature(self, X):
        r = SurfaceOptimizer.MSPN_from_X(X)
        tau = np.linspace(0,1,41)
        MC = []
        for q in self.planar_PHcurves:
            p = SpacePHCurve(r, q)
            frenet = FrenetFrame(p)
            max_curvature = np.max(frenet.curvature()(tau))
            MC.append(max_curvature)
        return np.max(MC)

