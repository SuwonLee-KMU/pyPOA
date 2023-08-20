import pyPOA.POA.Curves.SpacePHCurve as SpacePHCurve
import pyPOA.Optimizer.SurfaceOptimizer as SurfaceOptimizer
import numpy as np

class SpaceCurveConstraints():
    def __init__(self, Ld, qs):
        self.length_desired = Ld
        self.planar_PHcurves= qs # shape = nwps

    def __call__(self, X, verbose=False):
        L  = self.eval_length(X)
        Le = np.abs(L-self.length_desired)
        if verbose:
            print(f"> Length error: \t{Le}   (Desired / Actual = {self.length_desired}, {L})")
        return Le

    def eval_length(self, X):
        nwps = np.shape(self.planar_PHcurves)[0]
        r    = SurfaceOptimizer.MSPN_from_X(X)
        Ls   = []
        tau  = np.linspace(0,1,100)
        for i in range(nwps):
            q = self.planar_PHcurves[i]
            p = SpacePHCurve(r, q)
            p_eval = p(tau)
            Ls.append(SpaceCurveConstraints.compute_length(p_eval))
        L = np.sum(Ls)
        return L
        
    # 3차원 곡선의 길이 수치적분하는 함수
    @staticmethod
    def compute_length(P):
        P = np.reshape(P, newshape=[3,-1])
        T = np.diff(P, axis=1)
        T = np.concatenate([T,np.reshape(T[:,0],newshape=[3,1])],axis=1)
        N = np.linalg.norm(T,axis=0)
        L = sum(N)
        return L

