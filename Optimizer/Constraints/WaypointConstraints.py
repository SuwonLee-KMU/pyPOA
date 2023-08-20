import pyPOA.Optimizer.SurfaceOptimizer as SurfaceOptimizer
import pyPOA.POA.Curves.SpacePHCurve as SpacePHCurve
import numpy as np

class WaypointConstraints():
    def __init__(self, Ps, Ts, qs, w=[0.5, 0.5]):
        self.space_waypoints_pos = Ps # shape = nwps+1 x 3
        self.space_waypoints_tan = Ts # shape = nwps+1 x 3
        self.planar_PHcurves     = qs # shape = nwps
        self.weights             = np.array(w)

    def __call__(self, X, verbose=False):
        pos_error_norm = self.eval_pos_error(X)
        tan_error_norm = self.eval_tan_error(X)
        if verbose:
            print(f"> Position error norm: \t{pos_error_norm}")
            print(f"> Tangent error norm: \t{tan_error_norm}")
        return self.weights[0]*pos_error_norm + self.weights[1]*tan_error_norm

    def eval_pos_error(self, X):
        nwps          = np.shape(self.space_waypoints_pos)[0]-1
        r             = SurfaceOptimizer.MSPN_from_X(X)
        pos_err_array = []
        for i in range(nwps):
            q = self.planar_PHcurves[i]
            p = SpacePHCurve(r, q)
            pos_error = np.array(p(0)) - np.array(self.space_waypoints_pos[i])
            pos_err_array.append(pos_error)
        else:
            q = self.planar_PHcurves[i]
            p = SpacePHCurve(r, q)
            pos_error = np.array(p(1)) - np.array(self.space_waypoints_pos[i+1])
            pos_err_array.append(pos_error)
        error_norm_array = np.linalg.norm(pos_err_array,axis=1)
        return np.sum(error_norm_array)

    def eval_tan_error(self, X):
        nwps = np.shape(self.space_waypoints_pos)[0]-1
        r    = SurfaceOptimizer.MSPN_from_X(X)
        G    = []
        E    = []
        for i in range(nwps):
            q = self.planar_PHcurves[i]
            p = SpacePHCurve(r, q)
            grad = np.array(p(0.001)) - np.array(p(0))
            G.append(grad)

        T_u = np.array([v/np.linalg.norm(v) for v in 
            self.space_waypoints_tan[:-1]])
        G_u = np.array([v/np.linalg.norm(v) for v in G])
        # G_u, T_u = [v / np.linalg.norm(v) for v in [G, self.space_waypoints_tan[:-1]]]
        # print(G_u, T_u, sep='\t')
        for i in range(nwps):
            ang_btw  = np.arccos(np.clip(np.dot(G_u[i], T_u[i]), -1.0, 1.0))
            E.append(ang_btw)   # array of directional error (rad)
        return np.sum(np.abs(E))
