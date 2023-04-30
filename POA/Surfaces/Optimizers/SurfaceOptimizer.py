# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:42:54 PM GMT+9
# Code: 

import numpy as np
from scipy import optimize
from POA.Polynomials import PreimageComplexPolynomial
from POA.Surfaces import MSPN
from POA.Curves import SpacePHCurve, FrenetFrame

class SurfaceOptimizer():
    def __init__(self, wc, cc, crvc = None, w = [0.4, 0.3, 0.3]):
        self.WaypointConstraintsObj   = wc
        self.SpaceCurveConstraintsObj = cc
        self.CurvatureConstraintsObj = crvc
        self.weights                  = w

        self.WaypointConstraintsObj.weights = self.weights[0:2]

    def cost(self, X):
        return self.WaypointConstraintsObj(X) + self.SpaceCurveConstraintsObj(X)*self.weights[-1]

    def optimize(self, X0):
        res = optimize.minimize(self.cost, x0=X0)
        print(res['success'])
        Xopt = res['x']
        opt_MSPN = SurfaceOptimizer.MSPN_from_X(Xopt)
        return res, opt_MSPN

    def cost_test(self, X):
        return self.SpaceCurveConstraintsObj(X)*self.weights[-1]
    
    def optimize_test(self, X0, options=None):
        constraint_pos = optimize.NonlinearConstraint(fun=self.WaypointConstraintsObj.eval_pos_error, lb=0, ub=0.001) 
        constraint_tan = optimize.NonlinearConstraint(fun=self.WaypointConstraintsObj.eval_tan_error, lb=0, ub=0.001) 
        if self.CurvatureConstraintsObj is None:
            constraints = [constraint_pos, constraint_tan]
        else:
            constraint_crv = optimize.NonlinearConstraint(
                fun=self.CurvatureConstraintsObj.eval_maximum_curvature, 
                lb=0, 
                ub=self.CurvatureConstraintsObj.max_curvature
            )
            constraints = [constraint_pos,constraint_tan,constraint_crv]
        if not options:
            options = dict(maxiter=1000, disp=True)
        res = optimize.minimize(self.cost_test, x0=X0, options=options, constraints=constraints)
        print(res['success'])
        Xopt = res['x']
        opt_MSPN = SurfaceOptimizer.MSPN_from_X(Xopt)
        return res, opt_MSPN  

    def check_constraints(self, X):
        self.WaypointConstraintsObj(X, verbose=True)
        self.SpaceCurveConstraintsObj(X, verbose=True)
        return None

    def check_constraints_test(self, X):
        self.WaypointConstraintsObj(X, verbose=True)
        self.SpaceCurveConstraintsObj(X, verbose=True)
        if self.CurvatureConstraintsObj:
            self.CurvatureConstraintsObj(X, verbose=True)
        return None

    @staticmethod
    def X2Params(X):
        r0            = X[:3]                                # constant term of the space curve
        preimage_coef = np.reshape(X[3:], newshape=[4,-1])   # poly. coefficients of complex functions, Phi.
        u_coef_real = preimage_coef[0,:]
        u_coef_imag = preimage_coef[1,:]
        v_coef_real = preimage_coef[2,:]
        v_coef_imag = preimage_coef[3,:]
        u_coef = u_coef_real + u_coef_imag*1j
        v_coef = v_coef_real + v_coef_imag*1j
        return r0, u_coef, v_coef

    @staticmethod
    def initialize_params(n:int):
        if n < 1:
            raise ValueError('n should be larger than 0.')
        nparam = 4*n + 7
        X0 = np.random.randn(nparam)
        return X0

    @staticmethod
    def MSPN_from_X(X):
        r0, u_coef, v_coef = SurfaceOptimizer.X2Params(X)
        u = PreimageComplexPolynomial(u_coef)
        v = PreimageComplexPolynomial(v_coef)
        w = PreimageComplexPolynomial([1+0j])
        r = MSPN(u, v, w, r0)
        return r

    def eval_cost(self):
        pass

    def compute_cost(self, X):
        pass


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

