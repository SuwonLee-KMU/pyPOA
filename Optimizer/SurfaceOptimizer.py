from pyPOA.POA.Polynomials.ComplexPolynomial import *
from pyPOA.POA.Surfaces.MSPN import MSPN
from scipy import optimize
import numpy as np

class SurfaceOptimizer():
    def __init__(self, wc, cc=None, crvc = None, w = [0.4, 0.3, 0.3]):
        self.WaypointConstraintsObj   = wc
        self.SpaceCurveConstraintsObj = cc
        self.CurvatureConstraintsObj  = crvc
        self.weights                  = w

        self.WaypointConstraintsObj.weights = self.weights[0:2]

    def cost(self, X):
        if self.SpaceCurveConstraintsObj is None:
            return self.WaypointConstraintsObj(X)
        else:
            return self.WaypointConstraintsObj(X) + self.SpaceCurveConstraintsObj(X)*self.weights[-1]

    def optimize(self, X0, method=None):
        res = optimize.minimize(self.cost, x0=X0, method=method)
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
        if self.SpaceCurveConstraintsObj is not None:
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