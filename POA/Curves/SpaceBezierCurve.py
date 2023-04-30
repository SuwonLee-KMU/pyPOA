# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:17:29 PM GMT+9

import numpy as np
from pyPOA.POA.Curves.Bernstein import *

class SpaceBezierCurve():
    def __init__(self, coef_x, coef_y, coef_z):
        self.coef_x = coef_x
        self.coef_y = coef_y
        self.coef_z = coef_z

    def __call__(self, tau):
        func = self.get_curve_handle()
        return func(tau)
    
    def get_curve_handle(self):
        def BC(tau):
            fun_x = Bernstein_polynomial(self.coef_x)
            fun_y = Bernstein_polynomial(self.coef_y)
            fun_z = Bernstein_polynomial(self.coef_z)
            return np.array([fun_x(tau), fun_y(tau), fun_z(tau)])
        return BC
    
    def get_curve_Bernstein_coef(self): # aka control points
        return [[cx, cy, cz] for cx, cy, cz in zip(self.coef_x, self.coef_y, self.coef_z)]

