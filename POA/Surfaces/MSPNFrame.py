# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:41:43 PM GMT+9
# Code: 

import numpy as np
from POA.utilities import map_listed

class MSPNFrame():
    def __init__(self, SpacePHCurveObj):
        self.SpacePHCurveObj = SpacePHCurveObj
    
    def normal(self):
        @map_listed
        def fun(tau):
            aux = self.SpacePHCurveObj.q(tau, form="complex")
            normal = self.SpacePHCurveObj.r.normal(aux)
            return normal / np.linalg.norm(normal)
        return fun
    
    def tangent(self):
        tan_handle = self.SpacePHCurveObj.binomial_tangent_handle()
        @map_listed
        def fun(tau):
            tan = tan_handle(tau)
            return tan / np.linalg.norm(tan)
        return fun

    def binomal(self):
        fun_normal = self.normal()
        fun_tangent = self.tangent()
        @map_listed
        def fun(tau):
            normal = fun_normal(tau)
            tangent = fun_tangent(tau)
            binomal = np.cross(normal, tangent)
            return binomal / np.linalg.norm(binomal)
        return fun
