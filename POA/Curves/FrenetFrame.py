# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:15:11 PM GMT+9
# Code: 

import numpy as np

class FrenetFrame():
    """
    매개변수화된 3차원 공간곡선의 Frenet Frame을 계산하는 클래스.
    
    Property
        curve : 매개변수화된 3차원 곡선 객체.
            curve객체는 `.derivative(tau)`, `.second_derivative(tau)`, `third_derivative(tau)` 메서드를 가지고 있어야 한다.
    """
    def __init__(self, curve=None):
        self.curve = curve
    
    def tangent(self):
        @map_listed
        def fun(tau):
            r_prime = self.curve.derivative(tau)
            return r_prime / np.linalg.norm(r_prime)
        return fun
        
    def normal(self):
        @map_listed
        def fun(tau):
            r_prime = self.curve.derivative(tau)
            r_2prime = self.curve.second_derivative(tau)
            num = np.cross(np.cross(r_prime, r_2prime), r_prime)
            den = np.linalg.norm(np.cross(r_prime, r_2prime)) * np.linalg.norm(r_prime)
            return num / den
        return fun
    
    def binomal(self):
        @map_listed
        def fun(tau):
            r_prime = self.curve.derivative(tau)
            r_2prime = self.curve.second_derivative(tau)
            num = np.cross(r_prime, r_2prime)
            den = np.linalg.norm(num)
            return num / den
        return fun
    
    def curvature(self):
        @map_listed
        def fun(tau):
            r_prime = self.curve.derivative(tau)
            r_2prime = self.curve.second_derivative(tau)
            num = np.linalg.norm(np.cross(r_prime, r_2prime))
            den = np.linalg.norm(r_prime)**3
            return num / den
        return fun
    
    def torsion(self):
        @map_listed
        def fun(tau):
            r_prime = self.curve.derivative(tau)
            r_2prime = self.curve.second_derivative(tau)
            r_3prime = self.curve.third_derivative(tau)
            aux = np.cross(r_prime, r_2prime)
            if np.linalg.norm(aux) > 1e-2:
                num = np.dot(aux, r_3prime)
                den = np.linalg.norm(aux)**2
                return num / den
            else:
                return 0
        return fun

