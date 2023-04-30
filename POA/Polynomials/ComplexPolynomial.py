# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:38:29 PM GMT+9
# Code: 

import numpy as np
import sympy
from pyPOA.POA.utilities import map_listed

class ComplexPolynomial():
    """
    복소다항식을 다루는 클래스.
    Properties:
        coefficients : 계수 리스트. 낮은 차수(0차)부터 쓴다. 계수는 복소수이다.        
    """
    def __init__(self, coef):
        self.coefficients = coef

    def __call__(self, z):
        func = self.function_handle()
        return func(z)

    def function_handle(self):
        @map_listed
        def func(z):
            dim = len(list(self.coefficients))
            zbasis = np.array([z**i for i in range(dim)])
            return np.dot(self.coefficients, zbasis)
        return func
    

class PreimageComplexPolynomial(ComplexPolynomial):
    def __init__(self, coef):
        self.coefficients = coef # Array of complex numbers. 낮은 차수(0차)부터 쓴다.

    def get_linear_complex_polynomial_handle(self):
        """
        입력한 계수에 대해 복소다항식 함수핸들을 출력하는 함수.
        """
        def func(z=0):
            dim = len(list(self.coefficients))
            zbasis = np.array([z**i for i in range(dim)])
            return np.dot(self.coefficients, zbasis)
        return func

    def __call__(self, z):
        func = self.get_linear_complex_polynomial_handle()
        return func(z)

    def symbolic(self):
        z = sympy.Symbol('z')
        func = self.get_linear_complex_polynomial_handle()
        return sympy.collect(func(z),z)
