# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:40:48 PM GMT+9
# Code: 

import numpy as np
import sympy

class PlanarPHCurve():
    """
    2차원 평면에서 정의된 PH곡선.

    Attributes:
    H : HermiteInterpolation 객체
    
    Methods:
    tangent(tau) : 곡선의 매개변수 tau에서의 접선벡터를 계산하여 출력한다.
    get_curve_handle(form) : 곡선의 매개변수를 주면 곡선의 위치를 출력하는 함수 핸들을 출력하는 메서드.
                            form은 함수 핸들이 출력하는 형식을 복소수나 실수 벡터 중 하나로 지정해 준다.
    get_derivative_handle(form) : 곡선의 매개변수를 주면 곡선의 접선벡터를 출력하는 함수 핸들을 출력하는 메서드.
    get_curve_coef : 곡선의 계수를 출력하는 함수. 다항식의 계수이며 높은차수부터 낮은차수 순서이다.
    symbolic : 곡선을 심볼릭 표현으로 출력한다.
    """
    def __init__(self, HermiteInterpolation):
        self.H = HermiteInterpolation

    def __call__(self, tau, form="complex"):
        func = self.get_curve_handle(form=form)
        return func(tau)
    
    def tangent(self, tau, form="complex"):
        fun = self.get_derivative_handle(form=form)
        return fun(tau)
        
    def get_curve_handle(self, form="complex"):
        if form == "complex":
            def curve(tau):
                q_real, q_imag = self.H.eval(tau)
                return q_real + q_imag*1j
        elif form == "real":
            def curve(tau):
                return self.H.eval(tau)
        return curve

    def get_curve_coef(self):
        return np.flip(PlanarPHCurve.BernsteinPoly_to_power_coef(self.H.cps),axis=1) # 높은차수부터 낮은차수 순서로 출력한다.

    def symbolic(self):
        tau = sympy.Symbol('tau')
        func = self.get_curve_handle(form="real")
        return func(tau)

    def get_derivative_handle(self, form="complex"):
        """
        곡선의 매개변수에 대한 이차 미분 함수 핸들을 출력하는 함수. 복소수 표현식이거나 실수 표현식 중 하나를 선택할 수 있다.
        """
        if form == "complex":
            def derivative(tau):
                x_prime, y_prime = self.H.tangent(tau)
                return x_prime + y_prime*1j
        elif form == "real":
            def derivative(tau):
                return self.H.tangent(tau)
        return derivative
                
    def get_second_derivative_handle(self, form="complex"):
        if form == "complex":
            def second_derivative(tau):
                x_2prime, y_2prime = self.H.second_derivative(tau)
                return x_2prime + y_2prime*1j
        elif form == "real":
            def second_derivative(tau):
                return self.H.second_derivative(tau)
        return second_derivative
    
    def get_third_derivative_handle(self, form="complex"):
        if form == "complex":
            def fun(tau):
                x_3prime, y_3prime = self.H.third_derivative(tau)
                return x_3prime + y_3prime*1j
        elif form == "real":
            def fun(tau):
                return self.H.third_derivative(tau)
        return fun

    def derivative(self, tau):
        fun = self.get_derivative_handle(form="real")
        return fun(tau)

    def second_derivative(self, tau):
        fun = self.get_second_derivative_handle(form="real")
        return fun(tau)
    

class RefCurve(PlanarPHCurve):
    def __call__(self, tau, **kwargs):
        return super().__call__(tau, form="real")
    
