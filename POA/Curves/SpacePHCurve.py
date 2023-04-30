# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:15:00 PM GMT+9
# Code: 

import numpy as np
from POA.utilities import map_listed
import sympy

class SpacePHCurve():
    def __init__(self, MSPN, planarPHcurve):
        self.r = MSPN
        self.q = planarPHcurve

    def __call__(self, tau):
        func = self.get_curve_handle()
        return func(tau)

    def get_curve_handle(self):
        planar_curve_fun = self.q.get_curve_handle(form="complex")  
        surface_fun      = self.r.get_surface_handle(form="complex")
        def curve(tau):
            q_eval = planar_curve_fun(tau)
            r_eval = surface_fun(q_eval)
            return r_eval
        return curve

    def get_curve_poly_coef(self):
        q_coef = np.flip(self.q.get_curve_coef(), axis=1)                               # 높은차수부터 낮은차수 순서로 출력되므로 뒤집어 줌
        r_coef = np.flip(self.r.get_surface_complex_poly_coef(mode="numpy"), axis=1)    # 높은차수부터 낮은차수 순서로 출력되므로 뒤집어 줌
        return SpacePHCurve.compute_coef_of_composite(q_coef, r_coef)

    def get_normal_handle(self):
        planar_curve_fun = self.q.get_curve_handle(form="complex")  
        surface_fun      = self.r.get_normal_handle(form="complex")
        def normal(tau):
            q_eval = planar_curve_fun(tau)
            normal = surface_fun(q_eval)
            return normal
        return normal
    
    def tangent_handle(self):
        """
        공간곡선의 일차 미분을 계산하는 함수 핸들을 출력하는 함수.
        공간곡선은 평면과 평면곡선의 합성으로 정의되므로, 
        공간곡선의 미분은 평면(복소함수)의 복소미분과 
        평면곡선의 매개변수에 대한 미분의 연쇄법칙으로 계산된다.
        """
        qprime_fun = self.q.get_derivative_handle(form="complex")
        rprime_fun = self.r.nth_derivative_handle(n=1)
        def tangent(tau):
            q = self.q(tau,form="complex")
            qprime = qprime_fun(tau)
            rprime = rprime_fun(q)
            tangent = [np.real(rp*qprime) for rp in rprime]
            return tangent
        return tangent

    def binomial_tangent_handle(self):
        qprime_fun = self.q.get_derivative_handle(form="real")
        r_bp = self.r.binomial_polynomial()
        @map_listed
        def tangent(tau):
            q = self.q(tau, form="real")
            qprime = qprime_fun(tau)
            tgt = []
            for r_k in r_bp:
                Dxr = r_k.derivative(variable="x")(*q)
                Dyr = r_k.derivative(variable="y")(*q)
                tgt.append(Dxr*qprime[0] + Dyr*qprime[1])
            return tgt
        return tangent

    def binomial_second_derivative_handle(self):
        qprime_fun = self.q.get_derivative_handle(form="real")
        q2prime_fun = self.q.get_second_derivative_handle(form="real")
        r_bp = self.r.binomial_polynomial()
        @map_listed
        def second_derivative(tau):
            q = self.q(tau, form="real")
            qprime = qprime_fun(tau)
            q2prime = q2prime_fun(tau)
            second_deriv = []
            for r_k in r_bp:
                Dxr = r_k.derivative(variable="x")(*q)
                Dyr = r_k.derivative(variable="y")(*q)
                D2xr = r_k.derivative("x").derivative("x")(*q)
                D2yr = r_k.derivative("y").derivative("y")(*q)
                DxDyr = r_k.derivative("x").derivative("y")(*q)
                second_deriv.append(D2xr*qprime[0]**2 + Dxr*q2prime[0] + D2yr*qprime[1]**2 + Dyr*q2prime[1] + 2*DxDyr*qprime[0]*qprime[1])
            return second_deriv
        return second_derivative
    
    def binomial_third_derivative_handle(self):
        qprime_fun = self.q.get_derivative_handle(form="real")
        q2prime_fun = self.q.get_second_derivative_handle(form="real")
        q3prime_fun = self.q.get_third_derivative_handle(form="real")
        r_bp = self.r.binomial_polynomial()
        @map_listed
        def fun(tau):
            q = self.q(tau, form="real")
            qprime = qprime_fun(tau)
            q2prime = q2prime_fun(tau)
            q3prime = q3prime_fun(tau)
            third_deriv = []
            for r_k in r_bp:
                Dxr = r_k.derivative(variable="x")(*q)
                Dyr = r_k.derivative(variable="y")(*q)
                D2xr = r_k.derivative("x").derivative("x")(*q)
                D2yr = r_k.derivative("y").derivative("y")(*q)
                DxDyr = r_k.derivative("x").derivative("y")(*q)
                D3xr = r_k.derivative("x").derivative("x").derivative("x")(*q)
                D3yr = r_k.derivative("y").derivative("y").derivative("y")(*q)
                D2xDyr = r_k.derivative("x").derivative("x").derivative("y")(*q)
                DxD2yr = r_k.derivative("x").derivative("y").derivative("y")(*q)
                third_deriv.append(
                    (D3xr*qprime[0]+D2xDyr*qprime[1])*qprime[0]**2 + D2xr*2*qprime[0]*q2prime[0] + \
                    (D2xr*qprime[0]+DxDyr*qprime[1])*q2prime[0] + Dxr*q3prime[0] + \
                    2*(D2xDyr*qprime[0]+DxD2yr*qprime[1])*qprime[0]*qprime[1] + 2*DxDyr*(q2prime[0]*qprime[1]+qprime[0]*q2prime[1]) + \
                    (DxD2yr*qprime[0]+D3yr*qprime[1])*qprime[1]**2 + D2yr*2*qprime[1]*q2prime[1] + \
                    (DxDyr*qprime[0] + D2yr*qprime[1])*q2prime[1] + Dyr*q3prime[1]
                )
            return third_deriv
        return fun
    
    def derivative(self, tau):
        fun = self.binomial_tangent_handle()
        return fun(tau)
    
    def second_derivative(self, tau):
        fun = self.binomial_second_derivative_handle()
        return fun(tau)

    def third_derivative(self, tau):
        fun = self.binomial_third_derivative_handle()
        return fun(tau)

    def second_derivative_handle(self):
        """
        공간곡선의 2차 미분을 계산하는 함수 핸들을 출력하는 함수.
        연쇄법칙으로 계산된다.
        """
        rprime_fun = self.r.nth_derivative_handle(n=1)
        r2prime_fun = self.r.nth_derivative_handle(n=2)
        qprime_fun = self.q.get_derivative_handle(form="complex")
        q2prime_fun = self.q.get_second_derivative_handle(form="complex")
        def second_derivative(tau):
            q = self.q(tau,form="complex")
            qprime = qprime_fun(tau)
            q2prime = q2prime_fun(tau)
            rprime = rprime_fun(q)
            r2prime = r2prime_fun(q)
            second_deriv = np.real([r2p*qprime**2 + rp*q2prime for rp, r2p in zip(rprime, r2prime)])
            return second_deriv
        return second_derivative

    def get_tangent_handle_numerical(self):
        def tangent_fun(tau):
            eps = 0.001
            p1 = np.transpose(self.__call__(tau+eps))
            p0 = np.transpose(self.__call__(tau))
            if isinstance(tau, list):
                tan = [normalize(e) for e in p1-p0]
            else:
                tan = normalize(p1-p0)
            return tan
        return tangent_fun

    # 합성함수의 계수를 도출하는 코드 작성. 계수는 낮은 차수부터 높은 차수 순서로 쓴다.
    @staticmethod
    def compute_coef_of_composite(coef_q, coef_r): # 잘못됨. 버그 수정 필요.
        shape_q = np.shape(coef_q) # [2 x n+1]
        shape_r = np.shape(coef_r) # [3 x m+1]
        coef    = np.zeros([3,shape_q[1]+shape_r[1]], dtype=complex)
        for j in range(shape_r[1]):
            for k in range(shape_q[1]):
                power = j+k
                coef[0,power] += coef_r[0,j]*(coef_q[0,k] + coef_q[1,k]*1j)**j
                coef[1,power] += coef_r[1,j]*(coef_q[0,k] + coef_q[1,k]*1j)**j
                coef[2,power] += coef_r[2,j]*(coef_q[0,k] + coef_q[1,k]*1j)**j
        return coef

    def symbolic(self):
        tau = sympy.Symbol('tau', real=True)
        func = self.get_curve_handle()
        return func(tau)
    
    def get_poly_coefs(self):
        """
        심볼릭 표현식으로부터 합성함수 p=(r◦q)의 다항식 표현을 얻는다.
        계산 시간이 오래 걸리므로 최적화에 사용하기에는 부적합하다. 해석적인 방법으로 다시 구하는 중.
        """
        symb = self.symbolic() # 3차원(X,Y,Z) 심볼릭 표현식
        polycoef_complex = [sympy.Poly(sym).all_coeffs() for sym in symb] # 각 차원별 계수 추출
        polycoef = [[sympy.re(coef)for coef in polycoef] for polycoef in polycoef_complex] # 실수부 추출
        return np.array(polycoef).astype(float)

    def get_deriv_coefs(self, how_many):
        """
        심볼릭 표현식으로 얻은 공간곡선의 다항식을 미분하여 얻은 다항식의 계수를 구하는 함수.
        """
        polycoef = self.get_poly_coefs()
        deriv = lambda poly: [poly[i] * i for i in range(1, len(poly))]
        for i in range(how_many):
            polycoef = [deriv(pc) for pc in polycoef]
        return polycoef
    
    def get_tangent_handle(self):
        """
        공간곡선의 접선을 평가하는 함수를 출력하는 함수.
        심볼릭 표현식으로부터 계산하므로, 계산시간이 오래 걸린다.
        """
        deriv1_coef = self.get_deriv_coefs(1)
        def fun(tau):
            return np.transpose([np.polyval(coef, tau) for coef in deriv1_coef])
        return fun
    
    def get_curvature_handle(self):
        """
        공간곡선의 곡률을 계산하는 함수를 출력하는 함수.
        심볼릭 표현식으로부터 계산하므로, 계산시간이 오래 걸린다.
        """
        fun_r_prime = self.get_tangent_handle()
        deriv2_coef = self.get_deriv_coefs(2)
        def fun_r_2prime(tau):
            return np.transpose([np.polyval(coef, tau) for coef in deriv2_coef])
        def fun_curvature(tau):
            num = np.linalg.norm(np.array([np.cross(rp, rpp) for rp, rpp in zip(fun_r_prime(tau), fun_r_2prime(tau))]),axis=1)
            den = np.linalg.norm(fun_r_prime(tau),axis=1)
            return num/den
        return fun_curvature

    def curvature(self):
        """
        공간곡선의 곡률을 계산하는 함수를 출력한다.
        Binomial poly를 사용하여 계산한다.
        """
        fun_r_prime = self.binomial_tangent_handle()
        deriv2_coef = self.get_deriv_coefs(2)
        def fun_r_2prime(tau):
            return np.transpose([np.polyval(coef, tau) for coef in deriv2_coef])
        def fun_curvature(tau):
            num = np.linalg.norm(np.array([np.cross(rp, rpp) for rp, rpp in zip(fun_r_prime(tau), fun_r_2prime(tau))]),axis=1)
            den = np.linalg.norm(fun_r_prime(tau),axis=1)**(3/2)
            return num/den
        return fun_curvature
