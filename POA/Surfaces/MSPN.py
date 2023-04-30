# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:34:51 PM GMT+9
# Code: 

import numpy as np
import sympy
from POA.Polynomials import BinomialPolynomial
from POA.utilities import map_listed

class MSPN(): # Minimum Surface with Pythagorean Normal
    def __init__(self, preimage1, preimage2, preimage3, r0 = [0,0,0]):
        self.u  = preimage1
        self.v  = preimage2
        self.w  = preimage3
        self.r0 = np.array(r0)

    def __call__(self, z):
        func = self.get_surface_handle(form='complex')
        return func(z)

    def get_holomorphic_handles(self):
        """
        u, v, w 선형복소함수에 대해 Enneper–Weierstrass parameterization을 수행하여, 복소함수를 출력하는 함수.
        """
        u, v, w = self.u, self.v, self.w
        def Phi1(z):
            return w(z)*(u(z)**2+v(z)**2)
        def Phi2(z):
            return 2*w(z)*u(z)*v(z)
        def Phi3(z):
            return 1j*w(z)*(u(z)**2+v(z)**2)
        return Phi1, Phi2, Phi3

    def get_holomorphic_coefs(self):
        """
        다항식의 계수 u, v, w에 대해 Enneper–Weierstrass parameterization을 수행하여, 복소함수의 계수를 출력하는 함수.
        계수 리스트는 높은차수부터 낮은차수 순서로 쓴다.
        """
        u = np.flip(self.u.coefficients,axis=0)
        v = np.flip(self.v.coefficients,axis=0)
        w = np.flip(self.w.coefficients,axis=0)
        phi1 = np.polymul(w, np.polysub(np.polymul(u,u), np.polymul(v,v)))
        phi2 = np.polymul(2, np.polymul(w, np.polymul(u, v)))
        phi3 = np.polymul(1j, np.polymul(w, np.polyadd(np.polymul(u,u), np.polymul(v,v))))
        return phi1, phi2, phi3

    @staticmethod
    def poly_complex_integration(poly,c0=0):
        """
        복소함수(계수로 주어짐)를 복소평면상의 0에서 해당 복소수까지 선적분(부정적분)하여 얻은 복소함수의 계수를 출력하는 함수.
        """
        n = len(list(poly))
        aux = np.array([1/n for n in np.arange(start=n,stop=0,step=-1)])
        out = aux*np.array(poly)
        return np.append(out,c0)

    def get_surface_complex_poly_coef(self, mode='raw'): 
        """
        Enneper-Weierstrass parameterization을 통해 얻은 함수로부터 평면을 나타내는 복소함수의 계수를 출력하는 함수.
        높은차수부터 낮은차수 순서로 출력함. 논문의 결과와 비교해서 검증완료.
        """
        phi1_poly, phi2_poly, phi3_poly = self.get_holomorphic_coefs() 
        rx_poly = MSPN.poly_complex_integration(phi1_poly,self.r0[0])
        ry_poly = MSPN.poly_complex_integration(phi2_poly,self.r0[1])
        rz_poly = MSPN.poly_complex_integration(phi3_poly,self.r0[2])
        if mode == "raw":
            return rx_poly, ry_poly, rz_poly
        elif mode == "numpy":
            r_coef_raw = (rx_poly, ry_poly, rz_poly)
            mlen = max([len(x) for x in r_coef_raw])
            r_coef = []
            for i in range(3):
                num_zeros = mlen - len(r_coef_raw[i])
                r_coef.append(np.concatenate([np.zeros(num_zeros), r_coef_raw[i]],axis=0))
            r_coef = np.array(r_coef)
            return r_coef
        else:
            raise ValueError(f'invalid mode argument: {mode}')

    @staticmethod
    def poly_coef_to_derivative_poly_coef(poly):
        """
        다항식의 계수가 주어지면, 해당 다항식을 미분하여 얻은 다향식의 계수를 출력하는 함수.
        입출력 계수 모두, 높은 차수부터 낮은 차수 순서로 쓴다.
        """
        n = len(list(poly))-1
        coef = list(range(n,-1,-1))
        new_coef = [c*r for c, r in zip(coef, poly)][:-1]
        return new_coef        

    def nth_derivative_coef(self, n):
        """
        3차원 평면의 다항식을 n번 미분한 다항식의 계수를 출력하는 함수.
        높은차수부터 낮은차수 순서로 출력한다.
        """
        poly_list = self.get_surface_complex_poly_coef(mode="raw")
        coefs = []
        for poly in poly_list:
            for k in range(n):
                poly = MSPN.poly_coef_to_derivative_poly_coef(poly)
            coefs.append(poly)
        return coefs

    def get_surface_handle(self, form='complex'):
        """
        평면을 나타내는 복소함수의 계수들로부터 평면의 함수를 출력하는 함수.
        """
        rx_poly, ry_poly, rz_poly = self.get_surface_complex_poly_coef()
        if form == 'real':
            def r(x, y):
                z = x + y*1j
                rx = np.real(np.polyval(rx_poly, z)) + self.r0[0]
                ry = np.real(np.polyval(ry_poly, z)) + self.r0[1]
                rz = np.real(np.polyval(rz_poly, z)) + self.r0[2]
                return rx, ry, rz
        elif form == "complex":
            def r(z):
                rx = np.real(np.polyval(rx_poly, z)) + self.r0[0]
                ry = np.real(np.polyval(ry_poly, z)) + self.r0[1]
                rz = np.real(np.polyval(rz_poly, z)) + self.r0[2]
                return rx, ry, rz
        else:
            raise ValueError('Invalid form argument')
        return r

    def nth_derivative_handle(self, n):
        """
        평면의 복소다항식을 n회 복소미분한 함수 핸들.
        """
        coefs = self.nth_derivative_coef(n)
        def nderiv(z):
            derivative = [np.polyval(c, z) for c in coefs]
            return derivative
        return nderiv
        
    def symbolic(self): # Developing
        rx_poly, ry_poly, rz_poly = self.get_surface_complex_poly_coef()
        z = sympy.Symbol('z')
        # rx, ry, rz = sympy.symbols('r_x, r_y, r_z')
        rx = sympy.Poly(rx_poly, z).as_expr() + self.r0[0]
        ry = sympy.Poly(ry_poly, z).as_expr() + self.r0[1]
        rz = sympy.Poly(rz_poly, z).as_expr() + self.r0[2]
        return rx, ry, rz
        
    def normal(self, z):
        fun = self.get_normal_handle(form="complex")
        return fun(z)

    def get_normal_handle(self, form="complex"):
        @map_listed
        def get_normal(z):
            uR = np.real(self.u(z))
            uI = np.imag(self.u(z))
            vR = np.real(self.v(z))
            vI = np.imag(self.v(z))
            wR = np.real(self.w(z))
            wI = np.imag(self.w(z))
            gR = [
                wR*(uR**2 - uI**2 - vR**2 + vI**2) - 2*wI*(uR*uI - vR*vI),
                2*wR*(uR*vR - uI*vI) - 2*wI*(uR*vI + uI*vR),
                -2*wR*(uR*uI + vR*vI) - wI*(uR**2 - uI**2 + vR**2 - vI**2)
            ]
            gI = [
                2*wR*(uR*uI - vR*vI) + wI*(uR**2 - uI**2 - vR**2 + vI**2),
                2*wR*(uR*vI - uI*vR) + 2*wI*(uR*vR - uI*vI),
                wR*(uR**2 - uI**2 + vR**2 - vI**2) - 2*wI*(uR*uI + vR*vI)
            ]
            rx, ry = gR, [-e for e in gI]
            N = normalize(np.cross(rx, ry))
            return N
        
        if form == "real":
            def fun(x, y):
                z = x + y*1j
                if isinstance(z,list):
                    N = [get_normal(v) for v in z]
                else:
                    N = get_normal(z)
                return N
        elif form == "complex":
            def fun(z):
                if isinstance(z,list):
                    N = [get_normal(v) for v in z]
                else:
                    N = get_normal(z)
                return N
            
        return fun

    def binomial_polynomial(self):
        """
        곡면의 함수를 이변수 다항식으로 표현하여 출력하는 함수.
        """
        r_coefs = self.get_surface_complex_poly_coef(mode="numpy")
        r_coefs = [rc[::-1] for rc in r_coefs]
        r_bc = [] #dimension은 u, v, w의 최대차수보다 2배+1 이상으로 설정해야 함.
        for i, rc in enumerate(r_coefs):
            bp = BinomialPolynomial(dimension=len(r_coefs[i]))
            bp.update_from_coef_vector(rc, part="real")
            r_bc.append(bp)
        return r_bc
