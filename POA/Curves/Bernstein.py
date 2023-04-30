# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:18:19 PM GMT+9
# Code: 

import numpy as np
import math
from math import comb

def Bernstein_basis_polynomial(n, nu):
    def b_n_nu(x):
        return comb(n, nu) * x**nu * (1-x)**(n-nu)
    return b_n_nu

def Bernstein_basis_derivative_polynomial(n, nu, p):
    """
    Bernstein 기저 다항식의 미분 함수 핸들을 출력하는 함수. 
    Arguments: 
        n : 미분할 다항식의 차수
        nu : 기저 다항식의 인덱스
        p : 미분 횟수
    """
    def Dpb_n_nu(x):
        if p <= n:
            coef = math.factorial(n)/math.factorial(n-p)
            k_list = range(max(0, nu+p-n),min(nu,p)+1)
            return sum([
                coef 
                * (-1)**(k+p)*comb(p, k) 
                * Bernstein_basis_polynomial(n-p,nu-k)(x) 
                for k in k_list
                ])
        else:
            return np.zeros_like(x)
    return Dpb_n_nu

# Bernstein Polynomial
def Bernstein_polynomial(betas):
    n = np.size(betas) - 1 # n: degree of the polynomial
    def Bn(x):
        basis = np.array([Bernstein_basis_polynomial(n,nu)(x) for nu in range(n+1)])
        return np.dot(betas, basis)
    return Bn

def Bernstein_derivative_polynomial(betas, p):
    """
    Bernstein 다항식을 미분한 함수 핸들을 출력한다.
    Arguments:
        betas : 번스타인 다항식의 계수
        p : 미분한 횟수
    """
    n = np.size(betas) - 1
    def DpBn(x):
        basis = np.array([Bernstein_basis_derivative_polynomial(n,nu,p)(x) for nu in range(n+1)])
        return np.dot(betas, basis)
    return DpBn

def Bernsteinbasis_to_power_coef(n, k):
    coef = []
    for i in range(0,n+1):
        if i < k:
            coef.append(0)
        else:
            coef.append((-1)**(i-k)*comb(n,i)*comb(i,k))
    return coef

def BernsteinPoly_to_power_coef(CPs): # 낮은차수의 계수부터 높은차수의 계수 순서로 출력한다.
    shape = np.shape(CPs) # [n-1, dim]
    n   = shape[0] - 1
    dim = shape[1]
    coef = None
    for k in range(n+1):
        power_coef = np.array(Bernsteinbasis_to_power_coef(n, k))
        if coef is None:
            coef = np.array([power_coef*cp for cp in CPs[k,:]])
        else:
            coef = coef + np.array([power_coef*cp for cp in CPs[k,:]])
    return coef

