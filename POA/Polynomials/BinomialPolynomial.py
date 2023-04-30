# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:36:25 PM GMT+9
# Code: 

import numpy as np
from math import comb

class BinomialPolynomial():
    """
    이변수 다항식을 구현하는 클래스.
    이변수 다항식은 x, y의 power 항들의 합으로 구성된다. 다항식의 계수는 행렬 형태로 주어진다.
    행렬을 사용하여 다항식 f(x,y) 의 값을 계산하면 다음과 같이 표현된다. 
    f(x, y) = X.T @ M @ Y
    여기서 X=[1, x, x**2 ..., x**n]이고, Y도 비슷하게 정의된다.
    
    Properties
        _dimension : 계수행렬의 차원이다. 각 변수의 최대차수가 n일 때, _dimension = n+1이다.
        coefficient : 계수행렬
    """
    def __init__(self, dimension:int):
        self._dimension = dimension
        self.coefficient = [[0]*dimension for _ in range(dimension)]

    def __call__(self, x:float, y:float):
        pow = self._dimension - 1
        X = np.reshape([x**k for k in range(pow+1)], newshape=[pow+1,1])
        Y = np.reshape([y**k for k in range(pow+1)], newshape=[pow+1,1])
        return np.reshape(X.T @ np.array(self.coefficient) @ Y, newshape=[])

    def derivative(self, variable="x"):
        """
        함수의 (x, y)점에서의 편미분을 계산하여, BinomialPolynomial 클래스로 반환한다.
        variable : 어떤 변수로 편미분할지 결정한다. "x", "y" 중 하나의 값을 가질 수 있다.
        """
        if variable == "x":
            M = BinomialPolynomial.derivative_matrix(self._dimension).T @ np.array(self.coefficient)
        elif variable == "y":
            M = np.array(self.coefficient) @ BinomialPolynomial.derivative_matrix(self._dimension)
        derivative_BP = BinomialPolynomial(self._dimension)
        derivative_BP.coefficient = M
        return derivative_BP
    
    def update_from_coef_vector(self, coef_vector, part="real"):
        coef_real, coef_imag = self.vector_to_matrix(coef_vector)
        if part == "real":
            self.coefficient = coef_real
        else:
            self.coefficient = coef_imag

    def vector_to_matrix(self, coef_vector):
        """
        벡터(1차원 리스트)형태로 주어진 복소다항식의 계수로부터 이변수 다항식의 계수행렬을 생성하는 함수.
        Arguments
            coef_vector : 복소다항식의 계수. 낮은 차수(0차)부터 높은 차수 순서로 쓴다.
        """
        if (cvl:=len(list(coef_vector))) > self._dimension:
            raise ValueError(f"coef_vector의 길이(={cvl})는 계수행렬의 크기(={self._dimension})보다 작거나 같아야 합니다.")
        
        coef_R = [a.real for a in coef_vector]
        coef_I = [a.imag for a in coef_vector]
        M_R, M_I = [], []
        for k in range(self._dimension):
            M_R_k, M_I_k = BinomialPolynomial.binomial_basis(k,self._dimension)
            M_R_k, M_I_k = np.array(M_R_k), np.array(M_I_k)
            M_R.append(M_R_k)
            M_I.append(M_I_k)
        coef_real = sum([cR*MR - cI*MI for cR, cI, MR, MI in zip(coef_R, coef_I, M_R, M_I)])
        coef_imag = sum([cI*MR + cR*MI for cR, cI, MR, MI in zip(coef_R, coef_I, M_R, M_I)])
        return coef_real, coef_imag
    
    @staticmethod
    def binomial_basis(power, dim):
        if power >= dim:
            raise ValueError(f"power(={power})는 dim(={dim})보다 작아야 합니다.")
        M = [[0]*dim for _ in range(dim)]
        N = [[0]*dim for _ in range(dim)]
        coef_re = []
        coef_im = []
        for k in range(power+1):
            if k%4==0:
                coef_re.append(comb(power,k))
                coef_im.append(0)
            elif k%4==1:
                coef_re.append(0)
                coef_im.append(comb(power,k))
            elif k%4==2:
                coef_re.append(-comb(power,k))
                coef_im.append(0)
            else:
                coef_re.append(0)
                coef_im.append(-comb(power,k))
        rindex = range(power,-1,-1)
        cindex = range(0,power+1,1)
        for r,c,re,im in zip(rindex,cindex,coef_re,coef_im):
            M[r][c] = re
            N[r][c] = im
        return M,N
    
    @staticmethod
    def derivative_matrix(dim:int):
        """
        다항식 f의 계수와 곱해서 f를 미분한 다항식의 계수를 얻을 수 있는 행렬 D를 계산하는 함수.
        f = a0 + a1*x + a2*x**2 + ... + an*x**n = A.T @ X 일 때, f의 미분은 다음과 같다.
        Df = a1 + 2*a2*x + 3*a3*x**2 + ... + n*an-1*x**n-1 + 0*x**n = DA.T @ X
        여기서 DA.T @ X = A.T @ D @ X 가 성립하는 D를 계산한다.
        dim : 정방행렬의 크기. dim = n+1.
        """
        D = [[0]*dim for _ in range(dim)]
        for k in range(1,dim):
            D[k][k-1] = k
        return np.array(D)
