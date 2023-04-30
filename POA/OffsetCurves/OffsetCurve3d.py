# Generated on : Apr 18, 2023 6:00:59 PM
# Author : Suwon Lee from Kookmin University

import numpy as np
from POA.OffsetCurves import OffsetCurve
from POA.OffsetCurves import SamplePoints
from POA.utilities import map_listed

class SamplePoints3d(SamplePoints):
    verbose = True
    """
    매개변수화된 곡선에 대한 오프셋 곡선을 정의하기 위해 필요한 샘플링포인트 클래스.
    매개변수화된 곡선상의 점을 샘플링하고, 샘플포인트들에 대해 최소거리를 지정해주면 
    이 최소거리를 만족하도록 하는 오프셋곡선을 정의할 수 있다.
    `SamplePoints3d` 클래스에서는 이를 위한 인터페이스 등을 제공한다.
    최소거리는 법선벡터와 종법선벡터의 선형결합 벡터 방향을 갖는다.
    두 벡터의 선형결합 계수가 본 클래스의 인스턴스 속성으로 들어간다.
    계수를 제곱해서 더했을 때 1이 되어야 하므로, 인스턴스속성은 법선벡터계수만으로 충분하다.
    (종법선벡터계수는 법선벡터계수가 결정되면 그에 따라 계산된다)
    """
    def __init__(self, points=np.linspace(0,1,21)):
        self.points = points  # sampling points
        self.DLB = None     # Distance Lower Bound
        self.coef_normal = None # coefficient of the normal vector

    @property 
    def coef_binormal(self):
        return np.array([np.sqrt(1-an**2) for an in self.coef_normal])

class OffsetCurve3d(OffsetCurve):
    verbose = True
    def __init__(self, Reference=None, SamplePoints=None, RBF=None):
        super().__init__(Reference, SamplePoints, RBF)

    def __call__(self, tau):
        pos_fun = self.get_position_function()
        return np.transpose(np.reshape(pos_fun(tau), newshape=[-1,3]))
    
    def get_position_function(self):
        """
        레퍼런스가 3차원 곡선인 경우, 오프셋 곡선을 계산하는 함수이다.
        """
        @map_listed
        def fun(tau):
            ref_pos = np.reshape(self.__Reference(tau,form="real"),newshape=[3,1])
            phi_val = self.RBFsum_function(tau)
            offset_vector = np.array([[0,-1],[1,0]])@ np.reshape(self.__Reference.tangent(tau, form="real"), newshape=[2,1])
            
            return np.reshape(ref_pos + offset_vector*phi_val, newshape=[1,3])
        return fun