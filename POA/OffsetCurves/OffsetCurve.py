# Author: Suwon Lee from Kookmin University
# Date: April 13, 2023 3:21:36 PM GMT+9
# Code: Offset curve 클래스를 개발한다.

from POA.utilities import map_listed
import numpy as np

class DLBComputer():
    verbose = True
    """
    SamplePoints 클래스와 Obstacle클래스가 주어져있을 때, 샘플지점과 장애물의 관계를 파악하는 클래스.
    """
    def __init__(self, Reference=None, SamplePointsObj=None, step_length=0.01):
        self.Reference = Reference
        self.SamplePoints = SamplePointsObj
        self.step_length = step_length
        self.stepThreshold = 100

    def get_obstacle_index(self, ObstacleObj):
        """
        주어진 샘플지점들 중 장애물 내부에 속하는 샘플지점의 인덱스를 반환하는 함수.
        """
        points = self.SamplePoints.points
        whether_inside = []
        for pt in points:
            position = self.Reference(pt,form="real")
            whether_inside.append(ObstacleObj.check_inside(position))
        # print(whether_inside)
        return [i for i in range(len(whether_inside)) if whether_inside[i] is True]

    def check_nstep_point(self, ObstacleObj, sample_index, n):
        """
        주어진 샘플지점에 대해 n스텝 떨어진 지점이 여전히 장애물 내부에 속하는지 판단하는 함수
        """
        tau = self.SamplePoints.points[sample_index]
        ref_pos = np.reshape(self.Reference(tau,form="real"),newshape=[2,1])
        offset_vector = np.array([[0,-1],[1,0]])@ np.reshape(self.Reference.tangent(tau, form="real"), newshape=[2,1])
        nstep_pos = ref_pos + n*self.step_length*offset_vector
        return ObstacleObj.check_inside(nstep_pos)
    
    def determine_DLB_for_one_sample(self, ObstacleObj, sample_index):
        """
        주어진 샘플지점에 대해 한스텝씩 양방향으로 증가시켜가며 장애물에서 벗어나는 지점의 스텝수를 찾는 함수.
        최대 `stepThreshold` 스텝까지 검사하고, 찾지 못하면 포기한다.
        """
        step = 0
        while step < self.stepThreshold:
            step += 1
            still_inside_positive = self.check_nstep_point(ObstacleObj, sample_index, step)
            still_inside_negative = self.check_nstep_point(ObstacleObj, sample_index, -step)
            if not still_inside_positive:
                return step
            elif not still_inside_negative:
                return -step
        else:
            print("stepThreshold reached. Cannot avoid obstacle.")
            return 0

    def determine_DLBs(self, ObstacleObj):
        """
        모든 샘플지점에 대해, 주어진 장애물 내부에 속하지 않게 하는 DLB를 계산하는 함수.
        """
        obs_index = self.get_obstacle_index(ObstacleObj)
        steps = []
        DLBs = np.zeros_like(self.SamplePoints.points)
        if not obs_index: # empty list is False
            pass
        else:
            for id in obs_index:
                steps.append(self.determine_DLB_for_one_sample(ObstacleObj, id))
                DLBs[id] = steps[-1]*self.step_length
        return DLBs


class DLBComputer_multipleObstacle(DLBComputer):
    def get_obstacle_index(self, ObstacleObjList):
        index_set = set()
        for ObstacleObj in ObstacleObjList:
            index_set.update(set(super().get_obstacle_index(ObstacleObj)))
        return list(index_set)
    
    def determine_DLB_for_one_sample(self, ObstacleObjList, sample_index):
        step = 0
        stepThreshold = 100
        while step < stepThreshold:
            step += 1
            still_inside_positive, still_inside_negative = [], []
            for ObstacleObj in ObstacleObjList:
                still_inside_positive.append(self.check_nstep_point(ObstacleObj, sample_index, step))
                still_inside_negative.append(self.check_nstep_point(ObstacleObj, sample_index, -step))
            if not (True in still_inside_positive):
                return step
            elif not (True in still_inside_negative):
                return -step
        else:
            print("stepThreshold reached. Cannot avoid obstacle.")
            return 0
        
    def determine_DLBs(self, ObstacleObjList):
        return super().determine_DLBs(ObstacleObjList)


class SamplePoints():
    verbose = True
    """
    매개변수화된 곡선에 대한 오프셋 곡선을 정의하기 위해 필요한 샘플링포인트 클래스.
    매개변수화된 곡선상의 점을 샘플링하고, 샘플포인트들에 대해 최소거리를 지정해주면 
    이 최소거리를 만족하도록 하는 오프셋곡선을 정의할 수 있다.
    `SamplePoints` 클래스에서는 이를 위한 인터페이스 등을 제공한다.
    """
    def __init__(self, points=np.linspace(0,1,21)):
        self.points = points  # sampling points
        self.DLB = None     # Distance Lower Bound

    def uniform_sampling(self, num_samples):
        """
        매개변수화된 곡선의 매개변수[0,1]에 대해 균일한 지점을 샘플링하는 함수.
        인스턴스 속성에 저장한다.
        """
        self.points = np.linspace(0,1,num_samples)
        if SamplePoints.verbose:
            print(f"Sample points are updated using uniform sampling: {num_samples} sample points saved.")
        
    def set_constant_distance(self, dist):
        """
        전체 샘플포인트에 대해 상수값으로 오프셋을 준다.
        이 함수는 일반적인 용도로 사용되기보다는 테스트용도로 사용될 것으로 상정하고 개발되었다.
        """
        self.DLB = np.ones_like(self.points)*dist
        if SamplePoints.verbose:
            print(f"Constant distance(={dist}) is allocated for offset curve.")


class ExponentialRBF():
    def __init__(self, epsilon=1):
        self.epsilon = epsilon

    def __call__(self,r):
        return np.exp(-(self.epsilon*r)**2)
    
    def derivative(self, r):
        return -2*self.epsilon**2*r*self.__call__(r)


class OffsetCurve():
    verbose = True
    def __init__(self, Reference=None, SamplePoints=None, RBF=None):
        self.__Reference = Reference        # Reference curve object
        self.__SamplePoints = SamplePoints  # SamplePoints object
        self.__RBF = RBF                    # RBF object
        if self.__RBF is not None:
            self.RBFsum_function = self.get_RBFsum_function()
            self.RBF_drv_sum_function = self.get_RBFsum_function(calc_derivative=True)
            if OffsetCurve.verbose:
                print("RBFsum is ready for the OffsetCurve object.")

    def __call__(self, tau):
        pos_fun = self.get_position_function()
        return np.transpose(np.reshape(pos_fun(tau), newshape=[-1,2]))

    def derivative(self, tau):
        drv_fun = self.get_derivative_function()
        return np.transpose(np.reshape(drv_fun(tau), newshape=[-1,2]))

    def get_modified_RBF_function(self, mu, H, calc_derivative=False):
        if calc_derivative == True:
            function = lambda r:self.__RBF.derivative(r)
        else:
            function = lambda r:self.__RBF(r)
        def fun(tau):
            r = tau - mu
            return H*function(r)
        return fun
    
    def get_RBFsum_function(self, calc_derivative=False):
        mu_list = self.__SamplePoints.points
        H_list = self.__SamplePoints.DLB
        RBF_list = []
        for mu, H in zip(mu_list, H_list):
            RBF_list.append(self.get_modified_RBF_function(mu, H, calc_derivative=calc_derivative))
        def fun(tau):
            eval = 0
            for RBF in RBF_list:
                eval += RBF(tau)
            return eval
        return fun
    
    def get_position_function(self):
        """
        레퍼런스가 2차원 곡선인 경우, 오프셋 곡선을 계산하는 함수이다. 
        """
        @map_listed
        def fun(tau):
            ref_pos = np.reshape(self.__Reference(tau,form="real"),newshape=[2,1])
            phi_val = self.RBFsum_function(tau)
            offset_vector = np.array([[0,-1],[1,0]])@ np.reshape(self.__Reference.tangent(tau, form="real")/np.linalg.norm(self.__Reference.tangent(tau,form="real")), newshape=[2,1])
            return np.reshape(ref_pos + offset_vector*phi_val, newshape=[1,2])
        return fun

    def get_derivative_function(self):
        """
        오프셋 곡선의 미분값을 계산하는 함수이다.
        """
        @map_listed
        def fun(tau):
            rot_mtx = np.array([[0,-1],[1,0]])
            ref_drv = np.reshape(self.__Reference.derivative(tau),newshape=[2,1])
            ref_2nd_drv = np.reshape(self.__Reference.second_derivative(tau),newshape=[2,1])
            phi_val = self.RBFsum_function(tau)
            phi_drv_val = self.RBF_drv_sum_function(tau)
            off_vec = rot_mtx@ np.reshape(self.__Reference.tangent(tau, form="real")/np.linalg.norm(self.__Reference.tangent(tau,form="real")), newshape=[2,1])
            off_vec_drv = rot_mtx@(ref_2nd_drv/np.linalg.norm(ref_drv)-0.5*ref_drv/np.linalg.norm(ref_drv)**3)
            off_pos_drv = ref_drv + phi_drv_val*off_vec + phi_val*off_vec_drv
            return np.reshape(off_pos_drv,newshape=[1,2])
        return fun
    
    @property 
    def Reference(self):
        return self.__Reference
    @Reference.setter 
    def Reference(self, curveObj):
        self.__Reference = curveObj

    @property 
    def SamplePoints(self):
        return self.__SamplePoints
    @SamplePoints.setter 
    def SamplePoints(self, SpObj):
        self.__SamplePoints = SpObj
        self.RBFsum_function = self.get_RBFsum_function()
        self.RBF_drv_sum_function = self.get_RBFsum_function(calc_derivative=True)
    
    @property 
    def RBF(self):
        return self.__RBF
    @RBF.setter 
    def RBF(self, RBFObj):
        self.__RBF = RBFObj