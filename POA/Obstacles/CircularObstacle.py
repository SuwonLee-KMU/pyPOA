# Generated on : Apr 18, 2023 6:02:49 PM
# Author : Suwon Lee from Kookmin University

import numpy as np

class CircularObstacle():
    """
    2차원 평면상에 존재하는 원형 장애물 클래스.
    중심의 위치와 반지름을 지정할 수 있다.
    """
    def __init__(self, center=np.zeros([2,1]), radius=1, **kwargs):
        if "dimension" in kwargs:
            newshape = [kwargs["dimension"],1]
        else:
            newshape = [2,1]
        self.center = np.reshape(center, newshape=newshape)
        self.radius = radius

    def __str__(self):
        string = f"CircularObstacle with center={np.reshape(self.center,newshape=[2])}, radius={self.radius}"
        return string

    def check_inside(self, position, **kwargs):
        """
        어떤 2차원 위치가 주어지면, 해당 위치가 장애물 내부에 있는지 여부를 판단한다.
        kwargs 에 dimension 값이 주어지면 해당 차원으로 적용된다.
        """
        if "dimension" in kwargs:
            newshape = [kwargs["dimension"],1]
        else:
            newshape = [2,1]
        position = np.reshape(position, newshape=newshape)
        distance = np.linalg.norm(self.center - position)
        if distance <= self.radius:
            return True
        else:
            return False
        
    def add_patch(self, fig):
        x0 = self.center[0,0] - self.radius
        y0 = self.center[1,0] - self.radius
        x1 = self.center[0,0] + self.radius
        y1 = self.center[1,0] + self.radius
        fig.add_shape(type="circle",
            xref="x", yref="y",
            fillcolor="PaleTurquoise",
            x0=x0, y0=y0, x1=x1, y1=y1,
            line_color="LightSeaGreen",
            opacity=0.3
        )

class CollisionChecker():
    """
    두 원형 장애물 간 충돌이 있는지 여부를 판별하는 클래스.
    """
    def __init__(self, Obs1, Obs2, margin=0):
        self.Obs1 = Obs1
        self.Obs2 = Obs2
        self.margin = margin

    def check(self):
        """
        중심 간 거리와 반경의 합을 비교하여 안전하면 True, 충돌이 있다면 False를 출력하는 함수.
        """
        center_dist = np.linalg.norm(self.Obs1.center - self.Obs2.center)
        is_safe = center_dist > (self.Obs1.radius + self.Obs2.radius + self.margin)
        return is_safe
    
class MultipleCollisionChecker():
    """
    다수의 원형 장애물 간 충돌이 있는지 여부를 판별하는 클래스.
    """
    def __init__(self, Obs_list, margin=0):
        self.Obs_list = Obs_list
        self.margin = margin

    def check(self):
        num_obs = len(self.Obs_list)
        checker_list = []
        evaluations = []
        for i, obs in enumerate(self.Obs_list):
            try:
                for j in range(i+1, num_obs):
                    checker_list.append(CollisionChecker(obs, self.Obs_list[j], margin=self.margin))
                    evaluations.append(checker_list[-1].check())
            except:
                pass
        return evaluations, checker_list