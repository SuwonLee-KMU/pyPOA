# Author: Suwon Lee from Kookmin University
# Date: April 30, 2023 9:56:39 PM GMT+9
# Code: Repulsive VF wrt obstacle.

from math import sqrt
import pyPOA.PlanarVehicle as pv

class ObstacleRepulsiveVFG():
    def __init__(self, States):
        self.obs_list = []
        self.is_inside = []
        self.rep_offset = 0.1
        self.__States = States
        
    def check_inside(self):
        """
        현재위치가 장애물의 영향권내인지 여부를 검사한다.
        """
        x, y = self.States.x, self.States.y
        is_inside = []
        for obs in self.obs_list:
            cx, cy = obs.center
            dist_from_surface = (sqrt((cx-x)**2+(cy-y)**2)-obs.radius)
            if dist_from_surface > self.rep_offset:
                is_inside.append(False)
            else:
                is_inside.append(True)
        return is_inside
    
    def compute_vel_cmd(self):
        x, y = self.States.x, self.States.y

    @property 
    def States(self):
        return self.__States
    @States.setter 
    def States(self, States:pv.States):
        self.__States = States
        self.is_inside = self.check_inside()