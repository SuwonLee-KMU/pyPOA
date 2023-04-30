# Author: Suwon Lee from Kookmin University
# Date: April 30, 2023 9:56:39 PM GMT+9
# Code: Repulsive VF wrt obstacle.

from math import sqrt
import pyPOA.PlanarVehicle as pv
import pyPOA.POA as POA

class ObstacleRepulsiveVFG():
    def __init__(self, States:pv.States, rep_offset=0.1):
        self.obs_list = []
        self.is_inside = []
        self.dist_from_surface = []
        self.dir_vec_list = []
        self.rep_force_list = []
        self.__rep_offset = rep_offset
        self.__States = States
        
    def check_inside(self):
        """
        현재위치가 장애물의 영향권내인지 여부를 검사한다.
        """
        x, y = self.__States.x, self.__States.y
        is_inside = []
        dist_from_surface = []
        for obs in self.obs_list:
            cx, cy = obs.center
            dist_from_surface.append(sqrt((cx-x)**2+(cy-y)**2)-obs.radius)
            if dist_from_surface[-1] > self.__rep_offset:
                is_inside.append(False)
            else:
                is_inside.append(True)
        return is_inside, dist_from_surface
    
    def compute_vel_cmd(self):
        """
        장애물 정보를 바탕으로 현재위치에서의 장애물에 의한 velocity command를 계산한다.
        장애물에서 일정수준이하로 가까운 위치에서만 vel_cmd가 발생한다.
        vel_cmd는 dir_vec(방향), force(크기)를 별도로 계산해 준다. 
        """
        x, y = self.__States.x, self.__States.y
        stacked_force = []
        stacked_dir_vec = []
        for obs, dfs, inside in zip(self.obs_list, self.dist_from_surface, self.is_inside):
            if inside:
                cx, cy = obs.center
                force = 1/dfs
                dir_vec = [x-cx, y-cy] / sqrt((cx-x)**2+(cy-y)**2)
            else:
                force = 0
                dir_vec = [0,0]
            stacked_force.append(force)
            stacked_dir_vec.append(dir_vec)
        return stacked_dir_vec, stacked_force

    @property 
    def States(self):
        return self.__States
    @States.setter 
    def States(self, States:pv.States):
        self.__States = States
        self.is_inside, self.dist_from_surface = self.check_inside()
        self.dir_vec_list, self.rep_force_list = self.compute_vel_cmd()
    @property 
    def rep_offset(self):
        return self.__rep_offset
    @rep_offset.setter 
    def rep_offset(self, rep_offset:float):
        self.__rep_offset = rep_offset
        self.is_inside, self.dist_from_surface = self.check_inside()
        self.dir_vec_list, self.rep_force_list = self.compute_vel_cmd()