# Author: Suwon Lee from Kookmin University
# Date: May 1, 2023 10:54:21 AM GMT+9
# Code: 두 개 이상의 벡터필드(속도벡터) 커맨드를 하나로 융합하는 클래스.

from math import sqrt
import numpy as np

class VFGFusion():
    """
    인스턴스속성 vel_cmd_list와 force_cmd_list 값에 근거하여 속도벡터명령을 계산하는 클래스.
    update_cmd, add_cmd 메소드로 내부속성을 업데이트하고,
    fusion 메소드로 속도벡터명령을 출력한다.
    """
    def __init__(self):
        self.vel_cmd_list = []
        self.force_cmd_list = []
    
    def __call__(self):
        pass
                 
    def update_cmd(self, vel_cmd_list:list, force_cmd_list:list):
        self.vel_cmd_list = vel_cmd_list
        self.force_cmd_list = force_cmd_list

    def add_cmd(self, vel_cmd:list, force:float):
        self.vel_cmd_list.append(vel_cmd)
        self.force_cmd_list.append(force)

    def attach_cmd(self, attach_vel_cmd_list:list, attach_force_list:list):
        self.vel_cmd_list += attach_vel_cmd_list
        self.force_cmd_list += attach_force_list

    def fusion(self) -> list:
        cmd = [0,0]
        for vel_cmd, force in zip(self.vel_cmd_list, self.force_cmd_list):
            cmd[0] += vel_cmd[0]*force
            cmd[1] += vel_cmd[1]*force
        normalized_cmd = [c/sqrt(cmd[0]**2+cmd[1]**2) for c in cmd]
        return normalized_cmd
    
    @staticmethod
    def compute_latax_cmd(vel_cmd, vel, gain):
        error = np.arctan2(vel_cmd[1]*vel[0]-vel_cmd[0]*vel[1],vel_cmd[0]*vel[0]+vel_cmd[1]*vel[1]) # angle between 2D vectors
        latax_cmd = error*gain
        return latax_cmd
    