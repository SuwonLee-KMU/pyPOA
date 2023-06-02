# Author: Suwon Lee from Kookmin University
# Date: June 2, 2023 8:23:11 PM GMT+9
# Code: 곡선 근처의 무작위 장애물을 샘플링하는 모듈

import pyPOA.POA.Obstacles as obs
import numpy as np

# 기준곡선 근처의 무작위 점 샘플링
def sample_random_scalar(lb=0.3, ub=0.7):
    """
    기준곡선 근처의 무작위 점을 샘플링하기 위해 스칼라(예:곡선의 매개변수)를 샘플링한다.
    본 함수에서는 [lb, ub]사이의 값을 하나 무작위로 얻는다.
    """
    lower_bound = lb
    upper_bound = ub
    random_tau = np.random.random(size=[])*(upper_bound - lower_bound) + lower_bound
    return random_tau

def sample_random_pos_near(curveObj, sigma=1):
    """
    `curveObj` 곡선 근처의 무작위 점을 샘플링하는 함수.
    `sigma`값을 전달하여 곡선에서 얼마나 떨어진 점을 샘플링할지 결정할 수 있다.
    """
    rnd_tau = sample_random_scalar(lb=0.3, ub=0.7)
    random_vec = list(np.random.randn(2)*sigma)
    sample_pos = [p + r for p, r in zip(curveObj(rnd_tau), random_vec)]
    return sample_pos

def sample_random_obstacle_near(curveObj, num_obs=1, sigma=1, rbound=[0.05,0.3]):
    """
    `cureObj`곡선 근처의 무작위 위치에 무작위 반경을 갖는 원형 장애물을 샘플링하는 함수.
    `sigma`값을 전달하여 곡선에서 얼마나 떨어진 장애물을 샘플링할지 결정할 수 있다.
    `rbound`값을 전달하여 반경의 범위를 지정할 수 있다.
    """
    cobs = []
    cnt = 0
    while len(cobs) < num_obs:
        rnd_radius = sample_random_scalar(lb=rbound[0], ub=rbound[1])
        sample_pos = sample_random_pos_near(curveObj, sigma=sigma)
        cobs.append(obs.CircularObstacle(center=sample_pos, radius=rnd_radius))
        # 다수 장애물 간 겹침이 없도록 하는 기능 작동하지 않는 버그 있음. 추후 검토 요망 #
        mcc = obs.MultipleCollisionChecker(cobs, margin=0)
        evaluations, _ = mcc.check()
        if False in evaluations:
            cobs.pop()
            print("Invalid overlapping obstacles..")
        else:
            print("An obstacle is generated. ", end="\n")
        cnt += 1
        if cnt > 1000:
            print("\nCannot find any place to locate another obstacle. Stop generating...")
            break
        # ----------------------------------------------------------------- #
    return cobs