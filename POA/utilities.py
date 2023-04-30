# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:16:15 PM GMT+9
# Code: 


import numpy as np

def map_listed(func):
    """
    입력값이 리스트라면, 각 원소에 대해 함수를 매핑해주는 데코레이터.
    """
    def wrapper(arg):
        if isinstance(arg,list) or isinstance(arg,np.ndarray):
            return [func(a) for a in arg]
        else:
            return func(arg)
    return wrapper    


def normalize(vector):
    if len(np.shape(vector)) == 1:
        return vector / np.linalg.norm(vector)
    else:
        raise ValueError(f'Invalid shape of the input vector : {np.shape(vector)}')
