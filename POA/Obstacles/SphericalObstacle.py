# Generated on : Apr 18, 2023 6:05:13 PM
# Author : Suwon Lee from Kookmin University

import numpy as np
from .CircularObstacle import CircularObstacle

class SphericalObstacle(CircularObstacle):
    """ 
    3차원 공간상에 존재하는 구형 장애물 클래스.
    중심의 위치와 반지름을 지정할 수 있다.
    """
    def __init__(self, center=np.zeros([0,0,0]), radius=1):
        super().__init__(center, radius, dimension=3)

    def check_inside(self, position):
        return super().check_inside(position, dimension=3)
    
    def surface(self, fig):
        def makesphere(x, y, z, radius, resolution=10):
            """Return the coordinates for plotting a sphere centered at (x,y,z)"""
            u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
            X = radius * np.cos(u)*np.sin(v) + x
            Y = radius * np.sin(u)*np.sin(v) + y
            Z = radius * np.cos(v) + z
            return (X, Y, Z)

        X, Y, Z = makesphere(x=self.center[0], y=self.center[1], z=self.center[2], radius = self.radius)
        fig.add_surface(x=X,y=Y,z=Z,colorscale="viridis")
        return fig