# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 11:48:26 PM GMT+9
# Code: 

from pyPOA.POA.Curves import Dubins as db
import numpy as np
from scipy import io

class SpaceDubinsCurves():
    def __init__(self, Ps, Ts, turn_radius, delta=0.001):
        self.Ps = Ps
        self.Ts = Ts
        self.turn_radius = turn_radius
        self.delta = delta

    def get_paths(self):
        # User's waypoints: [x, y, heading (degrees)]
        Wptz = []
        for i in range(len(self.Ps)):
            angle = np.arctan2(self.Ts[i][0], self.Ts[i][1])*180/np.pi
            Wptz.append(db.Waypoint(self.Ps[i][0], self.Ps[i][1], angle))    
        
        paths = []
        for i in range(len(Wptz)-1):
            param = db.calcDubinsPath(Wptz[i], Wptz[i+1], turn_radius=0.2)
            path = db.dubins_traj(param,.001)
            paths.append(path)

        # path 리스트에 고도 정보 추가
        for i, path in enumerate(paths):
            npoints  = len(path)
            tmp      = np.reshape(np.linspace(self.Ps[i][2], self.Ps[i+1][2], npoints),newshape=[-1,1])
            paths[i] = np.concatenate([path,tmp], axis=1)

        return paths

    def save_matfile(self, fname):
        paths = self.get_paths()
        data = dict()
        curves = []
        for i, path in enumerate(paths):
            x, y, phi, z = np.transpose(path)
            space_curve = dict()
            space_curve['index'] = i
            space_curve['sample_x'] = x
            space_curve['sample_y'] = y
            space_curve['sample_z'] = z
            curves.append(space_curve)
        data['curves'] = curves
        io.savemat(fname, data)
        print('File saved!')
