# Author: Suwon Lee from Kookmin University
# Date: April 19, 2023 12:01:04 AM GMT+9
# Code: 평가된 곡선 객체를 만든다.

import numpy as np
import plotly.graph_objects as go

class Evaluated():
    def __init__(self, tau=np.linspace(0,1,101), x=np.zeros([101]), y=np.zeros([101]), z=np.zeros([101])):
        self.tau = tau
        self.x = x
        self.y = y
        self.z = z

    def show_path(self, fig=None):
        if fig is None:
            fig = go.Figure()
        fig.add_scatter3d(x=self.x, y=self.y, z=self.z)
        return fig    
    
    def show_x(self, fig=None):
        if fig is None:
            fig = go.Figure()
        fig.add_scatter(x=self.tau, y=self.x)
        return fig    
    
    def show_y(self, fig=None):
        if fig is None:
            fig = go.Figure()
        fig.add_scatter(x=self.tau, y=self.y)
        return fig    
    
    def show_z(self, fig=None):
        if fig is None:
            fig = go.Figure()
        fig.add_scatter(x=self.tau, y=self.z)
        return fig    
