# Author: Suwon Lee from Kookmin University
# Date: April 18, 2023 7:17:52 PM GMT+9
# Code: 

import numpy as np
from POA.utilities import normalize
from POA.Curves import Bernstein

class HermiteInterpolation():
    # PH Quintic을 사용해서 평면상의 Hermite Interpolation을 수행하는 클래스
    def __init__(self, Ps, Ts, scaling, tau=None, mode=0):
        Pi, Pf = Ps
        Ti, Tf = Ts
        if tau is None:
            ndata = 100
            tau = np.linspace(0,1,ndata)
        Ti = np.reshape(scaling*normalize(np.reshape(Ti,newshape=[2,])),newshape=[2])
        Tf = np.reshape(scaling*normalize(np.reshape(Tf,newshape=[2,])),newshape=[2])
        Delta_p0 = Ti
        Delta_p4 = Tf
        p1 = Pi + Ti
        p4 = Pf - Tf
        Delta_p4p0 = p4 - p1
        if mode == 0:
            sign1 = 1
            sign2 = 1
        elif mode == 1:
            sign1 = -1
            sign2 = 1
        elif mode == 2:
            sign1 = 1
            sign2 = -1
        else:
            sign1 = -1
            sign2 = -1
        self.sign1 = sign1
        self.sign2 = sign2
            
        self.us, self.vs    = HermiteInterpolation.get_uvs(Delta_p0, Delta_p4, Delta_p4p0, sign1, sign2)
        self.cps            = HermiteInterpolation.get_control_points(self.us, self.vs, Pi)
        self.sigma_coef     = self.get_sigma_coefficients()
        self.parametric_speed = Bernstein.Bernstein_polynomial(self.sigma_coef)
        self.Pi = Pi
        self.Pf = Pf
        self.Ti = Ti
        self.Tf = Tf
        self.tau = tau

        # Bezier curve
        betas_x = [self.cps[i,0] for i in range(0,6)]
        betas_y = [self.cps[i,1] for i in range(0,6)]
        self.ph_x = Bernstein.Bernstein_polynomial(betas_x)
        self.ph_y = Bernstein.Bernstein_polynomial(betas_y)

    def visualize(self, fig=None):
        if fig is None:
            fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=self.ph_x(self.tau),
            y=self.ph_y(self.tau),
            name='PH curve'
        ))
        betas_x = [self.cps[i,0] for i in range(0,6)]
        betas_y = [self.cps[i,1] for i in range(0,6)]
        fig.add_trace(go.Scatter(
            x = betas_x,
            y = betas_y,
            line = dict(dash='dash'),
            name='Control points'
        ))
        return fig
        
    @staticmethod
    def nonzero_sign(x):
        return 2*(x>=0).astype(int)-1

    # 식 (25.4)
    @staticmethod
    def get_uvs(Delta_p0, Delta_p4, Delta_p4p0, sign1, sign2):
        dx0, dy0 = Delta_p0
        dx4, dy4 = Delta_p4
        u0 = sign1*np.sqrt(5/2)*np.sqrt(np.linalg.norm(Delta_p0)+dx0)
        v0 = sign1*np.sqrt(5/2)*HermiteInterpolation.nonzero_sign(dy0)*np.sqrt(np.linalg.norm(Delta_p0)-dx0)
        u2 = sign2*np.sqrt(5/2)*np.sqrt(np.linalg.norm(Delta_p4)+dx4)
        v2 = sign2*np.sqrt(5/2)*HermiteInterpolation.nonzero_sign(dy4)*np.sqrt(np.linalg.norm(Delta_p4)-dx4)
        a = 9/16*(u0**2-v0**2+u2**2-v2**2)+5/8*(u0*u2-v0*v2)+15/2*(Delta_p4p0[0])
        b = 9/8*(u0*v0+u2*v2)+5/8*(u0*v2+u2*v0)+15/2*(Delta_p4p0[1])
        c = np.sqrt(a**2+b**2)
        u1 = -3/4*(u0+u2) + 1/np.sqrt(2)*np.sqrt(c+a)
        v1 = -3/4*(v0+v2) + np.sign(b)*1/np.sqrt(2)*np.sqrt(c-a)
        return [u0, u1, u2], [v0, v1, v2]

    # 식 (17.6)
    @staticmethod
    def get_control_points(us, vs, p0):
        p0 = np.reshape(p0, newshape=[2])
        u0, u1, u2 = us
        v0, v1, v2 = vs
        p1 = p0 + 1/5*np.array([u0**2-v0**2, 2*u0*v0])
        p2 = p1 + 1/5*np.array([u0*u1-v0*v1, u0*v1+u1*v0])
        p3 = p2 + 2/15*np.array([u1**2-v1**2, 2*u1*v1]) + 1/15*np.array([u0*u2-v0*v2,u0*v2+u2*v0])
        p4 = p3 + 1/5*np.array([u1*u2-v1*v2, u1*v2+u2*v1])
        p5 = p4 + 1/5*np.array([u2**2-v2**2, 2*u2*v2])
        return np.array([p0, p1, p2, p3, p4, p5])

    # 식 (17.12)
    def get_sigma_coefficients(self):
        u0, u1, u2 = self.us
        v0, v1, v2 = self.vs
        sigma       = [None]*5
        sigma[0]    = u0**2 + v0**2
        sigma[1]    = u0*u1 + v0*v1
        sigma[2]    = 2/3*(u1**2 + v1**2) + 1/3*(u0*u2 + v0*v2)
        sigma[3]    = u1*u2 + v1*v2
        sigma[4]    = u2**2 + v2**2
        return sigma

    def eval_sigma(self, tau=None):
        if tau is None:
            tau = self.tau
        return np.array(self.parametric_speed(tau))

    def eval(self, tau=None):
        if tau is None:
            tau = self.tau
        return np.array([self.ph_x(tau), self.ph_y(tau)])
    
    def tangent(self, tau):
        u_fun = Bernstein.Bernstein_polynomial(self.us)
        v_fun = Bernstein.Bernstein_polynomial(self.vs)
        x_prime_fun = lambda tau: (u_fun(tau)**2 - v_fun(tau)**2)
        y_prime_fun = lambda tau: 2*u_fun(tau)*v_fun(tau)
        return x_prime_fun(tau), y_prime_fun(tau)
    
    def eval_u(self, tau):
        u = Bernstein.Bernstein_polynomial(self.us)
        return u(tau)

    def first_derivative(self, tau):
        return self.tangent(tau)

    def second_derivative(self, tau):
        """
        곡선의 매개변수 tau위치에서의 2차 미분 벡터를 계산한다.
        """
        u_fun = Bernstein.Bernstein_polynomial(self.us)
        v_fun = Bernstein.Bernstein_polynomial(self.vs)
        u_prime_fun = Bernstein.Bernstein_derivative_polynomial(self.us,p=1)
        v_prime_fun = Bernstein.Bernstein_derivative_polynomial(self.vs,p=1)
        x_2prime_fun = lambda tau: 2*(u_fun(tau)*u_prime_fun(tau) - v_fun(tau)*v_prime_fun(tau))
        y_2prime_fun = lambda tau: 2*(u_prime_fun(tau)*v_fun(tau) - u_fun(tau)*v_prime_fun(tau))
        return x_2prime_fun(tau), y_2prime_fun(tau)

    def third_derivative(self, tau):
        """
        곡선의 매개변수 tau위치에서의 3차 미분 벡터를 계산한다.
        """
        u_fun = Bernstein.Bernstein_polynomial(self.us)
        v_fun = Bernstein.Bernstein_polynomial(self.vs)
        u_prime_fun = Bernstein.Bernstein_derivative_polynomial(self.us,p=1)
        v_prime_fun = Bernstein.Bernstein_derivative_polynomial(self.vs,p=1)
        u_2prime_fun = Bernstein.Bernstein_derivative_polynomial(self.us,p=2)
        v_2prime_fun = Bernstein.Bernstein_derivative_polynomial(self.vs,p=2)
        x_3prime_fun = lambda tau: 2*(u_prime_fun(tau)**2 + u_fun(tau)*u_2prime_fun(tau)-v_prime_fun(tau)**2-v_fun(tau)*v_2prime_fun(tau))
        y_3prime_fun = lambda tau: 2*(u_2prime_fun(tau)*v_fun(tau) + 2*u_prime_fun(tau)*v_prime_fun(tau) + u_fun(tau)*v_2prime_fun(tau))
        return x_3prime_fun(tau), y_3prime_fun(tau)
