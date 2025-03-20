# [1] Y. Lee, S.-H. Lee, and C. C. Chung, “Lpv H∞ control with disturbance estimation for permanent magnet synchronous motors,” IEEE Transactions on Industrial Electronics, vol. 65, no. 1, pp. 488–497, 2018.

import numpy as np

class Env:
    def __init__(self, dt):
        self.env = None

        self.dt = dt;

        self.x = np.zeros([2, 1])

        self.L = .66        # Inductance (mH)
        self.R = 2.51       # Resistance (Ohm)
        self.J = 3.24e-5    # Inertia (kg.m^2)
        self.Phi = 16.8e-3  # Flux (Wb)
        self.P = 4          # Pole pairs
        self.fv = 2e-3      # Viscous friction (N.m.s/rad)

        self.x_num = 4
        self.u_num = 2

        return
    
    def __str__(self):
        print("ENV INFORMATION")
        print("x: ", self.x)

        return
    
    def reset(self, x0):
        self.x = x0
        
        return
    
    def step(self, u):
        [th, w, ia, ib] = self.x.reshape(-1)
        [va, vb] = u.reshape(-1)

        L = self.L
        R = self.R
        J = self.J
        Phi = self.Phi
        P = self.P
        fv = self.fv

        trq = -(3/2)*P*Phi*np.sin(P*th)*ia + (3/2)*P*Phi*np.cos(P*th)*ib
        trqL = 0;

        grad = np.array([
            w,
            (1/J)*(trq - fv*w - trqL),
            (1/L)*(va - R*ia - P*Phi*w*np.sin(P*th)),
            (1/L)*(vb - R*ib + P*Phi*w*np.cos(P*th)),
        ]).reshape(4,1)

        self.x = self.x + grad * self.dt

        return
    
    def render(self):
        pass

    def close(self):
        pass