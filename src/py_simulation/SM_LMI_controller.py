from EnvPMSM import Env
import numpy as np
import cvxpy as cp

class CVSTEM(Env):
    def __init__(self, dt):
        super().__init__(dt)

        self.env = None

        self.dt = dt;

        self.x = np.zeros([4, 1])

        self.m_MAX = 1e-2
        # self.m_MIN = self.m_MAX *1e-2
        self.m_MIN = 1e-4
        self.b_MAX = 1e-3    

        self.mu = self.m_MAX    

        self.alpha = 1e3

        self.x_num = 4
        self.u_num = 2

        self.pre_W = np.eye(self.x_num) * 1e0
        # self.pre_W = np.zeros((self.x_num, self.x_num))
        R = np.diag([1, 1]) * 1e3
        self.inv_R = np.linalg.inv(R)

        return
    
    def getSytemGradient(self, x, xd, ud):
        R = self.R
        L = self.L
        J = self.J

        fv = self.fv
        P = self.P
        Phi = self.Phi
        
        th =    x[0,0]
        w =     x[1,0]
        ia =    x[2,0]
        ib =     x[3,0]

        thd =   xd[0,0]
        wd =    xd[1,0]
        iad =   xd[2,0]
        ibd =   xd[3,0]

        f1 = -(3/2)*P*Phi*np.sin(P*th)
        f2 = +(3/2)*P*Phi*np.cos(P*th)
        f3 = +P*Phi*np.sin(P*th)
        f4 = -P*Phi*np.cos(P*th)

        fd1 = -(3/2)*P*Phi*np.sin(P*thd)
        fd2 = +(3/2)*P*Phi*np.cos(P*thd)
        fd3 = +P*Phi*np.sin(P*thd)
        fd4 = -P*Phi*np.cos(P*thd)

        A = np.array([
            [0, 1, 0, 0],
            [0, -fv/J, f1/J, f2/J],
            [0, f3/L, -R/L, 0],
            [0, f4/L, 0, -R/L]
        ])

        Ad = np.array([
            [0, 1, 0, 0],
            [0, -fv/J, fd1/J, fd2/J],
            [0, fd3/L, -R/L, 0],
            [0, fd4/L, 0, -R/L]
        ])

        # SDC = A - Ad

        # SDC = np.eye(4) * -1e0
        # SDC = np.array([
        #     [0, 1, 0, 0],
        #     [0, -fv/J, 0, 0],
        #     [0, 0, -R/L, 0],
        #     [0, 0, 0, -R/L]
        # ])

        if abs(ia-iad) < 1e-6:
            delta_f1 = 0
        else:
            delta_f1 = (f1-fd1)/(ia-iad)/J
        if abs(ib-ibd) < 1e-6:
            delta_f2 = 0
        else:
            delta_f2 = (f2-fd2)/(ib-ibd)/J
        if abs(w-wd) < 1e-6:
            delta_f3 = 0
            delta_f4 = 0
        else:
            delta_f3 = (f3-fd3)/(w-wd)/L
            delta_f4 = (f4-fd4)/(w-wd)/L
    
        SDC = np.array([
            [0, 1, 0, 0],
            [0, -fv/J, delta_f1, delta_f2],
            [0, delta_f3, -R/L, 0],
            [0, delta_f4, 0, -R/L]
        ])


        # print(SDC)
        # A = self.getSDC(x, xd, ud)

        B = np.array([
            [0, 0],
            [0, 0],
            [1/L, 0],
            [0, 1/L]
        ])

        return SDC,B
    
    def getSDC(self, x, xd, ud):
        f, B = self.getSystemFunction(x)
        fd, Bd = self.getSystemFunction(xd)

        leftVec = f+B@ud - fd-Bd@ud
        rightVec = x-xd

        try:
            pInv = np.linalg.inv(rightVec@rightVec.T)
        except np.linalg.LinAlgError:
            pInv = np.eye(rightVec.shape[0]) * 0.001
            print(pInv)

        A = leftVec@rightVec.T@pInv

        return A

    def getSystemFunction(self, x):
        [th, w, ia, ib] = x.reshape(-1)

        L = self.L
        R = self.R
        J = self.J
        Phi = self.Phi
        P = self.P
        fv = self.fv

        trq = -(3/2)*P*Phi*np.sin(P*th)*ia + (3/2)*P*Phi*np.cos(P*th)*ib
        trqL = 0;

        f = np.array([
            w,
            (1/J)*(trq - fv*w - trqL),
            (1/L)*(- R*ia + P*Phi*w*np.sin(P*th)),
            (1/L)*(- R*ib - P*Phi*w*np.cos(P*th)),
        ]).reshape(4,1)

        B = np.array([
            [0, 0],
            [0, 0],
            [1/L, 0],
            [0, 1/L]
        ])

        return f, B
    
    def LMIsolver(self, A, B):
        pre_W = self.pre_W
        n = self.x_num
        mu = self.mu
        dt = self.dt
        inv_R = self.inv_R  

        X = cp.Variable()
        W = cp.Variable((self.x_num, self.x_num))

        objective = cp.Minimize(X)
        constraints = [
            W == W.T,
            -mu*(W-pre_W)/dt + mu*(W.T@A.T + A@W) - mu*2*B@inv_R@B.T <= -mu*2*self.alpha*W,
            np.eye(n) <= mu*W,
            mu*W <= X * np.eye(n),
            ]
        prob = cp.Problem(objective, constraints)

        result = prob.solve()

        if prob.status == cp.OPTIMAL:
            W_opt = W.value
            X_opt = X.value
        else:
            print("Problem is not optimal")
            print("status:", prob.status)
            print("optimal value", result)
            print(pre_W)
            W_opt = pre_W
            X_opt = None

        self.pre_W = W_opt

        return W_opt
    
    def getControl(self, x, xd, ud):
        (A, B) = self.getSytemGradient(x, xd, ud)

        W = self.LMIsolver(A, B)
        M = np.linalg.inv(W)

        u = ud - self.inv_R@B.T@M@(x - xd)

        return u, M