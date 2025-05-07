import numpy as np
import cvxpy as cp

class CVSTEM:
    def __init__(self, dt):
        self.env = None

        self.dt = dt;

        self.x = np.zeros([2, 1])

        self.m_MAX = 1.0
        self.m_MIN = 0.001
        self.b_MAX = 1    

        self.mu = self.m_MAX    

        self.alpha = .1

        self.u_num = 2

        self.pre_W = np.eye(self.u_num)
        self.R = np.diag([1, 1])
        self.inv_R = np.linalg.inv(self.R)

        return
    
    def getSytemGradient(self, x, xd):
        A = np.array([
            [.3, 1],
            [-2, -.3]
        ])
        B = np.array([
            [1, 0.1],
            [0.1, 1]
        ])

        return A,B
    
    def LMIsolver(self, A, B):
        pre_W = self.pre_W
        n = self.u_num
        mu = self.mu
        dt = self.dt
        inv_R = self.inv_R  

        X = cp.Variable()
        W = cp.Variable((self.u_num, self.u_num))

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
            W_opt = pre_W
            X_opt = None

        self.pre_W = W_opt

        return W_opt
    
    def getControl(self, x, xd, ud):
        (A, B) = self.getSytemGradient(x, xd)

        W = self.LMIsolver(A, B)
        M = np.linalg.inv(W)

        u = ud - self.inv_R@B.T@M@(x - xd)

        return u, M