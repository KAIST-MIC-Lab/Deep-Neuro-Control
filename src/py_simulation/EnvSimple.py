import numpy as np

class Env:
    def __init__(self, dt):
        self.env = None

        self.dt = dt;

        self.x = np.zeros([2, 1])

        self.A = np.array([
            [0, 1],
            [-2, -3]
        ])
        self.B = np.array([
            [0],
            [1]
        ])

        self.x_num = 2
        self.u_num = 1

        return
    
    def __str__(self):
        print("ENV INFORMATION")
        print("x: ", self.x)

        return
    
    def reset(self, x0):
        self.x = x0
        
        return
    
    def step(self, u):
        grad = np.dot(self.A, self.x) + np.dot(self.B, u)
        self.x = self.x + grad * self.dt

        return
    
    def render(self):
        pass

    def close(self):
        pass