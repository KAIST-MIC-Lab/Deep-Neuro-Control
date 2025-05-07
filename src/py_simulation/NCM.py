import torch
import numpy as np
from train_src.SimpleNN import SimpleNN

class NCM:
    def __init__(self, nn_path):
        self.R = np.diag([1, 1])
        self.inv_R = np.linalg.inv(self.R)

        self.loadNN(nn_path)

    def loadNN(self, nn_path):
        self.nn = SimpleNN()
        self.nn.load_state_dict(torch.load(nn_path))
        self.nn.eval()

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

    def getControl(self, x, xd, ud):
        # x = x.reshape(-1, 2, 1)
        # xd = xd.reshape(-1, 2, 1)
        # ud = ud.reshape(-1, 2, 1)

        # nn_inputs = np.concatenate((x, xd, ud), axis=1)
        nn_inputs = np.concatenate((x, xd), axis=1)

        nn_inputs = nn_inputs.reshape(-1, 4)

        M_chol = self.nn(torch.tensor(nn_inputs).float()).detach().numpy()
        print(M_chol.shape)
        M_sqrt = np.array([
            [M_chol[0][0], 0],
            [M_chol[0][1], M_chol[0][2]]
        ])
        M = M_sqrt @ M_sqrt.T

        A, B = self.getSytemGradient(x, xd)

        u = ud - self.inv_R@B.T@M@(x - xd)

        return u, M