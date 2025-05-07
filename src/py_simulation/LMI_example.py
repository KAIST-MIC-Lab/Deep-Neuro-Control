# example from https://www.cvxpy.org
# User Guide: https://www.cvxpy.org/tutorial/index.html


import cvxpy as cp
import numpy as np

# Problem data.
bar_d = 100
n = 2
mu = 5
alp = 0.02 # shall be determined
dt = 1e-2
np.random.seed(1)

pre_W = np.array([1, 0, 0, 1]).reshape((2,2))
dfdx = np.array([1, 0, 0, 1]).reshape((2,2))
# dfdx = np.array([10, 4, 1, 10]).reshape((2,2))

#1 Construct the problem.
X = cp.Variable()
W = cp.Variable((n,n))

objective = cp.Minimize(bar_d/alp*X)
# objective = cp.Minimize(cp.sum_squares(bar_d/alp*X))
# objective = cp.Minimize(X)
constraints = [
    W == W.T,
    mu*(W-pre_W)/dt - mu*(W.T@dfdx.T + dfdx@W) >= mu*2*alp*W,
    np.eye(n) <= mu*W,
    mu*W <= X * np.eye(n),
    ]
prob = cp.Problem(objective, constraints)

# The optimal objective value is returned by `prob.solve()`.
result = prob.solve()

print("status:", prob.status)
print("optimal value", result)
print("optimal X", X.value)
print("optimal W", W.value)

print(constraints[0].dual_value)
