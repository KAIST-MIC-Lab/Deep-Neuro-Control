##
#   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File:      djc1.py
#
#   Purpose: Demonstrates how to solve the problem with two disjunctions:
#
#      minimize    2x0 + x1 + 3x2 + x3
#      subject to   x0 + x1 + x2 + x3 >= -10
#                  (x0-2x1<=-1 and x2=x3=0) or (x2-3x3<=-2 and x1=x2=0)
#                  x0=2.5 or x1=2.5 or x2=2.5 or x3=2.5
##

from mosek.fusion import *
import mosek.fusion.pythonic
import sys

with Model('djc1') as M:
    # Create a variable of length 4
    x = M.variable('x', 4)

    # First disjunctive constraint
    M.disjunction( (x[0]-2*x[1]<=-1) & (x[2:4]==0) 
                  | 
                   (x[2]-3*x[3]<=-2) & (x[0:2]==0))
    
    # Second disjunctive constraint
    # Array of terms reading x_i = 2.5 for i = 0,1,2,3
    M.disjunction( [ x[i] == 2.5 for i in range(4) ])

    # The linear constraint
    M.constraint(Expr.sum(x) >= -10)

    # Objective
    M.objective(ObjectiveSense.Minimize, Expr.dot([2,1,3,1], x))

    # Useful for debugging
    M.writeTask("djc.ptf")          # Save to a readable file
    M.setLogHandler(sys.stdout)     # Enable log output

    # Solve 
    M.solve()

    if M.getPrimalSolutionStatus() == SolutionStatus.Optimal:
        print(f"Solution ={x.level()}")
    else:
        print("Another solution status")
