%   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%   File:      gp1.m
%
%   Purpose:   Demonstrates how to solve a simple Geometric Program (GP)
%              cast into conic form with exponential cones and log-sum-exp.
%
%              Example from
%                https://gpkit.readthedocs.io/en/latest/examples.html%maximizing-the-volume-of-a-box
function [xx,prosta,solsta] = gp1()
    Awall  = 200
    Afloor = 50
    alpha  = 2
    beta   = 10
    gamma  = 2
    delta  = 10
 
    % A model with variables [x,y,z,u1,u2]
    model = mosekmodel(...
        name = "gp1",...
        numvar = 5);
    model.objective("maximize", [ 1 1 1 0 0]);

    % u1 + u2 = 1
    model.appendcons(F = [ 0 0 0 1 1 ], domain = mosekdomain("equal",     rhs=1.0));
    
    % y + z <= log(Afloor)
    model.appendcons(F = [ 0 1 1 0 0 ], domain = mosekdomain("less than", rhs=log(Afloor)));
    
    % Two-sided bounds on x-y and z-y
    model.appendcons(F = [ 1 -1 0 0 0 ; ...
                           0 -1 1 0 0 ], domain = mosekdomain("greater than",dim=2, rhs=[log(alpha) log(gamma)]'));
    model.appendcons(F = [ 1 -1 0 0 0 ; ...
                           0 -1 1 0 0 ], domain = mosekdomain("less than",   dim=2, rhs=[log(beta)  log(delta)]'));

    % Conic constraints
    model.appendcons(F = sparse([1 3 3],[4 1 2],[1.0 1.0 1.0]), g = [0 1 log(alpha/Awall)]', domain = mosekdomain("exp"));
    model.appendcons(F = sparse([1 3 3],[5 1 3],[1.0 1.0 1.0]), g = [0 1 log(alpha/Awall)]', domain = mosekdomain("exp"));
    
    model.solve();
    if model.hassolution("interior")
        [xx,prosta,solsta] = model.getsolution("interior","x");
        x = xx(1:3);
        disp(exp(x));
    end
end
