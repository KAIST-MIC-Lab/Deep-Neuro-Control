%
%   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%   File:      ceo1.m
%
%   Purpose:   To demonstrate how to solve a small conic exponential
%              optimization problem using the MOSEK API.
%
function [xx,prosta,solsta] = ceo1()
    model = mosekmodel(...
        name = "ceo1",...
        numvar = 3);
 
    % Minimize x1+x2
    model.objective("minimize", [1 1 0]);
 
    % Subject to: sum(x) == 1
    model.appendcons(F = [1 1 1], domain = mosekdomain("equal", rhs = 1.0));
 
    % And subject to: (x1, x2, x3) \in ExpCone
    model.appendcons(F = speye(3), ...
                     domain = mosekdomain("exp", dim = 3));
    
    % Solve the model and retrieve the solution
    model.solve();
    if model.hassolution("interior")
        [xx,prosta,solsta] = model.getsolution("interior","x");
        disp(xx);
    end
end
