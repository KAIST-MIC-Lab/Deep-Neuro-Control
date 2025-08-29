%%
%   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%   File:      pow1.m
%
%   Purpose: Demonstrates how to solve the problem
%
%     maximize x^0.2*y^0.8 + z^0.4 - x
%           st x + y + 0.5z = 2
%              x,y,z >= 0
%
%%

function [xx,prosta,solsta] = pow1()
    model = mosekmodel(...
        name = "pow1",...
        numvar = 5);
    model.objective("maximize", [-1 0 0 1 1]);
    model.appendcons(F = [1 1 0.5], domain = mosekdomain("equal",rhs=2));
    model.appendcons(F = sparse([1 2 3],[1 2 4],[1.0 1.0 1.0]),       domain = mosekdomain("pow", dim=3, alpha=[0.2 0.8]));
    model.appendcons(F = sparse([1 3],[3 5],[1.0 1.0]), g = [0 1 0]', domain = mosekdomain("pow", dim=3, alpha=[0.4 0.6]));
    
    model.solve();
    if model.hassolution("interior")
        [xx,prosta,solsta] = model.getsolution("interior","x");
    end
end
