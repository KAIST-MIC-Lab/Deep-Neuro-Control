%%
%   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%   File:      cqo1.m
%
%   Purpose: Demonstrates how to solve the problem
%
%   minimize y1 + y2 + y3
%   such that
%            x1 + x2 + 2.0 x3  = 1.0
%                    x1,x2,x3 >= 0.0
%   and
%            (y1,x1,x2) in C_3,
%            (y2,y3,x3) in K_3
%
%   where C_3 and K_3 are respectively the quadratic and
%   rotated quadratic cone of size 3 defined as
%       C_3 = { z1,z2,z3 :      z1 >= sqrt(z2^2 + z3^2) }
%       K_3 = { z1,z2,z3 : 2 z1 z2 >= z3^2              }
%%

function [xx,prosta,solsta] = cqo1()
    model = mosekmodel(...
        name = "cqo1", ...
        numvar = 6);

    % Variable is [x, y]
    model.objective("minimize", [0 0 0 1 1 1]');

    % Linear constraint
    model.appendcons(F = [1 1 2 0 0 0], domain = mosekdomain("eq", rhs=1.0));
    
    % Bounds on x
    model.appendcons(F = sparse([eye(3) zeros(3)]), domain = mosekdomain("rplus", dim=3));

    % Quadratic cone
    model.appendcons(F = sparse([1,2,3],[4,1,2],[1,1,1]), domain = mosekdomain("qcone",dim=3));
    
    % Rotated quadratic cone
    model.appendcons(F = sparse([1,2,3],[5,6,3],[1,1,1]), domain = mosekdomain("rqcone",dim=3));

    model.solve();
    if model.hassolution("interior")
        [xx,prosta,solsta] = model.getsolution("interior","x");
        fprintf("Solution x1,x2,x3: [ %s ]\n", sprintf("%.2g ", xx(1:3)));
        fprintf("Solution y1,y2,y3: [ %s ]\n", sprintf("%.2g ", xx(4:6)));
    end

end
