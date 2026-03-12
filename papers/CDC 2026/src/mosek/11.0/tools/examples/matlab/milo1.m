%%
%    Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%    File:    milo1.m
%
%    Purpose:  Demonstrates how to solve a small mixed
%              integer linear optimization problem.
%%

function [xx,prosta,solsta] = milo1()
    model = mosekmodel(name="milo1", numvar=2, ...
                       intvars=[1 2]');    % Specify indices of integer variables

    model.objective("maximize", [ 1.0 0.64 ]);
    model.appendcons(F = speye(2),    domain = mosekdomain("rplus",dim=2));
    model.appendcons(F = [50.0 31.0], domain = mosekdomain("less than", rhs=250.0));
    model.appendcons(F = [4.0  -2.0], domain = mosekdomain("greater than", rhs=-4));

    model.solve();

    % Access the integer solution
    if model.hassolution("integer")
        [xx,prosta,solsta] = model.getsolution("integer", "x");
        fprintf("Solution: %s\n", sprintf("%g ", xx));
    else
        disp("Solve failed");
    end
end
