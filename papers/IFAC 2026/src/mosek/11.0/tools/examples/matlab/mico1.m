%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File      : mico1.m
%
%  Purpose   : Demonstrates how to solve a small mixed
%              integer conic optimization problem.
%
%              minimize    x^2 + y^2
%              subject to  x >= e^y + 3.8
%                          x, y - integer
%

function mico1()

    %The full variable is [t; x; y]
    model = mosekmodel(name = "mi-conic", ...
                       numvar = 3, objective = [1 0 0], objsense = "min", ...
                       F = [ speye(3); ...
                             0 1 0;    ...
                             0 0 0;    ...
                             0 0 1 ],  ...
                       g = [0 0 0 -3.8 1 0], ...
                       domain = [mosekdomain("qcone", dim = 3), ...
                                 mosekdomain("exp")], ...
                       intvars = [2, 3]);     % Specify the indices of integer variables

    % It is as always possible (but not required) to input an initial solution
    % to start the mixed-integer solver. 
    model.setsolution("x", [100, 9, -1]);

    model.solve();
    x = model.getsolution("integer", "x");
    fprintf("Optimal (x,y) = (%.2e, %.2e)\n", x(2), x(3));
end
