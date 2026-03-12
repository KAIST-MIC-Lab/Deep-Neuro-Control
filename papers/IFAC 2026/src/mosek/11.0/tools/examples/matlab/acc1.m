%%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      acc1.m
%
%  Purpose :   Tutorial example for affine conic constraints.
%              Models the problem:
%
%              maximize c^T x
%              subject to  sum(x) = 1
%                          gamma >= |Gx+h|_2
%%

function [xx,prosta,solsta] = acc1()
    % Example data
    c = [2 3 -1]';
    gamma = 0.03;
    G = [1.5 0.1 0; ...
         0.3 0   2.1];
    h = [0 0.1]';

    n = size(c, 1);
    k = size(h, 1);

    % Initialize the model
    model = mosekmodel(name = "acc1",...
                       numvar = n);
    
    % Set objective vector
    model.objective("maximize", c);

    % The constraint sum(x) = 1
    model.appendcons(F = ones(1,n), domain = mosekdomain("equal", rhs=1));

    % The conic quadratic constraint Fx+g \in \Quad with k+1 rows
    model.appendcons(F = sparse([zeros(1,n); G]), ...
                     g = [gamma; h], ...
                     domain = mosekdomain("quadratic cone", dim=k+1));

    % Solve the problem
    model.solve();

    % Check if solution is available
    [hassol, prosta, solsta] = model.hassolution("interior"); 

    if hassol && solsta == "OPTIMAL"
        % Get primal solution
        xx = model.getsolution("interior", "x");

        % Get dual solution
        y = model.getsolution("interior", "y");

        disp("Primal solution");
        disp(xx);
    end

end
