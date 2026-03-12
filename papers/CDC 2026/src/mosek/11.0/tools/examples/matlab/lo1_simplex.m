%%
%  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File:      lo1_simplex.m
%
%  Purpose: Demonstrates how to solve the problem
%
%  maximize 3*x0 + 1*x1 + 5*x2 + x3
%  such that
%           3*x0 + 1*x1 + 2*x2        = 30,
%           2*x0 + 1*x1 + 3*x2 + 1*x3 > 15,
%                  2*x1 +      + 3*x3 < 25
%  and
%           x0,x1,x2,x3 > 0,
%           0 < x1 < 10
%
%  using the linear/simplex toolbox.
%%
function [xx,prosta,solsta] = lo1_simplex()
    model = moseklinmodel(...
                  name = "lo1", ...
                  numvar = 6, ...
                  objsense = "maximize", ...
                  objective = [ 3 1 5 1 0 0 ]', ...
                  A = [ 3 1 2 0  0  0 ; ...
                        2 1 3 1 -1  0 ; ...
                        0 2 0 3  0 -1 ], ...
                  b = [ 30 15 25 ]', ...
                  blx = [ 0.0  0.0 0.0 0.0   0.0 -inf ]', ...
                  bux = [ inf 10.0 inf inf   inf  0.0 ]', ...
                  varnames = ["x1" "x2" "x3" "x4" "y1" "y2"]',...
                  conname = ["c1" "c2" "c3"]');

    % Choose the simplex optimizer to solve with
    model.solve(param = [ "MSK_IPAR_OPTIMIZER", "MSK_OPTIMIZER_FREE_SIMPLEX" ]);

    % Check the solution status and fetch the solution
    [hassol, prosta, solsta] = model.hassolution();
    if hassol && solsta == "OPTIMAL"
        fprintf("Solution, problem status = %s, solution status = %s:\n", prosta, solsta);

        % Primal variables
        result = model.getsolution("x");
        x = result(1:4);
        s = result(5:6);
        fprintf("  x[%d] = %.4f\n",[1:4 ; x']);
        fprintf("Slacks:\n");
        fprintf("  s[%d] = %.4f\n",[1:2 ; s']);

        % Other elements of the solution
        fprintf("Duals:\n");
        fprintf("  y[%d] = %.4f\n",[1:3 ; model.getsolution("y")']);
        fprintf("slx[%d] = %.4f\n",[1:6 ; model.getsolution("slx")']);
        fprintf("sux[%d] = %.4f\n",[1:6 ; model.getsolution("sux")']);

        fprintf("Status keys of constraints:\n");
        fprintf("skc[%s] = %s\n",[1:3 ; model.getsolution("skc")']);

        fprintf("Status keys of variables:\n");
        fprintf("skx[%s] = %s\n",[1:6 ; model.getsolution("skx")']);
    end
end

