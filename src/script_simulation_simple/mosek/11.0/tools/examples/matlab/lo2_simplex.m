%%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      lo2_simplex.m
%
%  Purpose :   Demonstrates how to solve small linear
%              optimization problem using the MOSEK API for Matlab
%              in multiple calls.
%%
function [xx,prosta,solsta] = lo2_simplex()
    % Create an empty model
    model = moseklinmodel(name = "lo2");

    % Add some variables with bounds
    model.appendvars(6,...
                     bl = [ 0.0  0.0 0.0 0.0   0.0 -inf ]', ...
                     bu = [ inf 10.0 inf inf   inf  0.0 ]');
    model.varname([1:6],["x1" "x2" "x3" "x4" "y1" "y2"]);
    
    % Set the objective
    model.objective("maximize", [ 3 1 5 1 0 0]', objfixterm = 0.0);

    % Constraint 3x1 + x2 + 2x3 == 30
    model.appendcons([ 3 1 2 0 0 0], [30]);

    % The two remaining constraints
    model.appendcons([ 2 1 3 1 -1  0 ; ...
                       0 2 0 3  0 -1 ], ...
                     [ 15 25 ]');

    % Constraint names
    model.conname([1:3],["c1" "c2" "c3"]);
    
    % Choose the simplex optimizer
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
