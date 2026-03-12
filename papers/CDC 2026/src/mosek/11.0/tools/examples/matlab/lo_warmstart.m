%%
%  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File:      lo_warmstart.m
%
%  Purpose: Demonstrates warm-starting the simplex optimizer
%           using the problem
%
%  maximize 3*x0 + 1*x1 + 5*x2 + x3
%  such that
%           3*x0 + 1*x1 + 2*x2        = 30,
%           2*x0 + 1*x1 + 3*x2 + 1*x3 > 15,
%                  2*x1 +      + 3*x3 < 25
%  and
%           x0,x1,x2,x3 > 0,
%           0 < x1 < 10
%%
function [] = lo_warmstart()
    % We set up the linear model
    % together with:
    % - approximate primal solution
    % - status keys
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
                  conname = ["c1" "c2" "c3"]', ...
                  solution_x = [0 0 15 8.3 0 0]', ...
                  solution_skx = ["LOW"  "LOW"  "BAS"  "BAS"  "BAS"  "UPR"]', ...
                  solution_skc = ["FIX"  "FIX", "FIX"]);

    % Choose the simplex optimizer to solve with
    model.solve(param = [ "MSK_IPAR_OPTIMIZER", "MSK_OPTIMIZER_FREE_SIMPLEX" ]);

    % Check the solution status and fetch the solution
    [hassol, prosta, solsta] = model.hassolution();
    if hassol && solsta == "OPTIMAL"
        fprintf("Solution, problem status = %s, solution status = %s:\n", prosta, solsta);

        % Show solution
        disp("Solution, status keys and duals");
        disp([(1:6)' model.getsolution("x") model.getsolution("skx") model.getsolution("slx") model.getsolution("sux")]);
        disp("Constraint status and duals");
        disp([(1:3)' model.getsolution("skc") model.getsolution("y")]);

        % Solver statistics
        fprintf("Optimizer time %.3f\n", model.info.MSK_DINF_OPTIMIZER_TIME);
        fprintf("#iterations    %ld\n",  model.info.MSK_LIINF_SIMPLEX_ITER);
    end
end

