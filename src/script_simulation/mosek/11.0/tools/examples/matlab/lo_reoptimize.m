%%
%  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File:      lo_reoptimize.m
%
%  Purpose: Demonstrates warm-starting the simplex optimizer
%           during reoptimizations of problems starting from
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
function [] = lo_reoptimize()
    % We set up and initially solve a sample linear model
    model = moseklinmodel(...
                  numvar = 6, ...
                  objsense = "maximize", ...
                  objective = [ 3 1 5 1 0 0 ]', ...
                  A = [ 3 1 2 0  0  0 ; ...
                        2 1 3 1 -1  0 ; ...
                        0 2 0 3  0 -1 ], ...
                  b = [ 30 15 25 ]', ...
                  blx = [ 0.0  0.0 0.0 0.0   0.0 -inf ]', ...
                  bux = [ inf 10.0 inf inf   inf  0.0 ]');
    
    model.solve(param = [ "MSK_IPAR_OPTIMIZER", "MSK_OPTIMIZER_FREE_SIMPLEX" ], quiet = true);    

    checksolution(model);

    % Introduce a stricter upper bound for the 3-rd constraint
    model.setb([30, 15, 22]');

    % Solve, using dual simplex
    model.solve(param = [ "MSK_IPAR_OPTIMIZER", "MSK_OPTIMIZER_DUAL_SIMPLEX" ]);

    checksolution(model);

    % Change an objective coefficient
    model.setc([3 2 0 2]');

    % Solve, using primal simplex
    model.solve(param = [ "MSK_IPAR_OPTIMIZER", "MSK_OPTIMIZER_PRIMAL_SIMPLEX" ]);

    checksolution(model);
end

function [] = checksolution(model)
    [hassol, prosta, solsta] = model.hassolution();
    if hassol && solsta == "OPTIMAL"
        fprintf("Solution, problem status = %s, solution status = %s:\n", prosta, solsta);
        disp("Solution, status keys and duals");
        disp([(1:6)' model.getsolution("x") model.getsolution("skx") model.getsolution("slx") model.getsolution("sux")]);
        disp("Constraint status and duals");
        disp([(1:3)' model.getsolution("skc") model.getsolution("y")]);
        fprintf("Optimizer time %.3f\n", model.info.MSK_DINF_OPTIMIZER_TIME);
        fprintf("#iterations    %ld\n",  model.info.MSK_LIINF_SIMPLEX_ITER);
    else
        fprintf("Solution, unexpected problem status = %s, solution status = %s:\n", prosta, solsta);
    end
end
