%%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      mioinitsol.m
%
%   Purpose:   To demonstrate how to solve a MIP with a start
%   guess.
%%

function mioinitsol()

    model = mosekmodel(name = "init sol", ...
                       numvar = 4, objective = [7 10 1 5], objsense = "max", ...
                       F = [1 1 1 1; eye(4)], ...
                       domain = [mosekdomain("eq", rhs = 2.5), ...
                                 mosekdomain("gt", dim = 4, rhs = zeros(4, 1))], ...
                       intvars = [1, 2, 3]);     % Specify the indices of integer variables

    % Specify start guess for the integer variables.
    model.setsolution("x", [1 1 0 nan]);

    % Request constructing the solution from integer variable values
    model.solve(param = ["MSK_IPAR_MIO_CONSTRUCT_SOL", "1"]);

    % Display the optimal solution.
    if model.hassolution("integer")
        x = model.getsolution("integer", "x");
        fprintf("x = [%s]\n", sprintf("%g ", x));

        % Was the initial solution used ?
        fprintf("Construct solution used?      %d\n", model.info.MSK_IINF_MIO_CONSTRUCT_SOLUTION)
        fprintf("Construct solution objective: %f\n", model.info.MSK_DINF_MIO_CONSTRUCT_SOLUTION_OBJ);
    end
end
