%%
%  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File:      response.m
%
%  Purpose:   This example demonstrates the workflow for proper response handling
%

function response(args)
    % Ceate a simple model for the sake of the example
    model = mosekmodel(objsense = "min", ...
                       objective = [ 1, 2 ]', ...
                       numvar = 2, ...
                       F = [ 1 0 ; ...
                             0 1 ; ...
                             3 -1 ], ...
                       domain = [ mosekdomain("greater than", rhs = 2), ... 
                                  mosekdomain("greater than", rhs = 3), ...
                                  mosekdomain("less than", rhs = 1) ]);
    % Solve the model
    try
        model.solve();
    catch ME
        warning("An error during optimization; handle it here.");
        rethrow(ME);
    end
    
    % We check if the interior-point solution exists and what status it has
    [exists, prosta, solsta] = model.hassolution("interior");

    if exists
        disp("Solved the problem with statuses:");
        disp(prosta);
        disp(solsta);

        switch solsta
            case "OPTIMAL"
                disp("Optimal solution found:");
                x = model.getsolution("interior");
                disp(x);
            case "PRIM_INFEAS_CER"
                disp("The problem is primal infeasible.");
            case "DUAL_INFEAS_CER"
                disp("The problem is dual infeasible.");
            case "UNKNOWN"
                disp("Solution status UNKNOWN. This could indicate numerical issues");
            default
                disp("Another solution status:")
                disk(solsta)
        end
    else
        warning("Solution does not exists");
    end

end
