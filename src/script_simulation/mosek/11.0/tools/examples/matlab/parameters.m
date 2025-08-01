%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      parameters.m
%
%  Purpose:    Demonstrates setting solver parameters and
%              retrieving information items.
%

function parameters()
    % Set up a small problem as an example
    model = mosekmodel(objsense = "min", ...
                       objective = [ 1, 2 ]', ...
                       numvar = 2, ...
                       F = [ 1 0 ; ...
                             0 1 ; ...
                             3 -1 ], ...
                       domain = [ mosekdomain("greater than", rhs = 2), ... 
                                  mosekdomain("greater than", rhs = 3), ...
                                  mosekdomain("less than", rhs = 1) ]);

    % Solve with a list of parameters
    model.solve(param = ["MSK_IPAR_LOG", "1", ...                           % Set log level (integer parameter)
                         "MSK_IPAR_CACHE_LICENSE", "MSK_OFF", ...           % Do not keep the license (integer parameter)
                         "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", "1.0e-7" ]);     % Set relative gap tolerance (double parameter)                      

    % Demonstrate information items after optimization
    fprintf("Optimizer time %.3f\n", model.info.MSK_DINF_OPTIMIZER_TIME);
    fprintf("#iterations    %d\n",   model.info.MSK_IINF_INTPNT_ITER);

end
