%%
%  File : pinfeas.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Purpose: Demonstrates how to fetch a primal infeasibility certificate
%           for a linear problem
%%

function pinfeas()

% In this example we set up a simple problem
[model] = testProblem();

% Perform the optimization.
% For comparison with the following code, we will also
% ask MOSEK to print the infeasibility certificate to the log
% following optimization.
model.solve(param = ["MSK_IPAR_INFEAS_REPORT_AUTO", "MSK_ON"]);

% Check problem status
[hassol, prosta, solsta] = model.hassolution("interior");
if hassol && prosta == "PRIM_INFEAS"
    % Set the tolerance at which we consider a dual value as essential
    eps = 1e-7;

    % Obtain the dual values (contianing certificate)
    y = model.getsolution("interior", "y");

    % List all certificate entries with (sufficiently) nonzero dual values
    disp("Constraint rows participating in infeasibility: ");
    cert = find(abs(y) > eps);
    disp(cert);
else 
    disp("The problem is not primal infeasible, no certificate to show");
end

% Set up a simple linear problem from the manual for test purposes
function [model] = testProblem()
    model = mosekmodel(name = "pinfeas", numvar = 7, ...
                       objsense = "minimize", objective = [1 2 5 2 1 2 1]');
    model.appendcons(name = "s", ...
                     F = sparse([1,1,2,2,3,3,3], [1,2,3,4,5,6,7], [1,1,1,1,1,1,1]), ...
                     domain = mosekdomain("less than", dim = 3, rhs=[200, 1000, 1000]));
    model.appendcons(name = "d", ...
                     F = sparse([1,1,2,3,3,4,4], [1,5,2,3,6,4,7], [1,1,1,1,1,1,1]), ...
                     domain = mosekdomain("equals", dim = 4, rhs=[1100, 200, 500, 500]));
    model.appendcons(name = "low", F = speye(7), domain = mosekdomain("rplus", dim = 7));

