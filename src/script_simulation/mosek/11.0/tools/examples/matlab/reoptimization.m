%%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      reoptimization.m
%
%  Purpose:   Demonstrates how to solve a  linear
%             optimization problem using the Conic Toolbox for Matlab
%             and modify and re-optimize the problem.
%%
function reoptimization()
    model = moseklinmodel(...
        name      = "reoptimization",...
        objsense  = "maximize", ... 
        numvar    = 7, ...
        blx = [ 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]', ...
        bux = [ inf inf inf inf inf inf inf ]', ...
        b   = [ 100000.0 50000.0 60000.0 ]', ...
        ... %         x1  x2  x3  s1  s2  s3  s4
        objective = [ 1.5 2.5 3.0 0.0 0.0 0.0 0.0 ]', ...
        A         = [ 2.0 4.0 3.0 1.0 0.0 0.0 0.0 ; ...
                      3.0 2.0 3.0 0.0 1.0 0.0 0.0 ; ...
                      2.0 3.0 2.0 0.0 0.0 1.0 0.0 ]);
   
    model.solve();

    [result,prosta,solsta] = model.getsolution("x"); x = result';
    display(x);

    % Alter A and reoptimize
    fprintf("Update first row non-zeros\n");
    model.setrows([ 3.0 4.0 3.0 0.0 0.0 0.0 0.0 ], first=1);
    model.solve();

    [result,prosta,solsta] = model.getsolution("x"); x = result';
    display(x);

    % Add variable
    fprintf("Add variable\n");
    model.appendvars(1, ... 
                     bl = [ 0.0 ]', ...
                     A  = [ 4.0 0.0 1.0 ]', ...
                     c  = [ 1.0 ]);
    
    model.solve(param = [ "MSK_IPAR_OPTIMIZER","MSK_OPTIMIZER_FREE_SIMPLEX" ]);
    [result,prosta,solsta] = model.getsolution("x"); x = result';
    display(x);
   
    % Add constraint
    fprintf("Add constraints\n");
    %                  x1  x2  x3  s1  s2  s3  s4  x4
    model.appendcons([ 1.0 2.0 1.0 0.0 0.0 0.0 1.0 1.0 ], ...
                     [ 30000.0 ]);

    model.solve(param = [ "MSK_IPAR_OPTIMIZER","MSK_OPTIMIZER_FREE_SIMPLEX" ]);
    [result,prosta,solsta] = model.getsolution("x"); x = result';
    display(x);

    % Change bounds
    fprintf("Change bounds\n");
    model.setb([ 80000.0 40000.0 50000.0 22000.0 ]');
    
    model.solve(param = [ "MSK_IPAR_OPTIMIZER","MSK_OPTIMIZER_FREE_SIMPLEX" ]);
    [result,prosta,solsta] = model.getsolution("x"); x = result';
    display(x);
end
