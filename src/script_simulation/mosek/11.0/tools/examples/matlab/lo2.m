%%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File :      lo2.m
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
%%
function [xx,prosta,solsta] = lo2()
    model = mosekmodel(name = "lo2", ...
                       objsense = "maximize", ...
                       objective = [ 3 1 5 1 ]', ...
                       numvar = 4);

    model.appendcons(name="con-eq30", F = [ 3 1 2 0 ], domain = mosekdomain("equal",        rhs=30));
    model.appendcons(name="con-gt14", F = [ 2 1 3 1 ], domain = mosekdomain("greater than", rhs=15));
    model.appendcons(name="con-lt25", F = [ 0 2 0 3 ], domain = mosekdomain("less than",    rhs=25)); 
    model.appendcons(name="con-nneg", F = speye(4),    domain = mosekdomain("nonnegative", n=4));
    model.appendcons(name="con-lt10", F = [ 0 1 0 0 ], domain = mosekdomain("less than", rhs="10"));

    model.varname([1:4],["x","y","z","w"]);

    model.solve();
    if model.hassolution("interior")
        [xx,prosta,solsta] = model.getsolution("interior","x");
    end
end
