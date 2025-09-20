%%
%  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  File:      helloworld.m
%
%  The most basic example of how to get started with MOSEK.
%

function [xx,prosta,solsta] = helloworld()
    model = mosekmodel(...
                  name = "helloworld", ...
                  objsense = "min", ...
                  objective = [ 1, 2 ]', ...
                  numvar = 2, ...
                  F = [ 1 0 ; ...
                        0 1 ; ...
                        3 -1 ], ...
                  domain = [ mosekdomain("greater than", rhs = 2), ...   % x >= 2,    1st row of F
                             mosekdomain("greater than", rhs = 3), ...   % y >= 3,    2nd row of F
                             mosekdomain("less than", rhs = 1) ]);       % 3x-y <= 1, 3rd row of F
    model.solve();
    
    if model.hassolution("interior")
        [xx,prosta,solsta] = model.getsolution("any", "x");
        xx
    end

end
