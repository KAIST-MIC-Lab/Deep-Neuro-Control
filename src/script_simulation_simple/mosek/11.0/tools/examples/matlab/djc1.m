%%
%   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%   File:      djc1.m
%
%   Purpose: Demonstrates how to solve the problem with two disjunctions:
%
%      minimize    2x0 + x1 + 3x2 + x3
%      subject to   x0 + x1 + x2 + x3 >= -10
%                  (x0-2x1<=-1 and x2=x3=0) or (x2-3x3<=-2 and x1=x2=0)
%                  x0=2.5 or x1=2.5 or x2=2.5 or x3=2.5
%%

function [xx,prosta,solsta] = djc1()
    model = mosekmodel(name = "djc1", ...
                       objsense = "minimize", ...
                       objective = [ 2,1,3,1 ]', ...
                       numvar = 4);

    % A disjunction of two clauses
    model.appenddisjunction( [ model.clause(F = [1 -2 0 0 ;
                                                 0  0 1 0 ; 
                                                 0  0 0 1], ...
                                            domain=[ mosekdomain("less than", rhs = [-1]), ...                % 1st simple term of 1st clause
                                                     mosekdomain("equal",     dim = 2, rhs = [0 0]') ]), ...  % 2nd simple term of 1st clause

                               model.clause(F = [0 0 1 -3 ; 
                                                 1 0 0  0 ; 
                                                 0 1 0  0], ...
                                            domain = [ mosekdomain("less than", rhs = -2), ...                % 1st simple term of 2nd clause
                                                       mosekdomain("equal", dim = 2, rhs = [0 0]') ]) ], ...  % 2nd simple term of 2nd clause
                            name = "first-djc");

    % A disjunction of four clauses
    model.appenddisjunction([ model.clause(F = [1 0 0 0], domain = mosekdomain("equal", rhs = 2.5)),...
                              model.clause(F = [0 1 0 0], domain = mosekdomain("equal", rhs = 2.5)),...
                              model.clause(F = [0 0 1 0], domain = mosekdomain("equal", rhs = 2.5)),...
                              model.clause(F = [0 0 0 1], domain = mosekdomain("equal", rhs = 2.5)) ],...
                            name = "second-djc");

    % The standard liner constraint
    model.appendcons(name = "C", F = [1 1 1 1], domain = mosekdomain("greater than", rhs = -10));

    model.solve();

    if model.hassolution("integer")
        [xx,prosta,solsta] = model.getsolution("integer","x");
        fprintf("Solution : [%s]\n", sprintf("%g ", xx));
    end
end

