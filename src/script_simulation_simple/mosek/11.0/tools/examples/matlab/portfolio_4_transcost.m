%%
%  File : portfolio_4_transcost.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model
%                 with fixed setup costs and transaction costs
%                 as a mixed-integer problem.
% 
%%

function [x,prosta,solsta] = portfolio_4_transcost()

    %
    %        Description:
    %            Extends the basic Markowitz model with a market cost term.
    %
    %        Input:
    %            n: Number of assets
    %            mu: An n dimmensional vector of expected returns
    %            GT: A matrix with n columns so (GT')*GT  = covariance matrix
    %            x0: Initial holdings 
    %            w: Initial cash holding
    %            gamma: Maximum risk (=std. dev) accepted
    %            f: If asset j is traded then a fixed cost f_j must be paid
    %            g: If asset j is traded then a cost g_j must be paid for each unit traded
    %
    %        Output:
    %           Optimal portfolio     
    %

    n      = 8;
    w      = 1.0;
    mu     = [0.07197349 0.15518171 0.17535435 0.0898094 0.42895777 0.39291844 0.32170722 0.18378628]';
    x0     = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]';
    gamma  = 0.36;
    GT     = [  0.30758, 0.12146, 0.11341, 0.11327, 0.17625, 0.11973, 0.10435, 0.10638; ...
                0.     , 0.25042, 0.09946, 0.09164, 0.06692, 0.08706, 0.09173, 0.08506; ...
                0.     , 0.     , 0.19914, 0.05867, 0.06453, 0.07367, 0.06468, 0.01914; ...
                0.     , 0.     , 0.     , 0.20876, 0.04933, 0.03651, 0.09381, 0.07742; ...
                0.     , 0.     , 0.     , 0.     , 0.36096, 0.12574, 0.10157, 0.0571 ; ...
                0.     , 0.     , 0.     , 0.     , 0.     , 0.21552, 0.05663, 0.06187; ...
                0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.22514, 0.03327; ...
                0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.2202 ];

    f = 0.01*ones(n,1);
    g = 0.001*ones(n,1);

    % Upper bound on the traded amount
    u = w + sum(x0);

    % A model with variables (x, z, y)
    model = mosekmodel(name = "Impact", ...
                       numvar = n+n+n, ...
                       objsense = "maximize", ...
                       objective = [mu; zeros(n,1); zeros(n,1)], ...
                       intvars = 2*n+1:3*n);        % The variables y are integer

    % y is binary, ie. in [0,1]
    model.appendcons(F = [sparse(n,2*n) speye(n,n)], domain = mosekdomain("gt", dim = n, rhs = zeros(n,1)));
    model.appendcons(F = [sparse(n,2*n) speye(n,n)], domain = mosekdomain("lt", dim = n, rhs = ones(n,1)));

    % Linear constraints
    % [ e'  g'  f' ]   [ x ]  =   w + e'*x0
    % [ I  -I   0  ] * [ z ]  <=  x0
    % [ I   I   0  ]   [ y ]  >=  x0
    % [ 0   I  -U  ]          <=  0

    % Linear budget constraint
    model.appendcons(name = "budget", F = [ ones(1,n), g', f' ], ...
                     domain = mosekdomain("equal", rhs = u));

    % Linear absolute value model z >= |x-x0|
    model.appendcons(name = "abs1", F = [ speye(n) -speye(n) sparse(n,n) ], ...
                     domain = mosekdomain("lt", dim = n, rhs = x0));
    model.appendcons(name = "abs2", F = [ speye(n) speye(n) sparse(n,n) ], ...
                     domain = mosekdomain("gt", dim = n, rhs = x0));

    % "big-M" constraint z<=Uy
    model.appendcons(name = "bigM", F = [ sparse(n,n) speye(n) -u*speye(n) ], ...
                     domain = mosekdomain("lt", dim = n, rhs = zeros(n,1)));

    % No shortselling, x >= 0
    model.appendcons(F = speye(n, 3*n),  ...
                     domain = mosekdomain("greater than", dim = n, rhs = zeros(n,1)));

    % Affine conic constraint representing [ gamma, GT*x ] in quadratic cone
    model.appendcons(F = [sparse(1, 3*n); ...
                          GT, sparse(n, 2*n)], ...
                     g = [gamma; zeros(n, 1)], ...
                     domain = mosekdomain("qcone", dim = n + 1));

    model.solve();

    % Check if the solution is an optimal point
    % See https://docs.mosek.com/latest/matlabapi/accessing-solution.html about handling solution statuses.
    [ok,prosta,solsta] = model.hassolution("integer");
    if ~ ok || solsta ~= "INTEGER_OPTIMAL"
        error("Unexpected solution status");
    else
        xx = model.getsolution("integer","x");
        x = xx(1:n);
        y = xx(2*n+(1:n));
        z = xx(n+(1:n));
    end

    fprintf("\nMarkowitz portfolio optimization with transactions cost")
    fprintf("Expected return: %.4e Std. deviation: %.4e Transactions cost: %.4e", ...
                     mu'*x, gamma, f'*y+g'*z)
end
