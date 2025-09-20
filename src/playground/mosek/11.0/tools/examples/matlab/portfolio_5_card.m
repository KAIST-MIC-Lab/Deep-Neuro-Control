%%
%  File : portfolio_5_card.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model
%                 with cardinality constraints on number of assets traded.
% 
%%
function portfolio_5_card()

    n      = 8;
    w      = 1.0;
    mu     = [0.07197349 0.15518171 0.17535435 0.0898094 0.42895777 0.39291844 0.32170722 0.18378628]';
    x0     = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]';
    gamma  = 0.25;
    GT     = [  0.30758, 0.12146, 0.11341, 0.11327, 0.17625, 0.11973, 0.10435, 0.10638; ...
                0.     , 0.25042, 0.09946, 0.09164, 0.06692, 0.08706, 0.09173, 0.08506; ...
                0.     , 0.     , 0.19914, 0.05867, 0.06453, 0.07367, 0.06468, 0.01914; ...
                0.     , 0.     , 0.     , 0.20876, 0.04933, 0.03651, 0.09381, 0.07742; ...
                0.     , 0.     , 0.     , 0.     , 0.36096, 0.12574, 0.10157, 0.0571 ; ...
                0.     , 0.     , 0.     , 0.     , 0.     , 0.21552, 0.05663, 0.06187; ...
                0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.22514, 0.03327; ...
                0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.2202 ];


    %
    %    Description:
    %        Extends the basic Markowitz model with cardinality constraints.
    %
    %    Input:
    %        n: Number of assets
    %        mu: An n dimensional vector of expected returns
    %        GT: A matrix with n columns so (GT')*GT  = covariance matrix
    %        x0: Initial holdings 
    %        w: Initial cash holding
    %        gamma: Maximum risk (=std. dev) accepted
    %        k: Maximum number of assets on which we allow to change position.
    %
    %    Output:
    %       Optimal portfolio     
    %

    fprintf("\nMarkowitz portfolio optimization with cardinality constraints\n");

    for k=1:n
        % Upper bound on the traded amount
        u = w + sum(x0);

        % A model with variables (x, z, y)
        model = mosekmodel(name = "Transaction costs", ...
                           numvar = n+n+n, ...
                           objsense = "maximize", ...
                           objective = [mu; zeros(n,1); zeros(n,1)], ...
                           intvars = 2*n+1:3*n);        % The variables y are integer

        % Linear constraints
        % [ e'  0   0  ]           =   w + e'*x0
        % [ I  -I   0  ]   [ x ]  <=  x0
        % [ I   I   0  ] * [ z ]  >=  x0
        % [ 0   I  -U  ]   [ y ]  <=  0
        % [ 0   0   e' ]          <=  k

        % Linear budget constraint
        model.appendcons(name = "budget", F = [ ones(1,n), zeros(1, 2*n) ], ...
                         domain = mosekdomain("equal", rhs = u));

        % Linear absolute value model z >= |x-x0|
        model.appendcons(name = "abs1", F = [ speye(n) -speye(n) sparse(n,n) ], ...
                         domain = mosekdomain("lt", dim = n, rhs = x0));
        model.appendcons(name = "abs2", F = [ speye(n) speye(n) sparse(n,n) ], ...
                         domain = mosekdomain("gt", dim = n, rhs = x0));

        % "big-M" constraint z<=Uy
        model.appendcons(name = "bigM", F = [ sparse(n,n) speye(n) -u*speye(n) ], ...
                         domain = mosekdomain("lt", dim = n, rhs = zeros(n,1)));

        % Cardinality bound sum(y) <= k
        model.appendcons(name = "card", F = [ zeros(1, 2*n) ones(1, n) ], ...
                         domain = mosekdomain("lt", rhs = k));

        % No shortselling, x >= 0
        model.appendcons(F = speye(n, 3*n),  ...
                         domain = mosekdomain("greater than", dim = n, rhs = zeros(n,1)));

        % Affine conic constraint representing [ gamma, GT*x ] in quadratic cone
        model.appendcons(F = [sparse(1, 3*n); ...
                              GT, sparse(n, 2*n)], ...
                         g = [gamma; zeros(n, 1)], ...
                         domain = mosekdomain("qcone", dim = n + 1));

        model.solve(quiet = true);

        % Check if the solution is an optimal point
        % See https://docs.mosek.com/latest/matlabapi/accessing-solution.html about handling solution statuses.
        [ok,prosta,solsta] = model.hassolution("integer");
        if ~ ok || solsta ~= "INTEGER_OPTIMAL"
            error("Unexpected solution status");
        else
            xx = model.getsolution("integer","x");
            x = xx(1:n);
            fprintf("Bound: %d   Expected return: %.4e  Solution: [ %s ]\n", k, mu'*x, sprintf("%.2g ", x));
        end
    end
end
