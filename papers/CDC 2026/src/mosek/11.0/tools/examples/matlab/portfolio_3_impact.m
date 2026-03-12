%%
%  File : portfolio_3_impact.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Description :  Implements a basic portfolio optimization model
%                 with transaction costs of order x^(3/2).
% 
%%

function [x,prosta,solsta] = portfolio_3_impact()

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
    %            m: It is assumed that  market impact cost for the j'th asset is
    %               m_j|x_j-x0_j|^3/2
    %
    %        Output:
    %           Optimal portfolio (x)
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

    % Somewhat arbirtrary choice of m
    m = 1.0e-2*ones(n,1);

    % A model with variables (x, t)
    model = mosekmodel(name = "Impact", ...
                       numvar = n+n, ...
                       objsense = "maximize", ...
                       objective = [mu; zeros(n,1)]);

    % Linear constraint
    % [ e' m' ]  * [ x; t ]  =   w + e'*x0
    model.appendcons(name = "budget", F = [ ones(1,n), m' ], ...
                     domain = mosekdomain("equal", rhs = w + sum(x0)));

    % No shortselling, x >= 0 and no other bounds
    model.appendcons(F = speye(n, 2*n),  ...
                     domain = mosekdomain("greater than", dim = n, rhs = zeros(n,1)));

    % Affine conic constraint representing [ gamma, GT*x ] in quadratic cone
    model.appendcons(F = [sparse(1, 2*n); ...
                          GT, sparse(n, n)], ...
                     g = [gamma; zeros(n, 1)], ...
                     domain = mosekdomain("qcone", dim = n + 1));


    % Power cone constraints representing t(i) >= |x(i)-x0(i)|^1.5
    fi = [];
    fj = [];
    g  = [];
    fv = repmat([1; 1], n, 1);
    for k=1:n
      fi = [fi; 3*k-2; 3*k];
      fj = [fj; n+k; k];
      g  = [g; 0; 1; -x0(k)];
    end

    model.appendcons(name = "impact", ...
                     F = sparse(fi, fj, fv, 3*n, 2*n), g = g, ...
                     domain = mosekdomain("power cone", dim = 3, n = n, alpha = [2.0, 1.0]));

    model.solve();

    % Check if the solution is an optimal point
    % See https://docs.mosek.com/latest/matlabapi/accessing-solution.html about handling solution statuses.
    [ok,prosta,solsta] = model.hassolution("interior");
    if ~ ok || solsta ~= "OPTIMAL"
        error("Unexpected solution status");
    else
        xx = model.getsolution("interior","x");
        x = xx(1:n);
    end
end
