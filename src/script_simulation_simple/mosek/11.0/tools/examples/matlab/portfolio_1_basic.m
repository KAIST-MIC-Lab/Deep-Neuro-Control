%
%  File : portfolio_1_basic.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Purpose :   Implements a basic portfolio optimization model.
%
function [xx,prosta,solsta] = portfolio_1_basic()
    % Purpose:
    %     Computes the optimal portfolio for a given risk
    %
    % Input:
    %     n: Number of assets
    %     mu: An n dimmensional vector of expected returns
    %     GT: A matrix with n columns so (GT')*GT  = covariance matrix
    %     x0: Initial holdings
    %     w: Initial cash holding
    %     gamma: Maximum risk (=std. dev) accepted
    %
    % Output:
    %     Optimal expected return and the optimal portfolio

    n      = 8;
    w      = 59;
    mu     = [0.07197349 0.15518171 0.17535435 0.0898094  0.42895777 0.39291844 0.32170722 0.18378628]';
    x0     = [8.0 5.0 3.0 5.0 2.0 9.0 3.0 6.0]';
    gamma  = 36;
    GT     = [ 0.30758 0.12146 0.11341 0.11327 0.17625 0.11973 0.10435 0.10638; ...
               0.      0.25042 0.09946 0.09164 0.06692 0.08706 0.09173 0.08506; ...
               0.      0.      0.19914 0.05867 0.06453 0.07367 0.06468 0.01914; ...
               0.      0.      0.      0.20876 0.04933 0.03651 0.09381 0.07742; ...
               0.      0.      0.      0.      0.36096 0.12574 0.10157 0.0571 ; ...
               0.      0.      0.      0.      0.      0.21552 0.05663 0.06187; ...
               0.      0.      0.      0.      0.      0.      0.22514 0.03327; ...
               0.      0.      0.      0.      0.      0.      0.      0.2202 ];
    m = size(GT,1);

    model = mosekmodel(...
        name="Basic Markowitz", ...
        numvar = n, ...
        objsense = "maximize", ...
        objective = mu);

    % Bound on the variables (holdings). Shortselling is not allowed.
    model.appendcons(F = eye(n), domain = mosekdomain("nonnegative", n = n));
        
    % The amount invested  must be identical to initial wealth
    model.appendcons(F = ones(1, n), domain = mosekdomain("equal", rhs=w+sum(x0)));

    % Imposes a bound on the risk
    model.appendcons(F = [zeros(1,n); GT], g = [gamma; zeros(m,1)],...
                     domain = mosekdomain("qcone", dim = m+1));

    model.solve();

    % Check if the interior point solution is an optimal point
    % See https://docs.mosek.com/latest/matlabapi/accessing-solution.html about handling solution statuses.
    [ok,prosta,solsta] = model.hassolution("interior");
    if ~ ok || solsta ~= "OPTIMAL"
        error("Unexpected solution status");
    else
        xx = model.getsolution("interior","x");

        disp("\n-----------------------------------------------------------------------------------");
        disp("Basic Markowitz portfolio optimization");
        disp("-----------------------------------------------------------------------------------\n");
        fprintf("Expected return: %.4e Std. deviation: %.4e\n", mu' * xx, gamma);
    end

end
