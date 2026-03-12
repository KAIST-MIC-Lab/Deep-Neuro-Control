%
%  File : portfolio_2_frontier.m
%
%  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
%
%  Purpose :   Implements a basic portfolio optimization model.
%              Computes points on the efficient frontier.
%
function [xx,prosta,solsta] = portfolio_2_frontier()
    %  Purpose:
    %      Computes several portfolios on the optimal portfolios by

    %          for alpha in alphas:
    %              maximize   expected return - alpha * variance
    %              subject to the constraints

    %  Input:
    %      n: Number of assets
    %      mu: An n dimmensional vector of expected returns
    %      GT: A matrix with n columns so (GT')*GT  = covariance matrix
    %      x0: Initial holdings
    %      w: Initial cash holding
    %      alphas: List of the alphas

    %  Output:
    %      The efficient frontier as list of tuples (alpha, expected return, variance)
    
    n      = 8;
    w      = 1.0;
    mu     = [0.07197, 0.15518, 0.17535, 0.08981, 0.42896, 0.39292, 0.32171, 0.18379]';
    x0     = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    GT     = [ 0.30758, 0.12146, 0.11341, 0.11327, 0.17625, 0.11973, 0.10435, 0.10638;
               0.     , 0.25042, 0.09946, 0.09164, 0.06692, 0.08706, 0.09173, 0.08506;
               0.     , 0.     , 0.19914, 0.05867, 0.06453, 0.07367, 0.06468, 0.01914;
               0.     , 0.     , 0.     , 0.20876, 0.04933, 0.03651, 0.09381, 0.07742;
               0.     , 0.     , 0.     , 0.     , 0.36096, 0.12574, 0.10157, 0.0571 ;
               0.     , 0.     , 0.     , 0.     , 0.     , 0.21552, 0.05663, 0.06187;
               0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.22514, 0.03327;
               0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.2202 ]; 

    % Some predefined alphas are chosen
    alphas = [ 0.0, 0.01, 0.1, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 10.0 ]';
    niter = size(alphas,1);
    frontier_mux = zeros(niter,1);
    frontier_s   = zeros(niter,1);

    for i = 1:niter
        alpha = alphas(i);
        model = mosekmodel(name = "Efficient frontier", numvar = n+1);

        % Defines the variables (holdings). Shortselling is not allowed.
        model.appendcons(F = speye(n), domain = mosekdomain("greater than", dim = n, rhs = zeros(n,1)));

        x_idxs = [1:n];
        s_idxs = n+1;

        model.appendcons(name = "budget", F = [ones(1,n) 0], domain = mosekdomain("equal", rhs = w + sum(x0)));

        % Computes the risk
        model.appendcons(name="variance", ...
                         F=[ sparse(1,n+1,[1.0],[1],[s_idxs]) ; ...
                             sparse(1,n+1) ; ...
                             GT zeros(n,1)], ...
                         g=[0; 0.5; zeros(n,1)], ...
                         domain = mosekdomain("rqcone", dim = n+2));

        %  Define objective as a weighted combination of return and variance
        model.objective("maximize", [ mu ; -alpha ]);

        % Solve the problem for the current alpha
        model.solve();

        % Check if the solution is an optimal point
        % See https://docs.mosek.com/latest/matlabapi/accessing-solution.html about handling solution statuses.
        [ok,prosta,solsta] = model.hassolution("interior");
        if ~ ok || solsta ~= "OPTIMAL"
            error("Unexpected solution status");
        else
            xx = model.getsolution("interior","x");
            x = xx(x_idxs,1);
            s = xx(s_idxs,1);

            frontier_mux(i,1) = mu' * x;
            frontier_s(i,1)   = s;
        end
    end

    fprintf("\n-----------------------------------------------------------------------------------\n");
    fprintf("Efficient frontier") ;
    fprintf("\n-----------------------------------------------------------------------------------\n");
    fprintf("\t%-12s  %-12s  %-12s\n", "alpha", "return", "std. dev.");

    sqrt_frontier_s = sqrt(frontier_s);
    for i = 1:niter
        fprintf("\t%-12.4f  %-12.4e  %-12.4e\n", alphas(i,1), frontier_mux(i,1), frontier_s(i,1));
    end

end
